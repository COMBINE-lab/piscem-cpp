#pragma once

#include "../include/parallel_hashmap/phmap.h"
#include "../include/reference_index.hpp"
#include "../include/CanonicalKmerIterator.hpp"
#include "../include/query/streaming_query_canonical_parsing.hpp"
#include "../include/projected_hits.hpp"
#include "../include/FastxParser.hpp"

#include "../include/hit_searcher.hpp"

#include <vector>
#include <cassert>
#include <fstream>
#include <cmath>  // for std::ceil on linux
#include <numeric>

namespace mapping {

namespace util {

constexpr int32_t invalid_frag_len = std::numeric_limits<int32_t>::min();
constexpr int32_t invalid_mate_pos = std::numeric_limits<int32_t>::min();

struct simple_hit {
    bool is_fw{false};
    bool mate_is_fw{false};
    int32_t pos{-1};
    float score{0.0};
    uint32_t num_hits{0};
    uint32_t tid{std::numeric_limits<uint32_t>::max()};
    int32_t mate_pos{std::numeric_limits<int32_t>::max()};
    int32_t fragment_length{std::numeric_limits<int32_t>::max()};
    inline bool valid_pos(int32_t read_len, uint32_t txp_len, int32_t max_over) {
        int32_t signed_txp_len = static_cast<int32_t>(txp_len);
        return (pos > -max_over) and ((pos + read_len) < (signed_txp_len + max_over));
    }
    inline bool has_mate() const { return mate_pos != invalid_mate_pos; }
    inline bool mate_is_mapped() const { return mate_pos != invalid_mate_pos; }
    inline int32_t frag_len() const { return (fragment_length != invalid_frag_len) ? fragment_length : 0; }
};

enum class MappingType : uint8_t { UNMAPPED, SINGLE_MAPPED , MAPPED_FIRST_ORPHAN, MAPPED_SECOND_ORPHAN, MAPPED_PAIR };

enum class HitDirection : uint8_t { FW, RC, BOTH };

struct sketch_hit_info {
    // add a hit to the current target that occurs in the forward
    // orientation with respect to the target.
    bool add_fw(int32_t ref_pos, int32_t read_pos, int32_t rl, int32_t k, int32_t max_stretch,
                float score_inc) {
        (void)rl;
        (void)k;
        bool added{false};

        // since hits are collected by moving _forward_ in the
        // read, if this is a fw hit, it should be moving
        // forward in the reference. Only add it if this is
        // the case.  This ensure that we don't
        // double-count a k-mer that might occur twice on
        // this target.
        if (ref_pos > last_ref_pos_fw and read_pos > last_read_pos_fw) {
            if (last_read_pos_fw == -1) {
                approx_pos_fw = ref_pos - read_pos;
            } else {
                if ((ref_pos - approx_pos_fw) > max_stretch) { return false; }
            }
            // if (last_ref_pos_fw > -1 and (ref_pos > last_ref_pos_fw + 15)) { return false; }
            last_ref_pos_fw = ref_pos;
            last_read_pos_fw = read_pos;
            fw_score += score_inc;
            ++fw_hits;
            added = true;
        }
        return added;
    }

    // add a hit to the current target that occurs in the forward
    // orientation with respect to the target.
    bool add_rc(int32_t ref_pos, int32_t read_pos, int32_t rl, int32_t k, int32_t max_stretch,
                float score_inc) {
        bool added{false};
        // since hits are collected by moving _forward_ in the
        // read, if this is an rc hit, it should be moving
        // backwards in the reference. Only add it if this is
        // the case.
        // This ensures that we don't double-count a k-mer that
        // might occur twice on this target.
        if (ref_pos < last_ref_pos_rc and read_pos > last_read_pos_rc) {
            approx_pos_rc = (ref_pos - (rl - (read_pos + k)));
            if (last_read_pos_rc == -1) {
                approx_end_pos_rc = ref_pos + read_pos;
            } else {
                if (approx_end_pos_rc - approx_pos_rc > max_stretch) { return false; }
            }
            // if (last_ref_pos_rc > -1 and ref_pos < last_ref_pos_rc - 15) { return false; }
            last_ref_pos_rc = ref_pos;
            last_read_pos_rc = read_pos;
            rc_score += score_inc;
            ++rc_hits;
            added = true;
        }
        return added;
    }

    inline uint32_t max_hits_for_target() { return std::max(fw_hits, rc_hits); }

    // true if forward, false if rc
    // second element is score
    inline HitDirection best_hit_direction() {
        int32_t fw_minus_rc = static_cast<int32_t>(fw_hits) - static_cast<int32_t>(rc_hits);
        return (fw_minus_rc > 0) ? HitDirection::FW
                                 : ((fw_minus_rc < 0) ? HitDirection::RC : HitDirection::BOTH);
    }

    inline simple_hit get_fw_hit() {
        return simple_hit{true,     false,   approx_pos_fw,
                          fw_score, fw_hits, std::numeric_limits<uint32_t>::max()};
    }

    inline simple_hit get_rc_hit() {
        return simple_hit{false,    false,   approx_pos_rc,
                          rc_score, rc_hits, std::numeric_limits<uint32_t>::max()};
    }

    inline simple_hit get_best_hit() {
        auto best_direction = best_hit_direction();
        return (best_direction != HitDirection::RC)
                   ? simple_hit{true,     false,   approx_pos_fw,
                                fw_score, fw_hits, std::numeric_limits<uint32_t>::max()}
                   : simple_hit{false,    false,   approx_pos_rc,
                                rc_score, rc_hits, std::numeric_limits<uint32_t>::max()};
    }

    inline std::string to_string() {
        std::stringstream ss;
        ss << "fw_hits: " << fw_hits << ", fw_score : " << fw_score
           << ", fw_pos : " << approx_pos_fw << " || rc_hits: " << rc_hits
           << ", rc_score: " << rc_score << ", rc_pos: " << approx_pos_rc;
        return ss.str();
    }

    int32_t last_read_pos_fw{-1};
    int32_t last_read_pos_rc{-1};

    int32_t last_ref_pos_fw{-1};
    int32_t last_ref_pos_rc{std::numeric_limits<int32_t>::max()};

    int32_t approx_pos_fw{-1};
    int32_t approx_pos_rc{-1};
    int32_t approx_end_pos_rc{-1};

    uint32_t fw_hits{0};
    uint32_t rc_hits{0};
    float fw_score{0.0};
    float rc_score{0.0};
};

struct mapping_cache_info {
public:
    mapping_cache_info(mindex::reference_index& ri) : k(ri.k()), q(ri.get_dict()), hs(&ri) {}

    inline void clear() {
        map_type = mapping::util::MappingType::UNMAPPED;
        q.start();
        hs.clear();
        hit_map.clear();
        accepted_hits.clear();
        has_matching_kmers = false;
    }

    // will store how the read mapped
    mapping::util::MappingType map_type{mapping::util::MappingType::UNMAPPED};

    // map from reference id to hit info
    phmap::flat_hash_map<uint32_t, mapping::util::sketch_hit_info> hit_map;
    std::vector<mapping::util::simple_hit> accepted_hits;

    // map to recall the number of unmapped reads we see
    // for each barcode
    phmap::flat_hash_map<uint64_t, uint32_t> unmapped_bc_map;

    size_t max_occ_default = 200;
    size_t max_occ_recover = 1000;
    const bool attempt_occ_recover = (max_occ_recover > max_occ_default);
    size_t alt_max_occ = 2500;
    size_t k{0};

    // to perform queries
    sshash::streaming_query_canonical_parsing q;
    // implements the PASC algorithm
    mindex::hit_searcher hs;
    size_t max_chunk_reads = 5000;
    // regardless of having full mappings, did any k-mers match
    bool has_matching_kmers{false};
};

inline bool map_read(std::string* read_seq, mapping_cache_info& map_cache) {
    map_cache.clear();
    // rebind map_cache variables to
    // local names
    auto& q = map_cache.q;
    auto& hs = map_cache.hs;
    auto& hit_map = map_cache.hit_map;
    auto& accepted_hits = map_cache.accepted_hits;
    auto& map_type = map_cache.map_type;
    const bool attempt_occ_recover = map_cache.attempt_occ_recover;
    auto k = map_cache.k;

    map_cache.has_matching_kmers = hs.get_raw_hits_sketch(*read_seq, q, true, false);
    bool early_stop = false;

    // if there were hits
    if (map_cache.has_matching_kmers) {
        uint32_t num_valid_hits{0};
        uint64_t total_occs{0};
        uint64_t largest_occ{0};
        auto& raw_hits = hs.get_left_hits();

        // SANITY
        decltype(raw_hits[0].first) prev_read_pos = -1;
        // the maximum span the supporting k-mers of a
        // mapping position are allowed to have.
        // NOTE this is still > read_length b/c the stretch is measured wrt the
        // START of the terminal k-mer.
        int32_t max_stretch = static_cast<int32_t>(read_seq->length() * 1.0);

        // a raw hit is a pair of read_pos and a projected hit

        // the least frequent hit for this fragment.
        uint64_t min_occ = std::numeric_limits<uint64_t>::max();

        // this is false by default and will be set to true
        // if *every* collected hit for this fragment occurs
        // max_occ_default times or more.
        bool had_alt_max_occ = false;
        int32_t signed_rl = static_cast<int32_t>(read_seq->length());
        auto collect_mappings_from_hits =
            [&max_stretch, &min_occ, &hit_map, &num_valid_hits, &total_occs, &largest_occ,
             &early_stop, signed_rl, k](auto& raw_hits, auto& prev_read_pos, auto& max_allowed_occ,
                                        auto& had_alt_max_occ) -> bool {
            for (auto& raw_hit : raw_hits) {
                auto& read_pos = raw_hit.first;
                auto& proj_hits = raw_hit.second;
                auto& refs = proj_hits.refRange;

                uint64_t num_occ = static_cast<uint64_t>(refs.size());
                min_occ = std::min(min_occ, num_occ);
                had_alt_max_occ = true;

                bool still_have_valid_target = false;
                prev_read_pos = read_pos;
                if (num_occ <= max_allowed_occ) {
                    total_occs += num_occ;
                    largest_occ = (num_occ > largest_occ) ? num_occ : largest_occ;
                    float score_inc = 1.0;

                    for (auto v : refs) {
                        // uint64_t v = *pos_it;
                        const auto& ref_pos_ori = proj_hits.decode_hit(v);
                        uint32_t tid = sshash::util::transcript_id(v);
                        int32_t pos = static_cast<int32_t>(ref_pos_ori.pos);
                        bool ori = ref_pos_ori.isFW;
                        auto& target = hit_map[tid];

                        // Why >= here instead of == ?
                        // Because hits can happen on the same target in both the forward
                        // and rc orientations, it is possible that we start the loop with
                        // the target having num_valid_hits hits in a given orientation (o)
                        // we see a new hit for this target in oriention o (now it has
                        // num_valid_hits + 1) then we see a hit for this target in
                        // orientation rc(o).  We still want to add / consider this hit, but
                        // max_hits_for_target() > num_valid_hits. So, we must allow for
                        // that here.
                        if (target.max_hits_for_target() >= num_valid_hits) {
                            if (ori) {
                                target.add_fw(pos, static_cast<int32_t>(read_pos), signed_rl, k,
                                              max_stretch, score_inc);
                            } else {
                                target.add_rc(pos, static_cast<int32_t>(read_pos), signed_rl, k,
                                              max_stretch, score_inc);
                            }

                            still_have_valid_target |=
                                (target.max_hits_for_target() >= num_valid_hits + 1);
                        }

                    }  // DONE: for (auto &pos_it : refs)

                    ++num_valid_hits;

                    // if there are no targets reaching the valid hit threshold, then break
                    // early
                    if (!still_have_valid_target) { return true; }

                }  // DONE : if (num_occ <= max_allowed_occ)
            }      // DONE : for (auto& raw_hit : raw_hits)

            return false;
        };

        bool _discard = false;
        auto mao_first_pass = map_cache.max_occ_default - 1;
        early_stop = collect_mappings_from_hits(raw_hits, prev_read_pos, mao_first_pass, _discard);

        // If our default threshold was too stringent, then fallback to a more liberal
        // threshold and look up the k-mers that occur the least frequently.
        // Specifically, if the min occuring hits have frequency < max_occ_recover (2500 by
        // default) times, then collect the min occuring hits to get the mapping.
        if (attempt_occ_recover and (min_occ >= map_cache.max_occ_default) and
            (min_occ < map_cache.max_occ_recover)) {
            prev_read_pos = -1;
            uint64_t max_allowed_occ = min_occ;
            early_stop = collect_mappings_from_hits(raw_hits, prev_read_pos, max_allowed_occ,
                                                    had_alt_max_occ);
        }

        uint32_t best_alt_hits = 0;
        // int32_t signed_read_len = static_cast<int32_t>(record.seq.length());

        for (auto& kv : hit_map) {
            auto best_hit_dir = kv.second.best_hit_direction();
            // if the best direction is FW or BOTH, add the fw hit
            // otherwise add the RC.
            auto simple_hit = (best_hit_dir != mapping::util::HitDirection::RC)
                                  ? kv.second.get_fw_hit()
                                  : kv.second.get_rc_hit();

            if (simple_hit.num_hits >= num_valid_hits) {
                simple_hit.tid = kv.first;
                accepted_hits.emplace_back(simple_hit);
                // if we had equally good hits in both directions
                // add the rc hit here (since we added the fw)
                // above if the best hit was either FW or BOTH
                if (best_hit_dir == mapping::util::HitDirection::BOTH) {
                    auto second_hit = kv.second.get_rc_hit();
                    second_hit.tid = kv.first;
                    accepted_hits.emplace_back(second_hit);
                }
            } else {
                // best_alt_score = simple_hit.score > best_alt_score ? simple_hit.score :
                // best_alt_score;
                best_alt_hits =
                    simple_hit.num_hits > best_alt_hits ? simple_hit.num_hits : best_alt_hits;
            }
        }

        // alt_max_occ = had_alt_max_occ ? accepted_hits.size() : max_occ_default;

        /*
         * This rule; if enabled, allows through mappings missing a single hit, if there
         * was no mapping with all hits. NOTE: this won't work with the current early-exit
         * optimization however.
        if (accepted_hits.empty() and (num_valid_hits > 1) and (best_alt_hits >=
        num_valid_hits
        - 1)) { for (auto& kv : hit_map) { auto simple_hit = kv.second.get_best_hit(); if
        (simple_hit.num_hits >= best_alt_hits) {
              //if (simple_hit.valid_pos(signed_read_len, transcripts[kv.first].RefLength,
        10)) { simple_hit.tid = kv.first; accepted_hits.emplace_back(simple_hit);
              //}
            }
          }
        }
        */
    }  // DONE : if (rh)

    // If the read mapped to > maxReadOccs places, discard it
    if (accepted_hits.size() > map_cache.alt_max_occ) {
        accepted_hits.clear();
        map_type = mapping::util::MappingType::UNMAPPED;
    } else if (!accepted_hits.empty()) {
        map_type = mapping::util::MappingType::SINGLE_MAPPED;
    }

    return early_stop;
}

inline void merge_se_mappings(mapping_cache_info& map_cache_left,
                              mapping_cache_info& map_cache_right, int32_t left_len,
                              int32_t right_len, mapping_cache_info& map_cache_out) {
    map_cache_out.clear();
    auto& accepted_left = map_cache_left.accepted_hits;
    auto& accepted_right = map_cache_right.accepted_hits;

    size_t had_matching_kmers_left = map_cache_left.has_matching_kmers;
    size_t had_matching_kmers_right = map_cache_right.has_matching_kmers;

    size_t num_accepted_left = accepted_left.size();
    size_t num_accepted_right = accepted_right.size();

    if ((num_accepted_left > 0) and (num_accepted_right > 0)) {
        // look for paired end mappings
        // so we have to sort our accepted hits
        struct {
            // sort first by orientation, then by transcript id, and finally by position
            bool operator()(const mapping::util::simple_hit& a,
                            const mapping::util::simple_hit& b) {
                if (a.is_fw != b.is_fw) { return a.is_fw > b.is_fw; }
                // orientations are the same
                if (a.tid != b.tid) { return a.tid < b.tid; }
                return a.pos < b.pos;
            }
        } simple_hit_less;
        std::sort(accepted_left.begin(), accepted_left.end(), simple_hit_less);
        std::sort(accepted_right.begin(), accepted_right.end(), simple_hit_less);

        const mapping::util::simple_hit smallest_rc_hit = {false, false, -1, 0.0, 0, 0};
        // start of forward sub-list
        auto first_fw1 = accepted_left.begin();
        // end of forward sub-list is first non-forward hit
        auto last_fw1 = std::lower_bound(accepted_left.begin(), accepted_left.end(),
                                         smallest_rc_hit, simple_hit_less);
        // start of rc list
        auto first_rc1 = last_fw1;
        // end of rc list
        auto last_rc1 = accepted_left.end();

        // start of forward sub-list
        auto first_fw2 = accepted_right.begin();
        // end of forward sub-list is first non-forward hit
        auto last_fw2 = std::lower_bound(accepted_right.begin(), accepted_right.end(),
                                         smallest_rc_hit, simple_hit_less);
        // start of rc list
        auto first_rc2 = last_fw2;
        // end of rc list
        auto last_rc2 = accepted_right.end();

        auto back_inserter = std::back_inserter(map_cache_out.accepted_hits);
        using iter_t = decltype(first_fw1);
        using out_iter_t = decltype(back_inserter);

        auto merge_lists = [left_len, right_len](iter_t first1, iter_t last1, iter_t first2, iter_t last2,
                              out_iter_t out) -> out_iter_t {
            // https://en.cppreference.com/w/cpp/algorithm/set_intersection
            while (first1 != last1 && first2 != last2) {
                if (first1->tid < first2->tid) {
                    ++first1;
                } else {
                    if (!(first2->tid < first1->tid)) {
                        // first1->tid == first2->tid have the same transcript.
                        int32_t pos_fw = first1->is_fw ? first1->pos : first2->pos;
                        int32_t pos_rc = first1->is_fw ? first2->pos : first1->pos;
                        int32_t frag_len = (pos_rc - pos_fw);
                        if ((-20 < frag_len) and (frag_len < 1000)) {
                            // if left is fw and right is rc then
                            // fragment length is (right_pos + right_len - left_pos) + 1
                            // otherwise it is (left_pos + left_len - right_pos) + 1
                            bool right_is_rc = !first2->is_fw;
                            int32_t tlen =
                                right_is_rc ? ((first2->pos + right_len - first1->pos) + 1)
                                            : ((first1->pos + left_len - first2->pos) + 1);
                            *out++ = {first1->is_fw, first2->is_fw, first1->pos, 0.0, 0, first1->tid, first2->pos, tlen};
                            ++first1;
                        }
                    }
                    ++first2;
                }
            }
            return out;
        };

        // find hits of form 1:fw, 2:rc
        merge_lists(first_fw1, last_fw1, first_rc2, last_rc2, back_inserter);
        // find hits of form 1:rc, 2:fw
        merge_lists(first_rc1, last_rc1, first_fw2, last_fw2, back_inserter);

        map_cache_out.map_type = (map_cache_out.accepted_hits.size() > 0) ?
          MappingType::MAPPED_PAIR : MappingType::UNMAPPED;
    } else if ((num_accepted_left > 0) and !had_matching_kmers_right) {
        // just return the left mappings
        std::swap(map_cache_left.accepted_hits, map_cache_out.accepted_hits);
        map_cache_out.map_type = (map_cache_out.accepted_hits.size() > 0) ?
          MappingType::MAPPED_FIRST_ORPHAN : MappingType::UNMAPPED;
    } else if ((num_accepted_right > 0) and !had_matching_kmers_left) {
        // just return the right mappings
        std::swap(map_cache_right.accepted_hits, map_cache_out.accepted_hits);
        map_cache_out.map_type = (map_cache_out.accepted_hits.size() > 0) ?
          MappingType::MAPPED_SECOND_ORPHAN : MappingType::UNMAPPED;
    } else {
        // return nothing
    }
}

}  // namespace util

}  // namespace mapping
