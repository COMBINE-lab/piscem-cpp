#pragma once

#include "../include/parallel_hashmap/phmap.h"
#include "../include/reference_index.hpp"
#include "../include/CanonicalKmerIterator.hpp"
#include "../include/query/streaming_query_canonical_parsing.hpp"
#include "../include/projected_hits.hpp"
#include "../include/FastxParser.hpp"
#include "../include/itlib/small_vector.hpp"

#include "../include/hit_searcher.hpp"

#include <algorithm>
#include <limits>
#include <vector>
#include <cassert>
#include <fstream>
#include <cmath>  // for std::ceil on linux
#include <numeric>
#include <type_traits>

namespace mapping {

namespace util_bin {

constexpr int32_t invalid_frag_len = std::numeric_limits<int32_t>::min();
constexpr int32_t invalid_mate_pos = std::numeric_limits<int32_t>::min();

enum class orientation_filter : uint8_t { NONE, FORWARD_ONLY, RC_ONLY };

inline std::pair<uint32_t, uint32_t> get_bin_id(int32_t pos, int32_t bin_size=20000, int32_t overlap=1000) {
    uint32_t bin1 = (pos+1)/bin_size; // 1 added since 0 based
    uint32_t bin2 = (pos+1) > (bin1+1)*bin_size-overlap ? (bin1+1) :
        std::numeric_limits<uint32_t>::max(); // std::numeric_limits<uint32_t>::max() indicates that the kmer does not belong to the overlapping region
    return {bin1, bin2};
}
struct simple_hit {
    bool is_fw{false};
    bool mate_is_fw{false};
    int32_t pos{-1};
    float score{0.0};
    uint32_t num_hits{0};
    uint32_t tid{std::numeric_limits<uint32_t>::max()};
    uint32_t bin_id{std::numeric_limits<uint32_t>::max()};
    int32_t mate_pos{std::numeric_limits<int32_t>::max()};
    int32_t fragment_length{std::numeric_limits<int32_t>::max()};
    inline bool valid_pos(int32_t read_len, uint32_t txp_len, int32_t max_over) {
        int32_t signed_txp_len = static_cast<int32_t>(txp_len);
        return (pos > -max_over) and ((pos + read_len) < (signed_txp_len + max_over));
    }
    inline bool has_mate() const { return mate_pos != invalid_mate_pos; }
    inline bool mate_is_mapped() const { return mate_pos != invalid_mate_pos; }
    inline int32_t frag_len() const {
        return (fragment_length != invalid_frag_len) ? fragment_length : 0;
    }
};

enum class MappingType : uint8_t {
    UNMAPPED = 0,
    SINGLE_MAPPED = 1,
    MAPPED_FIRST_ORPHAN = 2,
    MAPPED_SECOND_ORPHAN = 3,
    MAPPED_PAIR = 4
};

enum class HitDirection : uint8_t { FW, RC, BOTH };

constexpr uint8_t max_distortion = std::numeric_limits<uint8_t>::max();

struct chain_state {
    int32_t read_start_pos{-1};
    int32_t prev_pos{-1};
    int32_t curr_pos{-1};
    uint8_t num_hits{0};
    uint8_t min_distortion{max_distortion};
};

inline bool compare_chains(const chain_state& a, const chain_state& b) {
    return a.prev_pos < b.prev_pos;
}

struct sketch_hit_info {
    static constexpr size_t max_num_chains = 8;
    // add a hit to the current target that occurs in the forward
    // orientation with respect to the target.
    bool add_fw(int32_t ref_pos, int32_t read_pos, int32_t rl, int32_t k, int32_t max_stretch,
                float score_inc) {
        (void)rl;
        (void)k;
        bool added{false};
        int32_t approx_map_pos = ref_pos - read_pos;
        // std::cout << "map_pos" << approx_map_pos << std::endl;
        // If structural constraints have been disabled
        // then simply count the number of hits we see in
        // the given orientation (being careful to count
        // a k-mer of a given rank only one time).
        // std::cout<<"val of ignore struct constraints" <<  ignore_struct_constraints_fw << std::endl;
        if (ignore_struct_constraints_fw) {
            // std::cout << ignore_struct_constraints_fw << std::endl;
            if (read_pos > last_read_pos_fw) {
                if (last_read_pos_fw == -1) { approx_pos_fw = approx_map_pos; }
                last_ref_pos_fw = ref_pos;
                last_read_pos_fw = read_pos;
                fw_score += score_inc;
                ++fw_hits;
                added = true;
            }
            return added;
        }
        // NO STRUCTURAL CONSTRAINTS

        // If this is a k-mer of a new rank
        if (read_pos > last_read_pos_fw) {
            // if this is the first k-mer we are seeing
            if (last_read_pos_fw == -1) {
                approx_pos_fw = approx_map_pos;
            } else {
                // we are seeing a k-mer of a new rank.
                // at this point, we've seen all k-mers of the
                // previous rank, so copy over any valid chains
                // for the next search.
                compact_chains(fw_chains, fw_hits);
            }

            // update the current query position
            // (i.e. rank) of the seed we are processing
            last_read_pos_fw = read_pos;
            ++fw_rank;
        }

        // if this is still a hit for the *first*
        // k-mer of the chains — *note* the first
        // k-mer will be of rank 0 because fw_rank
        // is initialized to -1.
        if (fw_rank == 0) {
            process_rank0_hit(approx_map_pos, ref_pos, fw_chains, approx_pos_fw,
                              ignore_struct_constraints_fw, fw_hits, added);
        } else {
            // this is a hit for a k-mer of rank > 0, so we already
            // have a set of active chains.
            process_hit(true, approx_map_pos, ref_pos, max_stretch, fw_chains, fw_hits, added);
        }
        return added;

        ///// orig
        /*
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
            last_ref_pos_fw = ref_pos;
            last_read_pos_fw = read_pos;
            fw_score += score_inc;
            ++fw_hits;
            added = true;
        }
        return added;
        */
    }

    // add a hit to the current target that occurs in the forward
    // orientation with respect to the target.
    bool add_rc(int32_t ref_pos, int32_t read_pos, int32_t rl, int32_t k, int32_t max_stretch,
                float score_inc) {
        bool added{false};
        int32_t approx_map_pos = (ref_pos - (rl - (read_pos + k)));
        // std::cout << "map_pos" << approx_map_pos << std::endl;
        // NO STRUCTURAL CONSTRAINTS
        if (ignore_struct_constraints_rc) {
            if (read_pos > last_read_pos_rc) {
                approx_pos_rc = approx_map_pos;
                if (last_read_pos_rc == -1) {
                    approx_end_pos_rc = ref_pos + read_pos;
                    first_read_pos_rc = read_pos;
                }
                rc_score += score_inc;
                ++rc_hits;

                // new
                rightmost_bound_rc = last_ref_pos_rc;

                last_ref_pos_rc = ref_pos;
                last_read_pos_rc = read_pos;
                added = true;
            }
            return added;
        }
        // NO STRUCTURAL CONSTRAINTS

        // If this is a k-mer of a new rank
        if (read_pos > last_read_pos_rc) {
            approx_pos_rc = approx_map_pos;
            // if this is the first k-mer we are seeing
            if (last_read_pos_rc == -1) {
                approx_end_pos_rc = ref_pos + read_pos;
                first_read_pos_rc = read_pos;
            } else {
                // point, we've seen all k-mers of the previous rank, so copy over any valid chains
                // for the next search.
                compact_chains(rc_chains, rc_hits);
            }

            // update the current query position
            // (i.e. rank) of the seed we are processing
            ++rc_rank;
            // new
            rightmost_bound_rc = last_ref_pos_rc;

            last_ref_pos_rc = ref_pos;
            last_read_pos_rc = read_pos;
            added = true;
        }

        // if this is still a hit for the *first*
        // k-mer of the chains
        if (rc_rank == 0) {
            process_rank0_hit(approx_map_pos, ref_pos, rc_chains, approx_pos_rc,
                              ignore_struct_constraints_rc, rc_hits, added);
        } else {
            // this is a hit for a k-mer of rank > 0, so we already
            // have a set of active chains.
            process_hit(false, approx_map_pos, ref_pos, max_stretch, rc_chains, rc_hits, added);
        }
        return added;
        /*
        // since hits are collected by moving _forward_ in the
        // read, if this is an rc hit, it should be moving
        // backwards in the reference.
        // In general, we only add the hit if this is the case.
        // This ensures that we don't double-count a k-mer that
        // might occur twice on this target.

        // we have a special case here; what if the same exact
        // k-mer (i.e. not just the same sequence but same position
        // on the query) occurs more than one time on this refernece?
        //
        // In that case, the GENERAL case code will already have
        // processed and seen a k-mer with the read position
        // equal to `read_pos`.  In the case below, we see
        // a hit with the *same* read pos again (one or more times).
        //
        // Here, we swap out the previous hit having this read_pos
        // if the position of the current hit on the read is
        // the same and the position on the reference is greater
        // (this is a heuristic to help in the case of tandem repeats or
        // highly-repetitive subsequence).
        // NOTE: consider if a similar heuristic should be
        // adopted for the forward case.
        if ((read_pos == last_read_pos_rc) and (ref_pos > last_ref_pos_rc) and
            (ref_pos < rightmost_bound_rc)) {

            last_ref_pos_rc = ref_pos;
            // if the read_pos was the same as the first read pos
            // then also update the approx_end_pos_rc accordingly
            // NOTE: for the time being don't mess with this position
            // empirically this does better, but if we really want
            // to optimize this for accuracy we need a better general
            // heuristic.
            // if (read_pos == first_read_pos_rc) {
            //   approx_end_pos_rc = ref_pos + read_pos;
            //}

            return added;
        }

        // GENERAL case
        if (ref_pos < last_ref_pos_rc and read_pos > last_read_pos_rc) {
            approx_pos_rc = (ref_pos - (rl - (read_pos + k)));
            if (last_read_pos_rc == -1) {
                approx_end_pos_rc = ref_pos + read_pos;
                first_read_pos_rc = read_pos;
            } else {
                if (approx_end_pos_rc - approx_pos_rc > max_stretch) { return false; }
            }
            rc_score += score_inc;
            ++rc_hits;

            // new
            rightmost_bound_rc = last_ref_pos_rc;

            last_ref_pos_rc = ref_pos;
            last_read_pos_rc = read_pos;
            added = true;
        }
        return added;
        */
    }

    // for directly incrementing the number of hits
    // even when we are not building chains (e.g. in the case
    // of filtering based on occurrences of ambiguous seeds).
    inline void inc_fw_hits() { ++fw_hits; }
    inline void inc_rc_hits() { ++rc_hits; }

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

    int32_t tid{std::numeric_limits<int32_t>::max()};
    int32_t last_read_pos_fw{-1};
    int32_t last_read_pos_rc{-1};
    int32_t rightmost_bound_rc{std::numeric_limits<int32_t>::max()};

    // marks the read position (key) of the
    // first hit we see in the rc direction
    int32_t first_read_pos_rc{-1};

    int32_t last_ref_pos_fw{-1};
    int32_t last_ref_pos_rc{std::numeric_limits<int32_t>::max()};

    int32_t approx_pos_fw{-1};
    int32_t approx_pos_rc{-1};
    int32_t approx_end_pos_rc{-1};

    uint32_t fw_hits{0};
    uint32_t rc_hits{0};
    float fw_score{0.0};
    float rc_score{0.0};

    bool ignore_struct_constraints_fw{false};
    bool ignore_struct_constraints_rc{false};

    int32_t fw_rank{-1};
    int32_t rc_rank{-1};
    itlib::small_vector<chain_state, max_num_chains> fw_chains;
    itlib::small_vector<chain_state, max_num_chains> rc_chains;

private:
    inline void compact_chains(itlib::small_vector<chain_state, max_num_chains>& chains,
                               const uint32_t required_hits) {
        chains.erase(std::remove_if(chains.begin(), chains.end(),
                                    [required_hits](chain_state& s) -> bool {
                                        // remove this chain if it doesn't satisfy
                                        // the hit constraint.
                                        if (s.num_hits < required_hits) { return true; }
                                        // the current position becomes the
                                        // previous position for kmers of the
                                        // next rank, and the curr pos gets reset
                                        // to -1.
                                        s.prev_pos = s.curr_pos;
                                        s.curr_pos = -1;
                                        return false;
                                    }),
                     chains.end());
    }

    inline void process_rank0_hit(int32_t approx_map_pos, int32_t hit_pos,
                                  itlib::small_vector<chain_state, max_num_chains>& chains,
                                  int32_t& approx_pos_out, bool& ignore_struct_constraints,
                                  uint32_t& num_hits, bool& added) {
        // if there are too many possible chains, just punt
        // and turn off structural constraints for this
        // query, reference, orientation tuple.
        if (chains.size() == chains.capacity()) {
            ignore_struct_constraints = true;
            approx_pos_out = chains.front().read_start_pos;
            num_hits = chains.front().num_hits;
            added = true;
        } else {
            // otherise add this hit
            chains.push_back({approx_map_pos, -1, hit_pos, 1, max_distortion});
            added = true;
        }
        num_hits = 1;
    }

    inline void process_hit(
        bool is_fw_hit,  // we are processing a hit in the forward orientation (otherwise, RC)
        int32_t read_start_pos, int32_t next_hit_pos, int32_t max_stretch,
        itlib::small_vector<chain_state, max_num_chains>& chains, uint32_t& num_hits, bool& added) {
        (void)max_stretch;
        // find the chain that best matches this k-mer.
        chain_state predecessor_probe{read_start_pos, next_hit_pos, -1, 0, max_distortion};
        auto chain_pos =
            std::lower_bound(chains.begin(), chains.end(), predecessor_probe, compare_chains);

        // If this is a forward hit, lower bound will return
        // the first element >= the key (the current hit), so back up
        // by one element.
        // Otherwise, if this is a reverse complement hit, we actually
        // want the first element >= the current hit, so keep it.
        if (is_fw_hit) {
            if (chain_pos > chains.begin()) {
                chain_pos--;
            } else {
                added = false;
                return;
            }
        }

        // if we found a valid chain to extend
        if (chain_pos < chains.end()) {
            auto stretch = std::abs(chain_pos->read_start_pos - read_start_pos);
            stretch = std::min(stretch, static_cast<decltype(stretch)>(max_distortion));
            // and if it hasn't yet been extended
            if (chain_pos->curr_pos == -1) {
                if (stretch < 15) {
                    // then extend this chain
                    chain_pos->curr_pos = next_hit_pos;
                    chain_pos->min_distortion = static_cast<uint8_t>(stretch);
                    // and increment the number of hits
                    ++(chain_pos->num_hits);
                    uint32_t curr_max_hits = num_hits;
                    num_hits = std::max(curr_max_hits,
                                        static_cast<decltype(curr_max_hits)>(chain_pos->num_hits));
                    added = true;
                }
            } else if (stretch < chain_pos->min_distortion) {
                chain_pos->min_distortion = static_cast<uint8_t>(stretch);
                added = true;
            }
        }
    }
};

struct mapping_cache_info {
public:
    mapping_cache_info(mindex::reference_index& ri) : k(ri.k()), q(ri.get_dict()), hs(&ri) {}

    inline void clear() {
        map_type = mapping::util_bin::MappingType::UNMAPPED;
        q.start();
        hs.clear();
        hit_map.clear();
        accepted_hits.clear();
        has_matching_kmers = false;
        ambiguous_hit_indices.clear();
    }

    // will store how the read mapped
    mapping::util_bin::MappingType map_type{mapping::util_bin::MappingType::UNMAPPED};

    // map from reference id to hit info
    phmap::flat_hash_map<uint32_t, mapping::util_bin::sketch_hit_info> hit_map;
    std::vector<mapping::util_bin::simple_hit> accepted_hits;

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
    // holds the indices of k-mers too ambiguous to chain, but which
    // we might later want to check the existence of
    itlib::small_vector<uint32_t, 255> ambiguous_hit_indices;

    // max ec card
    uint32_t max_ec_card{256};
};

inline bool map_read(std::string* read_seq, mapping_cache_info& map_cache, bool& k_match, bool verbose = false) {
    map_cache.clear();
    // rebind map_cache variables to
    // local names
    auto& q = map_cache.q;
    auto& hs = map_cache.hs;
    auto& hit_map = map_cache.hit_map;
    auto& accepted_hits = map_cache.accepted_hits;
    auto& map_type = map_cache.map_type;
    const bool attempt_occ_recover = map_cache.attempt_occ_recover;
    const bool perform_ambig_filtering = map_cache.hs.get_index()->has_ec_table();
    auto k = map_cache.k;

    map_cache.has_matching_kmers = hs.get_raw_hits_sketch(*read_seq, q, true, false);
    bool early_stop = false;

    // if we are checking ambiguous hits, the maximum EC
    // size we will consider.
    const size_t max_ec_ambig = map_cache.max_ec_card;

    // if there were hits
    if (map_cache.has_matching_kmers) {
        k_match = true;
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
             &early_stop, signed_rl, k, &map_cache, perform_ambig_filtering,
             verbose](auto& raw_hits, auto& prev_read_pos, auto& max_allowed_occ,
                      auto& ambiguous_hit_indices, auto& had_alt_max_occ) -> bool {
            int32_t hit_idx{0};

            for (auto& raw_hit : raw_hits) { // for each queried kmer hit
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

                    for (auto v : refs) { // ref range over unitigs
                    // a unitig can appear multiple times, then however a single unitig 
                    // can map to multiple transcripts which does not seem to be the case
                    // based on this code they map to single transcript or it could be that the 
                    // position is w.r.t each color thus it makes sense
                        const auto& ref_pos_ori = proj_hits.decode_hit(v);
                        uint32_t tid = sshash::util::transcript_id(v);
                        int32_t pos = static_cast<int32_t>(ref_pos_ori.pos);
                        bool ori = ref_pos_ori.isFW;
                        
                        auto& target = hit_map[tid];
                        // std::cout << "tid is " << tid;
                        if (verbose) {
                            auto& tname = map_cache.hs.get_index()->ref_name(tid);
                            std::cerr << "\traw_hit [read_pos: " << read_pos << " ]:" << tname
                                      << ", " << pos << ", " << (ori ? "fw" : "rc") << "\n";
                        }
                        

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
                                    target.ignore_struct_constraints_fw = true;
                                    target.add_fw(pos, static_cast<int32_t>(read_pos), signed_rl, k,
                                              max_stretch, score_inc);
                            } else {
                                    target.ignore_struct_constraints_rc = true;
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

                } else if (perform_ambig_filtering) {  // HERE we have that num_occ >
                                                       // max_allowed_occ
                    ambiguous_hit_indices.push_back(hit_idx);
                }

                ++hit_idx;
            }  // DONE : for (auto& raw_hit : raw_hits)

            return false;
        };

        bool _discard = false;
        auto mao_first_pass = map_cache.max_occ_default - 1;
        early_stop = collect_mappings_from_hits(raw_hits, prev_read_pos, mao_first_pass,
                                                map_cache.ambiguous_hit_indices, _discard);

        // If our default threshold was too stringent, then fallback to a more liberal
        // threshold and look up the k-mers that occur the least frequently.
        // Specifically, if the min occuring hits have frequency < max_occ_recover (2500 by
        // default) times, then collect the min occuring hits to get the mapping.
        if (attempt_occ_recover and (min_occ >= map_cache.max_occ_default) and
            (min_occ < map_cache.max_occ_recover)) {
            map_cache.ambiguous_hit_indices.clear();
            prev_read_pos = -1;
            uint64_t max_allowed_occ = min_occ;
            early_stop =
                collect_mappings_from_hits(raw_hits, prev_read_pos, max_allowed_occ,
                                           map_cache.ambiguous_hit_indices, had_alt_max_occ);
        }

        // Further filtering of mappings by ambiguous k-mers
        if (perform_ambig_filtering and !hit_map.empty() and
            !map_cache.ambiguous_hit_indices.empty()) {
            phmap::flat_hash_set<uint64_t> observed_ecs;
            size_t min_cardinality_ec_size = std::numeric_limits<size_t>::max();
            uint64_t min_cardinality_ec = std::numeric_limits<size_t>::max();
            size_t min_cardinality_index = 0;
            size_t visited = 0;
            auto& ec_table = map_cache.hs.get_index()->get_ec_table();

            auto visit_ec = [&hit_map](uint64_t ent, bool fw_on_contig) -> bool {
                uint32_t tid = (ent >> 2);
                auto hm_it = hit_map.find(tid);
                bool found = false;
                if (hm_it != hit_map.end()) {
                    // we found this target, now:
                    // (1) check the orientation
                    uint32_t ori = (ent & 0x3);
                    // (2) add hits in the appropriate way
                    switch (ori) {
                        case 0:  // fw
                            (fw_on_contig) ? hm_it->second.inc_fw_hits()
                                           : hm_it->second.inc_rc_hits();
                            break;
                        case 1:  // rc
                            (fw_on_contig) ? hm_it->second.inc_rc_hits()
                                           : hm_it->second.inc_fw_hits();
                            break;
                        default:  // both
                            hm_it->second.inc_fw_hits();
                            hm_it->second.inc_rc_hits();
                    }
                    found = true;
                }
                return found;
            };

            // for each ambiguous hit
            for (auto hit_idx : map_cache.ambiguous_hit_indices) {
                auto& proj_hit = raw_hits[hit_idx].second;
                uint32_t contig_id = proj_hit.contig_id();
                bool fw_on_contig = proj_hit.hit_fw_on_contig();

                // put the combination of the eq and the k-mer orientation
                // into the map.
                uint64_t ec = ec_table.ec_for_tile(contig_id);
                uint64_t ec_key = ec | (fw_on_contig ? 0 : 0x8000000000000000);
                // if we've already seen this ec, no point in processing
                // it again.
                if (observed_ecs.contains(ec_key)) { continue; }
                // otherwise, insert it.
                observed_ecs.insert(ec_key);

                auto ec_entries = ec_table.entries_for_ec(ec);
                if (ec_entries.size() < min_cardinality_ec_size) {
                    min_cardinality_ec_size = ec_entries.size();
                    min_cardinality_ec = ec;
                    min_cardinality_index = hit_idx;
                }
                if (ec_entries.size() > max_ec_ambig) { continue; }
                ++visited;
                for (const auto& ent : ec_entries) {
                    visit_ec(ent, fw_on_contig);
                }  // all target oritentation pairs in this eq class
                ++num_valid_hits;
            }  // all ambiguous hits

            // if we haven't visited *any* equivalence classes (they were all)
            // too ambiguous, then make a last-ditch effort to just visit the
            // the one with smallest cardinality.
            if (visited == 0) {
                auto hit_idx = min_cardinality_index;
                auto& proj_hit = raw_hits[hit_idx].second;
                bool fw_on_contig = proj_hit.hit_fw_on_contig();

                uint64_t ec = min_cardinality_ec;
                auto ec_entries = ec_table.entries_for_ec(ec);
                for (const auto& ent : ec_entries) {
                    visit_ec(ent, fw_on_contig);
                }  // all target oritentation pairs in this eq class
                ++num_valid_hits;
            }  // done visiting the last-ditch ec

        }  // if we are processing ambiguous hits

        uint32_t best_alt_hits = 0;
        // int32_t signed_read_len = static_cast<int32_t>(record.seq.length());

        for (auto& kv : hit_map) {
            auto best_hit_dir = kv.second.best_hit_direction();

            // if the best direction is FW or BOTH, add the fw hit
            // otherwise add the RC.
            auto simple_hit = (best_hit_dir != mapping::util_bin::HitDirection::RC)
                                  ? kv.second.get_fw_hit()
                                  : kv.second.get_rc_hit();

            if (simple_hit.num_hits >= num_valid_hits) {
                simple_hit.tid = kv.first;
                accepted_hits.emplace_back(simple_hit);
                // if we had equally good hits in both directions
                // add the rc hit here (since we added the fw)
                // above if the best hit was either FW or BOTH
                if (best_hit_dir == mapping::util_bin::HitDirection::BOTH) {
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
        map_type = mapping::util_bin::MappingType::UNMAPPED;
    } else if (!accepted_hits.empty()) {
        map_type = mapping::util_bin::MappingType::SINGLE_MAPPED;
    }

    return early_stop;
}

inline bool map_atac_read(std::string* read_seq, mapping_cache_info& map_cache,
    bool verbose, bool& k_match, mindex::reference_index& ri, bool psc_off=false, bool ps_skip=true, float thr=1.0) {
    map_cache.clear();
    // rebind map_cache variables to
    // local names
    auto& q = map_cache.q;
    auto& hs = map_cache.hs;
    auto& hit_map = map_cache.hit_map;
    auto& accepted_hits = map_cache.accepted_hits;
    auto& map_type = map_cache.map_type;
    const bool attempt_occ_recover = map_cache.attempt_occ_recover;
    const bool perform_ambig_filtering = map_cache.hs.get_index()->has_ec_table();
    auto k = map_cache.k;

    
    if (ps_skip) {
        map_cache.has_matching_kmers = hs.get_raw_hits_sketch(*read_seq, q, true, verbose);
    }
    else {
        map_cache.has_matching_kmers = hs.get_raw_hits_sketch_everykmer(*read_seq, q, true, verbose);
    }
    
    bool early_stop = false;

    // if we are checking ambiguous hits, the maximum EC
    // size we will consider.
    const size_t max_ec_ambig = map_cache.max_ec_card;

    // if there were hits
    if (map_cache.has_matching_kmers) {
        k_match = true;
        float num_valid_hits{0};
        uint64_t total_occs{0};
        uint64_t largest_occ{0};
        auto& raw_hits = hs.get_left_hits();
        // std::cout << "raw hits" << raw_hits.size() << std::endl;
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
             &early_stop, signed_rl, k, &map_cache, perform_ambig_filtering,
             verbose, psc_off, ps_skip](auto& raw_hits, auto& prev_read_pos, auto& max_allowed_occ,
                      auto& ambiguous_hit_indices, auto& had_alt_max_occ) -> bool {
            int32_t hit_idx{0};
            // return false;

            for (auto& raw_hit : raw_hits) {
                auto& read_pos = raw_hit.first;
                auto& proj_hits = raw_hit.second;
                auto& refs = proj_hits.refRange;
                
                uint64_t num_occ = static_cast<uint64_t>(refs.size());
                min_occ = std::min(min_occ, num_occ);
                had_alt_max_occ = true;

                bool still_have_valid_target = false;
                prev_read_pos = read_pos;
                // std::cout << raw_hits.size() << "entered\n";
                if (num_occ <= max_allowed_occ) {
                    total_occs += num_occ;
                    largest_occ = (num_occ > largest_occ) ? num_occ : largest_occ;
                    float score_inc = 1.0;
                    
                    for (auto v : refs) {
                        const auto& ref_pos_ori = proj_hits.decode_hit(v);
                        uint32_t tid = sshash::util::transcript_id(v);
                        int32_t pos = static_cast<int32_t>(ref_pos_ori.pos);
                        bool ori = ref_pos_ori.isFW;
                        
                        auto& target = hit_map[tid];
                        
                        /*
                        if (verbose) {
                            auto& tname = map_cache.hs.get_index()->ref_name(tid);
                            std::cerr << "\traw_hit [read_pos: " << read_pos << " ]:" << tname
                                      << ", " << pos << ", " << (ori ? "fw" : "rc") << "\n";
                        }
                        */

                    //     // Why >= here instead of == ?
                    //     // Because hits can happen on the same target in both the forward
                    //     // and rc orientations, it is possible that we start the loop with
                    //     // the target having num_valid_hits hits in a given orientation (o)
                    //     // we see a new hit for this target in oriention o (now it has
                    //     // num_valid_hits + 1) then we see a hit for this target in
                    //     // orientation rc(o).  We still want to add / consider this hit, but
                    //     // max_hits_for_target() > num_valid_hits. So, we must allow for
                    //     // that here.
                        if (target.max_hits_for_target() >= num_valid_hits) {
                            
                            if (ori) {
                                if (psc_off) {
                                    target.ignore_struct_constraints_fw = true;
                                }
                                target.add_fw(pos, static_cast<int32_t>(read_pos), signed_rl, k,
                                            max_stretch, score_inc);
                            } else {
                                if (psc_off) {
                                    target.ignore_struct_constraints_rc = true;
                                }
                                    
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

                } else if (perform_ambig_filtering) {  // HERE we have that num_occ >
                                                       // max_allowed_occ
                    ambiguous_hit_indices.push_back(hit_idx);
                }

                ++hit_idx;
            }  // DONE : for (auto& raw_hit : raw_hits)

            return false;
        };

        auto collect_mappings_from_hits_thr =
            [&max_stretch, &min_occ, &hit_map, &num_valid_hits, &total_occs, &largest_occ,
             &early_stop, signed_rl, k, &map_cache, perform_ambig_filtering,
             verbose, psc_off, ps_skip, &thr, &ri](auto& raw_hits, auto& prev_read_pos, auto& max_allowed_occ,
                      auto& ambiguous_hit_indices, auto& had_alt_max_occ, auto& thr, auto &ri) -> bool {
            int32_t hit_idx{0};
            // return false;
            for (auto& raw_hit : raw_hits) {
                auto& read_pos = raw_hit.first;
                auto& proj_hits = raw_hit.second;
                auto& refs = proj_hits.refRange;
                
                uint64_t num_occ = static_cast<uint64_t>(refs.size());
                min_occ = std::min(min_occ, num_occ);
                had_alt_max_occ = true;

                bool still_have_valid_target = false;
                prev_read_pos = read_pos;
                
                // std::cout << "proj hits " << num_occ << std::endl;
                // std::cout << raw_hits.size() << "entered\n";
                if (num_occ <= max_allowed_occ) {
                    total_occs += num_occ;
                    float score_inc = 1.0;
                    
                    for (auto v : refs) {
                        const auto& ref_pos_ori = proj_hits.decode_hit(v);
                        uint32_t tid = sshash::util::transcript_id(v);
                        int32_t pos = static_cast<int32_t>(ref_pos_ori.pos);
                        std::pair<uint32_t, uint32_t> bins = get_bin_id(pos); 
                        bool ori = ref_pos_ori.isFW;
                        // std::cout << "tid " << tid << std::endl;
                        
                        auto& target1 = hit_map[bins.first];
                        target1.tid = tid;
                        auto& target2 = hit_map[bins.second];
                        /*
                        if (verbose) {
                            auto& tname = map_cache.hs.get_index()->ref_name(tid);
                            std::cerr << "\traw_hit [read_pos: " << read_pos << " ]:" << tname
                                      << ", " << pos << ", " << (ori ? "fw" : "rc") << "\n";
                        }
                        */

                        // Why >= here instead of == ?
                        // Because hits can happen on the same target in both the forward
                        // and rc orientations, it is possible that we start the loop with
                        // the target having num_valid_hits hits in a given orientation (o)
                        // we see a new hit for this target in oriention o (now it has
                        // num_valid_hits + 1) then we see a hit for this target in
                        // orientation rc(o).  We still want to add / consider this hit, but
                        // max_hits_for_target() > num_valid_hits. So, we must allow for
                        // that here.
                        
                            
                        if (ori) {
                            if (psc_off) {
                                target1.ignore_struct_constraints_fw = true;
                            }
                            // std::cout << "txp " << ri.ref_name(tid) << std::endl;
                            target1.add_fw(pos, static_cast<int32_t>(read_pos), signed_rl, k,
                                        max_stretch, score_inc);
                        } else {
                            if (psc_off) {
                                target1.ignore_struct_constraints_rc = true;
                            }
                            // std::cout << "txp " << ri.ref_name(tid) << std::endl;    
                            target1.add_rc(pos, static_cast<int32_t>(read_pos), signed_rl, k,
                                        max_stretch, score_inc);
                        }
                    
                        if (bins.second!=std::numeric_limits<uint32_t>::max()) {
                            target2.tid = tid;
                            if (ori) {
                                if (psc_off) {
                                    target2.ignore_struct_constraints_fw = true;
                                }
                                target2.add_fw(pos, static_cast<int32_t>(read_pos), signed_rl, k,
                                        max_stretch, score_inc);
                            } else {
                                if (psc_off) {
                                    target2.ignore_struct_constraints_rc = true;
                                }
                                // std::cout << "txp " << ri.ref_name(tid) << std::endl;    
                                target2.add_rc(pos, static_cast<int32_t>(read_pos), signed_rl, k,
                                            max_stretch, score_inc);
                            }       
                        }
                    }  // DONE: for (auto &pos_it : refs)
              
                    ++num_valid_hits;

                } else if (perform_ambig_filtering) {  // HERE we have that num_occ >
                                                       // max_allowed_occ
                    ambiguous_hit_indices.push_back(hit_idx);
                }

                ++hit_idx;
            }  // DONE : for (auto& raw_hit : raw_hits)

            return false;
        };

        bool _discard = false;
        auto mao_first_pass = map_cache.max_occ_default - 1;
        if (thr==1.0) {
            early_stop = collect_mappings_from_hits(raw_hits, prev_read_pos, mao_first_pass,
                                        map_cache.ambiguous_hit_indices, _discard);
            // std::cout << early_stop << "thr\n";                                        
            // If our default threshold was too stringent, then fallback to a more liberal
            // threshold and look up the k-mers that occur the least frequently.
            // Specifically, if the min occuring hits have frequency < max_occ_recover (2500 by
            // default) times, then collect the min occuring hits to get the mapping.
            if (attempt_occ_recover and (min_occ >= map_cache.max_occ_default) and
                (min_occ < map_cache.max_occ_recover)) {
                map_cache.ambiguous_hit_indices.clear();
                prev_read_pos = -1;
                uint64_t max_allowed_occ = min_occ;
                early_stop =
                    collect_mappings_from_hits(raw_hits, prev_read_pos, max_allowed_occ,
                                            map_cache.ambiguous_hit_indices, had_alt_max_occ);
            }            
        
        }

        else {
            early_stop = collect_mappings_from_hits_thr(raw_hits, prev_read_pos, mao_first_pass,
                                       map_cache.ambiguous_hit_indices, _discard, thr, ri);
            // std::cout << early_stop << "thr\n"                           ;
            // If our default threshold was too stringent, then fallback to a more liberal
            // threshold and look up the k-mers that occur the least frequently.
            // Specifically, if the min occuring hits have frequency < max_occ_recover (2500 by
            // default) times, then collect the min occuring hits to get the mapping.
            if (attempt_occ_recover and (min_occ >= map_cache.max_occ_default) and
                (min_occ < map_cache.max_occ_recover)) {
                map_cache.ambiguous_hit_indices.clear();
                prev_read_pos = -1;
                uint64_t max_allowed_occ = min_occ;
                early_stop =
                    collect_mappings_from_hits_thr(raw_hits, prev_read_pos, max_allowed_occ,
                                    map_cache.ambiguous_hit_indices, had_alt_max_occ, thr, ri);
            }
        }


        
        // Further filtering of mappings by ambiguous k-mers
        if (perform_ambig_filtering and !hit_map.empty() and
            !map_cache.ambiguous_hit_indices.empty()) {
            phmap::flat_hash_set<uint64_t> observed_ecs;
            size_t min_cardinality_ec_size = std::numeric_limits<size_t>::max();
            uint64_t min_cardinality_ec = std::numeric_limits<size_t>::max();
            size_t min_cardinality_index = 0;
            size_t visited = 0;
            auto& ec_table = map_cache.hs.get_index()->get_ec_table();

            auto visit_ec = [&hit_map](uint64_t ent, bool fw_on_contig) -> bool {
                uint32_t bin_id = (ent >> 2);
                auto hm_it = hit_map.find(bin_id);
                bool found = false;
                if (hm_it != hit_map.end()) {
                    // we found this target, now:
                    // (1) check the orientation
                    uint32_t ori = (ent & 0x3);
                    // (2) add hits in the appropriate way
                    switch (ori) {
                        case 0:  // fw
                            (fw_on_contig) ? hm_it->second.inc_fw_hits()
                                           : hm_it->second.inc_rc_hits();
                            break;
                        case 1:  // rc
                            (fw_on_contig) ? hm_it->second.inc_rc_hits()
                                           : hm_it->second.inc_fw_hits();
                            break;
                        default:  // both
                            hm_it->second.inc_fw_hits();
                            hm_it->second.inc_rc_hits();
                    }
                    found = true;
                }
                return found;
            };

            // for each ambiguous hit
            for (auto hit_idx : map_cache.ambiguous_hit_indices) {
                auto& proj_hit = raw_hits[hit_idx].second;
                uint32_t contig_id = proj_hit.contig_id();
                bool fw_on_contig = proj_hit.hit_fw_on_contig();

                // put the combination of the eq and the k-mer orientation
                // into the map.
                uint64_t ec = ec_table.ec_for_tile(contig_id);
                uint64_t ec_key = ec | (fw_on_contig ? 0 : 0x8000000000000000);
                // if we've already seen this ec, no point in processing
                // it again.
                if (observed_ecs.contains(ec_key)) { continue; }
                // otherwise, insert it.
                observed_ecs.insert(ec_key);

                auto ec_entries = ec_table.entries_for_ec(ec);
                if (ec_entries.size() < min_cardinality_ec_size) {
                    min_cardinality_ec_size = ec_entries.size();
                    min_cardinality_ec = ec;
                    min_cardinality_index = hit_idx;
                }
                if (ec_entries.size() > max_ec_ambig) { continue; }
                ++visited;
                for (const auto& ent : ec_entries) {
                    visit_ec(ent, fw_on_contig);
                }  // all target oritentation pairs in this eq class
                ++num_valid_hits;
            }  // all ambiguous hits

            // if we haven't visited *any* equivalence classes (they were all)
            // too ambiguous, then make a last-ditch effort to just visit the
            // the one with smallest cardinality.
            if (visited == 0) {
                auto hit_idx = min_cardinality_index;
                auto& proj_hit = raw_hits[hit_idx].second;
                bool fw_on_contig = proj_hit.hit_fw_on_contig();

                uint64_t ec = min_cardinality_ec;
                auto ec_entries = ec_table.entries_for_ec(ec);
                for (const auto& ent : ec_entries) {
                    visit_ec(ent, fw_on_contig);
                }  // all target oritentation pairs in this eq class
                ++num_valid_hits;
            }  // done visiting the last-ditch ec

        }  // if we are processing ambiguous hits

        uint32_t best_alt_hits = 0;
        // int32_t signed_read_len = static_cast<int32_t>(record.seq.length());
        // std::cout << "before" << num_valid_hits << "\n";
        if(thr != 1) {
            num_valid_hits = num_valid_hits*thr;
        }
        // std::cout << "after" << num_valid_hits << "\n";
        for (auto& kv : hit_map) {
            auto best_hit_dir = kv.second.best_hit_direction();

            // if the best direction is FW or BOTH, add the fw hit
            // otherwise add the RC.
            auto simple_hit = (best_hit_dir != mapping::util_bin::HitDirection::RC)
                                  ? kv.second.get_fw_hit()
                                  : kv.second.get_rc_hit();
            
            if (simple_hit.num_hits >= num_valid_hits) {
                // std::cout << "enter valid\n";
                simple_hit.bin_id = kv.first;
                simple_hit.tid = kv.second.tid;
                accepted_hits.emplace_back(simple_hit);
                // std::cout << "simple hits" << simple_hit.num_hits << " " << num_valid_hits << " tid " << simple_hit.tid << std::endl;
                // if we had equally good hits in both directions
                // add the rc hit here (since we added the fw)
                // above if the best hit was either FW or BOTH
                if (best_hit_dir == mapping::util_bin::HitDirection::BOTH) {
                    auto second_hit = kv.second.get_rc_hit();
                    second_hit.bin_id = kv.first;
                    second_hit.tid = kv.second.tid;
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
    // std::cout << accepted_hits.size() << " " << map_cache.alt_max_occ << "\n";
    if (accepted_hits.size() > map_cache.alt_max_occ) {
        accepted_hits.clear();
        map_type = mapping::util_bin::MappingType::UNMAPPED;
    } else if (!accepted_hits.empty()) {
        map_type = mapping::util_bin::MappingType::SINGLE_MAPPED;
    }

    return early_stop;
}

void print_hits(const std::vector<mapping::util_bin::simple_hit> &hits ) {
    for (const auto& hit : hits) {
        std::cout << "isFw: " << hit.is_fw << std::endl;
        std::cout << "pos: " << hit.pos << std::endl;
        std::cout << "num hits: " << hit.num_hits << std::endl;
        std::cout << "tid: " << hit.tid << std::endl;
        std::cout << "bin_id: " << hit.bin_id << std::endl;
        std::cout << "------------------------" << std::endl;
    }
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
    // std::cout << "num_accepted " << num_accepted_left << " " << 
    //     "num_accepted_right" << num_accepted_right << "\n";
    if ((num_accepted_left > 0) and (num_accepted_right > 0)) {
        // std::cout << "entered both\n";
        // print_hits(accepted_left);
        // std::cout << "left right\n";
        // print_hits(accepted_right);
        // look for paired end mappings
        // so we have to sort our accepted hits
        struct {
            // sort first by orientation, then by transcript id, and finally by position
            bool operator()(const mapping::util_bin::simple_hit& a,
                            const mapping::util_bin::simple_hit& b) {
                if (a.is_fw != b.is_fw) { return a.is_fw > b.is_fw; }
                // orientations are the same
                if (a.tid != b.tid) { return a.tid < b.tid; }
                return a.pos < b.pos;
            }
        } simple_hit_less;
        std::sort(accepted_left.begin(), accepted_left.end(), simple_hit_less);
        std::sort(accepted_right.begin(), accepted_right.end(), simple_hit_less);

        const mapping::util_bin::simple_hit smallest_rc_hit = {false, false, -1, 0.0, 0, 0};
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

        auto merge_lists = [left_len, right_len](iter_t first1, iter_t last1, iter_t first2,
                                                 iter_t last2, out_iter_t out) -> out_iter_t {
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
                        if (frag_len == 0) {
                            // std::cout << "0 fragment length";
                        }
                        // std::cout << first1->tid << " tids " << first2->tid << "pos " << pos_fw << " " << pos_rc << " frag len " << frag_len << "\n";
                        if ((-20 < frag_len) and (frag_len < 1000)) {
                            // if left is fw and right is rc then
                            // fragment length is (right_pos + right_len - left_pos) + 1
                            // otherwise it is (left_pos + left_len - right_pos) + 1
                            bool right_is_rc = !first2->is_fw;
                            int32_t tlen = right_is_rc
                                               ? ((first2->pos + right_len - first1->pos) + 1)
                                               : ((first1->pos + left_len - first2->pos) + 1);
                            *out++ = {first1->is_fw, first2->is_fw, first1->pos, 0.0, 0,
                                      first1->tid,   first2->pos,   tlen};
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

        map_cache_out.map_type = (map_cache_out.accepted_hits.size() > 0) ? MappingType::MAPPED_PAIR
                                                                          : MappingType::UNMAPPED;
        // std::cout << "map cache " << map_cache_out.accepted_hits.size() << std::endl;
    } else if ((num_accepted_left > 0) and !had_matching_kmers_right) {
        // just return the left mappings
        std::swap(map_cache_left.accepted_hits, map_cache_out.accepted_hits);
        map_cache_out.map_type = (map_cache_out.accepted_hits.size() > 0)
                                     ? MappingType::MAPPED_FIRST_ORPHAN
                                     : MappingType::UNMAPPED;
    } else if ((num_accepted_right > 0) and !had_matching_kmers_left) {
        // just return the right mappings
        std::swap(map_cache_right.accepted_hits, map_cache_out.accepted_hits);
        map_cache_out.map_type = (map_cache_out.accepted_hits.size() > 0)
                                     ? MappingType::MAPPED_SECOND_ORPHAN
                                     : MappingType::UNMAPPED;
    } else {
        // return nothing
    }
}

inline void merge_se_mappings(mapping_cache_info& map_cache_left,
                              mapping_cache_info& map_cache_right, int32_t left_len,
                              int32_t right_len, mapping_cache_info& map_cache_out, 
                              mindex::reference_index& ri) {
    map_cache_out.clear();
    auto& accepted_left = map_cache_left.accepted_hits;
    auto& accepted_right = map_cache_right.accepted_hits;

    size_t had_matching_kmers_left = map_cache_left.has_matching_kmers;
    size_t had_matching_kmers_right = map_cache_right.has_matching_kmers;

    size_t num_accepted_left = accepted_left.size();
    size_t num_accepted_right = accepted_right.size();
    // std::cout << "num hits " << num_accepted_left << " " << num_accepted_right << "\n";
    // std::cout << "matching kmers " << had_matching_kmers_left << " " << had_matching_kmers_right << "\n";
    if ((num_accepted_left > 0) and (num_accepted_right > 0)) {
        // std::cout << "entered both\n";
        // print_hits(accepted_left);
        // std::cout << "left right\n";
        // print_hits(accepted_right);
        // look for paired end mappings
        // so we have to sort our accepted hits
        struct {
            // sort first by orientation, then by transcript id, and finally by position
            bool operator()(const mapping::util_bin::simple_hit& a,
                            const mapping::util_bin::simple_hit& b) {
                if (a.is_fw != b.is_fw) { return a.is_fw > b.is_fw; }
                // orientations are the same
                if (a.bin_id != b.bin_id) { return a.bin_id < b.bin_id; }
                return a.pos < b.pos;
            }
        } simple_hit_less;
        std::sort(accepted_left.begin(), accepted_left.end(), simple_hit_less);
        std::sort(accepted_right.begin(), accepted_right.end(), simple_hit_less);

        const mapping::util_bin::simple_hit smallest_rc_hit = {false, false, -1, 0.0, 0, 0, 0};
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
        // auto merge_lists = [left_len, right_len](iter_t first1, iter_t last1, iter_t first2,
        //                                          iter_t last2, out_iter_t out) -> out_iter_t {
        auto merge_lists = [left_len, right_len, ri](iter_t first1, iter_t last1, iter_t first2,
                                                 iter_t last2, out_iter_t out) -> out_iter_t {
            // https://en.cppreference.com/w/cpp/algorithm/set_intersection
            while (first1 != last1 && first2 != last2) {
                if (first1->bin_id < first2->bin_id) {
                    ++first1;
                } else {
                    if (!(first2->bin_id < first1->bin_id)) {
                        // first1->tid == first2->tid have the same transcript.
                        int32_t pos_fw = first1->is_fw ? first1->pos : first2->pos;
                        int32_t pos_rc = first1->is_fw ? first2->pos : first1->pos;
                        int32_t frag_len = (pos_rc - pos_fw);
                        // std::cout << frag_len << " fragment length\n";
                        if (frag_len == 0) {
                            // std::cout << "0 fragment length";
                        }
                        // std::cout << ri.ref_name(first1->bin_id) << " tids " << ri.ref_name(first2->bin_id) << "pos " << pos_fw << " " << pos_rc << " frag len " << frag_len << "\n";
                        if ((-20 < frag_len) and (frag_len < 1000)) {
                            // if left is fw and right is rc then
                            // fragment length is (right_pos + right_len - left_pos) + 1
                            // otherwise it is (left_pos + left_len - right_pos) + 1
                            bool right_is_rc = !first2->is_fw;
                            int32_t tlen = right_is_rc
                                               ? ((first2->pos + right_len - first1->pos) + 1)
                                               : ((first1->pos + left_len - first2->pos) + 1);
                            *out++ = {first1->is_fw, first2->is_fw, first1->pos, 0.0, 0,
                                      first1->tid, first1->bin_id,   first2->pos,   tlen};
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

        map_cache_out.map_type = (map_cache_out.accepted_hits.size() > 0) ? MappingType::MAPPED_PAIR
                                                                          : MappingType::UNMAPPED;
        // std::cout << "map cache sizes" << map_cache_out.accepted_hits.size() << std::endl;
        // std::cout << "map cache2 " << map_cache_out.accepted_hits[0] << std::endl;
    } else if ((num_accepted_left > 0) and !had_matching_kmers_right) {
        // just return the left mappings
        std::swap(map_cache_left.accepted_hits, map_cache_out.accepted_hits);
        map_cache_out.map_type = (map_cache_out.accepted_hits.size() > 0)
                                     ? MappingType::MAPPED_FIRST_ORPHAN
                                     : MappingType::UNMAPPED;
    } else if ((num_accepted_right > 0) and !had_matching_kmers_left) {
        // just return the right mappings
        std::swap(map_cache_right.accepted_hits, map_cache_out.accepted_hits);
        map_cache_out.map_type = (map_cache_out.accepted_hits.size() > 0)
                                     ? MappingType::MAPPED_SECOND_ORPHAN
                                     : MappingType::UNMAPPED;
    } else {
        // return nothing
    }
    // std::cout << "hits right\n";
    // print_hits(map_cache_right.accepted_hits);
    // std::cout << "hits left\n";
    // print_hits(map_cache_left.accepted_hits);
}

}  // namespace util

}  // namespace mapping
