#include "../external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "../include/reference_index.hpp"
#include "../include/CanonicalKmerIterator.hpp"
//#include "../include/query/contig_info_query_canonical_parsing.cpp"
#include "../include/query/streaming_query_canonical_parsing.hpp"
#include "../include/projected_hits.hpp"
#include "../include/util.hpp"
#include "../include/mapping_util.hpp"
#include "../include/spdlog/spdlog.h"
#include "../include/spdlog/sinks/stdout_color_sinks.h"
#include "../include/parallel_hashmap/phmap.h"
#include "../include/cli11/CLI11.hpp"
#include "../include/FastxParser.hpp"
#include "FastxParser.cpp"
#include "hit_searcher.cpp"
#include "zlib.h"

#include <atomic>
#include <iostream>
#include <vector>
#include <numeric>
#include <cstdio>
#include <thread>
#include <sstream>

using namespace klibpp;

void print_header(mindex::reference_index& ri, std::string& cmdline) {
    std::cout << "@HD\tVN:1.0\tSO:unsorted\n";
    for (uint64_t i = 0; i < ri.num_refs(); ++i) {
        std::cout << "@SQ\tNS:" << ri.ref_name(i) << "\tLN:" << ri.ref_len(i) << "\n";
    }
    std::cout << "@PG\tID:mindex_map\tPN:mapper\tVN:0.0.1\t"
              << "CL:" << cmdline << "\n";
}

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

bool is_fw{false};
int32_t pos{-1};
float score{0.0};
uint32_t num_hits{0};
uint32_t tid{std::numeric_limits<uint32_t>::max()};

void merge_se_mappings(mapping_cache_info& map_cache_left, mapping_cache_info& map_cache_right,
                       mapping_cache_info& map_cache_out) {
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

        const mapping::util::simple_hit smallest_rc_hit = {false, -1, 0.0, 0, 0};
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

        auto merge_lists = [](iter_t first1, iter_t last1, iter_t first2, iter_t last2,
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
                            *out++ = {first1->is_fw, pos_fw, 0.0, 0, first1->tid};
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
    } else if ((num_accepted_left > 0) and !had_matching_kmers_right) {
        // just return the left mappings
        std::swap(map_cache_left.accepted_hits, map_cache_out.accepted_hits);
    } else if ((num_accepted_right > 0) and !had_matching_kmers_left) {
        // just return the right mappings
        std::swap(map_cache_right.accepted_hits, map_cache_out.accepted_hits);
    } else {
        // return nothing
    }
}

// single-end
bool map_fragment(fastx_parser::ReadSeq& record, mapping_cache_info& map_cache_left,
                  mapping_cache_info& map_cache_right, mapping_cache_info& map_cache_out) {
    (void)map_cache_left;
    (void)map_cache_right;
    return map_read(&record.seq, map_cache_out);
}

// paried-end
bool map_fragment(fastx_parser::ReadPair& record, mapping_cache_info& map_cache_left,
                  mapping_cache_info& map_cache_right, mapping_cache_info& map_cache_out) {
    bool early_exit_left = map_read(&record.first.seq, map_cache_left);
    bool early_exit_right = map_read(&record.second.seq, map_cache_right);

    merge_se_mappings(map_cache_left, map_cache_right, map_cache_out);

    return (early_exit_left or early_exit_right);
}

inline void write_sam_mappings(mapping_cache_info& map_cache_out,
                               fastx_parser::ReadSeq& record, std::string& workstr_left,
                               std::string& workstr_right, std::atomic<uint64_t>& global_nhits,
                               std::ostringstream& osstream) {
    (void)workstr_right;
    constexpr uint16_t is_secondary = 256;
    constexpr uint16_t is_rc = 16;

    if (!map_cache_out.accepted_hits.empty()) {
        ++global_nhits;
        bool secondary = false;
        for (auto& ah : map_cache_out.accepted_hits) {
            uint16_t flag = secondary ? is_secondary : 0;
            // flag += 2;
            flag += ah.is_fw ? 0 : is_rc;
            // flag += first_seg;

            std::string* sptr = nullptr;
            if (is_rc) {
                combinelib::kmers::reverseComplement(record.seq, workstr_left);
                sptr = &workstr_left;
            } else {
                sptr = &record.seq;
            }
            osstream << record.name << "\t" << flag << "\t" << map_cache_out.hs.get_index()->ref_name(ah.tid) << "\t"
                     << ah.pos + 1 << "\t255\t*\t*\t0\t" << record.seq.length() << "\t" << *sptr
                     << "\t*\n";
            secondary = true;
        }
    } else {
        osstream << record.name << "\t" << 4 << "\t"
                 << "*\t0\t0\t*\t*\t0\t0\t" << record.seq << "\t*\n";
    }
}

// paired-end
inline void write_sam_mappings(mapping_cache_info& map_cache_out,
                               fastx_parser::ReadPair& record, std::string& workstr_left,
                               std::string& workstr_right, std::atomic<uint64_t>& global_nhits,
                               std::ostringstream& osstream) {
    (void)workstr_right;
    constexpr uint16_t is_secondary = 256;
    constexpr uint16_t is_rc = 16;

    if (!map_cache_out.accepted_hits.empty()) {
        ++global_nhits;
        bool secondary = false;
        for (auto& ah : map_cache_out.accepted_hits) {
            uint16_t flag = secondary ? is_secondary : 0;
            // flag += 2;
            flag += ah.is_fw ? 0 : is_rc;
            // flag += first_seg;

            std::string* sptr = nullptr;
            if (is_rc) {
                combinelib::kmers::reverseComplement(record.first.seq, workstr_left);
                sptr = &workstr_left;
            } else {
                sptr = &record.first.seq;
            }
            osstream << record.first.name << "\t" << flag << "\t" << map_cache_out.hs.get_index()->ref_name(ah.tid) << "\t"
                     << ah.pos + 1 << "\t255\t*\t*\t0\t" << record.first.seq.length() << "\t" << *sptr
                     << "\t*\n";
            secondary = true;
        }
    } else {
        osstream << record.first.name << "\t" << 4 << "\t"
                 << "*\t0\t0\t*\t*\t0\t0\t" << record.first.seq << "\t*\n";
    }
}

template <typename FragT>
void do_map(mindex::reference_index& ri, fastx_parser::FastxParser<FragT>& parser,
            std::atomic<uint64_t>& global_nr, std::atomic<uint64_t>& global_nhits,
            std::mutex& iomut) {
    CanonicalKmer::k(ri.k());

    mapping_cache_info map_cache_left(ri);
    mapping_cache_info map_cache_right(ri);
    mapping_cache_info map_cache_out(ri);

    sshash::streaming_query_canonical_parsing q(ri.get_dict());
    mindex::hit_searcher hs(&ri);
    uint64_t read_num = 0;
    uint64_t processed = 0;
    uint64_t buff_size = 10000;
    
    // these don't really belong here
    std::string workstr_left;
    std::string workstr_right;

    std::ostringstream osstream;

    // Get the read group by which this thread will
    // communicate with the parser (*once per-thread*)
    auto rg = parser.getReadGroup();

    while (parser.refill(rg)) {
        // Here, rg will contain a chunk of read pairs
        // we can process.
        for (auto& record : rg) {
            ++global_nr;
            ++read_num;
            auto rctr = global_nr.load();
            auto hctr = global_nhits.load();
            if (rctr % 100000 == 0) {
                std::cerr << "readnum : " << rctr << ", num_hits : " << hctr << "\n";
            }

            // this *overloaded* function will just do the right thing.
            // If record is single-end, just map that read, otherwise, map both and look
            // for proper pairs.
            bool had_early_stop =
                map_fragment(record, map_cache_left, map_cache_right, map_cache_out);
            (void)had_early_stop;

            write_sam_mappings(map_cache_out, record, workstr_left, workstr_right, global_nhits,
                                osstream);

            if (processed >= buff_size) {
                std::string o = osstream.str();
                iomut.lock();
                std::cout << o;
                iomut.unlock();
                osstream.clear();
                osstream.str("");
                processed = 0;
            }
        }
    }

    // dump any remaining output
    std::string o = osstream.str();
    iomut.lock();
    std::cout << o;
    iomut.unlock();
    osstream.clear();
    // don't need this here because osstream goes away at end of scope
    // osstream.str("");
}

int main(int argc, char** argv) {
    /**
     * Mapper
     **/
    std::ios_base::sync_with_stdio(false);

    std::string index_basename;
    std::vector<std::string> left_read_filenames;
    std::vector<std::string> right_read_filenames;
    std::vector<std::string> single_read_filenames;
    size_t nthread{16};
    bool quiet{false};

    CLI::App app{"Mapper"};
    app.add_option("-i,--index", index_basename, "input index prefix")->required();

    auto ogroup = app.add_option_group("input reads", "provide input reads");

    CLI::Option* read_opt =
        ogroup->add_option("-r,--reads", single_read_filenames, "path to list of single-end files")
            ->delimiter(',');
    CLI::Option* paired_left_opt =
        ogroup->add_option("-1,--read1", left_read_filenames, "path to list of read 1 files")
            ->delimiter(',');
    CLI::Option* paired_right_opt =
        ogroup->add_option("-2,--read2", right_read_filenames, "path to list of read 2 files")
            ->delimiter(',');

    paired_left_opt->excludes(read_opt);
    paired_right_opt->excludes(read_opt);
    read_opt->excludes(paired_left_opt, paired_right_opt);
    paired_left_opt->needs(paired_right_opt);
    paired_right_opt->needs(paired_left_opt);

    ogroup->require_option(1, 2);

    app.add_option("-t,--threads", nthread,
                   "An integer that specifies the number of threads to use")
        ->default_val(16);
    app.add_flag("--quiet", quiet, "try to be quiet in terms of console output");

    CLI11_PARSE(app, argc, argv);

    auto input_filename = index_basename;
    auto read_filename = single_read_filenames;

    spdlog::drop_all();
    auto logger = spdlog::create<spdlog::sinks::stderr_color_sink_mt>("");
    logger->set_pattern("%+");
    if (quiet) { logger->set_level(spdlog::level::warn); }
    spdlog::set_default_logger(logger);

    mindex::reference_index ri(input_filename);

    std::string cmdline;
    size_t narg = static_cast<size_t>(argc);
    for (size_t i = 0; i < narg; ++i) {
        cmdline += std::string(argv[i]);
        cmdline.push_back(' ');
    }
    cmdline.pop_back();
    print_header(ri, cmdline);

    std::mutex iomut;

    // if we have paired-end data
    if (read_opt->empty()) {
        std::vector<std::thread> workers;
        uint32_t np = 1;
        if ((left_read_filenames.size() > 1) and (nthread >= 6)) {
            np += 1;
            nthread -= 1;
        }

        fastx_parser::FastxParser<fastx_parser::ReadPair> rparser(
            left_read_filenames, right_read_filenames, nthread, np);
        rparser.start();

        std::atomic<uint64_t> global_nr{0};
        std::atomic<uint64_t> global_nh{0};
        for (size_t i = 0; i < nthread; ++i) {
            workers.push_back(std::thread([&ri, &rparser, &global_nr, &global_nh, &iomut]() {
                do_map(ri, rparser, global_nr, global_nh, iomut);
            }));
        }

        for (auto& w : workers) { w.join(); }
        rparser.stop();
    } else {  // single-end
        std::vector<std::thread> workers;
        uint32_t np = 1;
        if ((single_read_filenames.size() > 1) and (nthread >= 6)) {
            np += 1;
            nthread -= 1;
        }

        fastx_parser::FastxParser<fastx_parser::ReadSeq> rparser(single_read_filenames, nthread,
                                                                 np);
        rparser.start();

        std::atomic<uint64_t> global_nr{0};
        std::atomic<uint64_t> global_nh{0};
        for (size_t i = 0; i < nthread; ++i) {
            workers.push_back(std::thread([&ri, &rparser, &global_nr, &global_nh, &iomut]() {
                do_map(ri, rparser, global_nr, global_nh, iomut);
            }));
        }

        for (auto& w : workers) { w.join(); }
        rparser.stop();
    }

    return 0;
}
