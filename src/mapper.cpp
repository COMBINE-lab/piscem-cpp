
#include "../external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "../include/reference_index.hpp"
#include "../include/kseq++.hpp"
#include "../include/CanonicalKmerIterator.hpp"
//#include "../include/query/contig_info_query_canonical_parsing.cpp"
#include "../include/projected_hits.hpp"
#include "../include/util.hpp"
#include "../include/mapping_util.hpp"
#include "../include/parallel_hashmap/phmap.h"
#include "hit_searcher.cpp"
#include "zlib.h"

#include <iostream>
#include <vector>
#include <numeric>

using namespace klibpp;



void print_header(mindex::reference_index& ri, std::string& cmdline) {
    std::cout << "@HD\tVN:1.0\tSO:unsorted\n";
    for (uint64_t i = 0; i < ri.num_refs(); ++i) {
        std::cout << "@SQ\tNS:" << ri.ref_name(i) << "\tLN:" << ri.ref_len(i) << "\n";
    }
    std::cout << "@PG\tID:mindex_map\tPN:mapper\tVN:0.0.1\t"
              << "CL:" << cmdline << "\n";
}

void do_map(mindex::reference_index& ri, const std::string& reads_filename) {
    CanonicalKmer::k(ri.k());

    KSeq record;
    gzFile fp = gzopen(reads_filename.c_str(), "r");
    auto ks = make_kstream(fp, gzread, mode::in);

    constexpr uint16_t is_secondary = 256;
    constexpr uint16_t is_unmapped = 4;
    constexpr uint16_t is_rc = 16;
    constexpr uint16_t first_seg = 64;

    // map from reference id to hit info
    phmap::flat_hash_map<uint32_t, mapping::util::sketch_hit_info> hit_map; 
    std::vector<mapping::util::simple_hit> accepted_hits;

    size_t max_occ_default = 200;
    size_t max_occ_recover = 1000;
    const bool attempt_occ_recover = (max_occ_recover > max_occ_default);
    // size_t alt_max_occ = 0;

    size_t num_hits = 0;
    sshash::contig_info_query_canonical_parsing q(ri.get_dict());
    mindex::hit_searcher hs(&ri);
    uint64_t read_num = 0;
    std::string workstr;
    int32_t k = static_cast<int32_t>(ri.k());
    while (ks >> record) {
        ++read_num;
        // std::cout << "readnum : " << read_num << "\n";
        if (read_num % 100000 == 0) {
            std::cerr << "readnum : " << read_num << ", num_hits : " << num_hits << "\n";
        }
        
        //alt_max_occ = 0;
        q.start();
        hs.clear();
        hit_map.clear();
        accepted_hits.clear();
        bool had_left_hit = hs.get_raw_hits_sketch(record.seq, q, true, false);
        bool early_stop = false;

        // if there were hits
        if (had_left_hit) {
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
            int32_t max_stretch = static_cast<int32_t>(record.seq.length() * 1.0);

            // a raw hit is a pair of read_pos and a projected hit

            // the least frequent hit for this fragment.
            uint64_t min_occ = std::numeric_limits<uint64_t>::max();

            // this is false by default and will be set to true
            // if *every* collected hit for this fragment occurs
            // max_occ_default times or more.
            bool had_alt_max_occ = false;
            int32_t signed_rl = static_cast<int32_t>(record.seq.length());
            auto collect_mappings_from_hits =
                [&max_stretch, &min_occ, &hit_map, &num_valid_hits, &total_occs,
                 &largest_occ, &early_stop, signed_rl, k](auto& raw_hits, auto& prev_read_pos, auto& max_allowed_occ,
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

                        for (auto& pos_it : refs) {
                            const auto& ref_pos_ori = proj_hits.decode_hit(pos_it);
                            uint32_t tid = pos_it.transcript_id();
                            int32_t pos = static_cast<int32_t>(ref_pos_ori.pos);
                            bool ori = ref_pos_ori.isFW;
                            auto& target = hit_map[tid];

                            // Why >= here instead of == ?
                            // Because hits can happen on the same target in both the forward
                            // and rc orientations, it is possible that we start the loop with
                            // the target having num_valid_hits hits in a given orientation (o)
                            // we see a new hit for this target in oriention o (now it has
                            // num_valid_hits + 1) then we see a hit for this target in orientation
                            // rc(o).  We still want to add / consider this hit, but
                            // max_hits_for_target() > num_valid_hits. So, we must allow for that
                            // here.
                            if (target.max_hits_for_target() >= num_valid_hits) {
                                if (ori) {
                                    target.add_fw(pos, static_cast<int32_t>(read_pos), 
                                                  signed_rl, k, max_stretch, score_inc);
                                } else {
                                    target.add_rc(pos, static_cast<int32_t>(read_pos), 
                                                  signed_rl, k, max_stretch, score_inc);
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
            auto mao_first_pass = max_occ_default - 1;
            early_stop =
                collect_mappings_from_hits(raw_hits, prev_read_pos, mao_first_pass, _discard);

            // If our default threshold was too stringent, then fallback to a more liberal
            // threshold and look up the k-mers that occur the least frequently.
            // Specifically, if the min occuring hits have frequency < max_occ_recover (2500 by
            // default) times, then collect the min occuring hits to get the mapping.
            if (attempt_occ_recover and (min_occ >= max_occ_default) and
                (min_occ < max_occ_recover)) {
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
                auto simple_hit = (best_hit_dir != mapping::util::HitDirection::RC) ? kv.second.get_fw_hit()
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
            if (accepted_hits.empty() and (num_valid_hits > 1) and (best_alt_hits >= num_valid_hits
            - 1)) { for (auto& kv : hit_map) { auto simple_hit = kv.second.get_best_hit(); if
            (simple_hit.num_hits >= best_alt_hits) {
                  //if (simple_hit.valid_pos(signed_read_len, transcripts[kv.first].RefLength, 10))
            { simple_hit.tid = kv.first; accepted_hits.emplace_back(simple_hit);
                  //}
                }
              }
            }
            */
        }  // DONE : if (rh)

        if (!accepted_hits.empty()) {
            bool secondary = false;
            for (auto& ah : accepted_hits) {
                uint16_t flag = secondary ? is_secondary : 0;
                //flag += 2;
                flag += ah.is_fw ? 0 : is_rc;
                //flag += first_seg;

                std::string* sptr = nullptr;
                if (is_rc) {
                    combinelib::kmers::reverseComplement(record.seq, workstr);
                    sptr = &workstr;
                } else {
                    sptr = &record.seq;
                }
                std::cout << record.name 
                    << "\t"
                    << flag 
                    << "\t" << ri.ref_name(ah.tid)
                    << "\t" << ah.pos + 1 
                    << "\t255\t*\t*\t0\t" << record.seq.length() 
                    << "\t" << *sptr << "\t*\n";
                secondary = true;
            }
        } else {
            std::cout << record.name  << "\t"
                << 4 << "\t"
                << "*\t0\t0\t*\t*\t0\t0\t"
                << record.seq << "\t*\n";
        }
    }
}

int main(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);

    /* mandatory arguments */
    parser.add("input_filename",
               "input index prefix.");
    parser.add("reads",
               "read filename.");

    if (!parser.parse()) return 1;

    auto input_filename = parser.get<std::string>("input_filename");
    auto read_filename = parser.get<std::string>("reads");

    mindex::reference_index ri(input_filename);

    std::string cmdline;
    size_t narg = static_cast<size_t>(argc);
    for (size_t i = 0; i < narg; ++i) {
        cmdline += std::string(argv[i]);
        cmdline.push_back(' ');
    }
    cmdline.pop_back();
    print_header(ri, cmdline);

    do_map(ri, read_filename);
    return 0;
}