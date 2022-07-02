#include "../external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "../include/reference_index.hpp"
#include "../include/CanonicalKmerIterator.hpp"
#include "../include/Kmer.hpp"
//#include "../include/query/contig_info_query_canonical_parsing.cpp"
#include "../include/query/streaming_query_canonical_parsing.hpp"
#include "../include/projected_hits.hpp"
#include "../include/util.hpp"
#include "../include/mapping_util.hpp"
#include "../include/parallel_hashmap/phmap.h"
#include "../include/FastxParser.hpp"
#include "../include/rad/rad_writer.hpp"
#include "../include/rad/rad_header.hpp"
#include "../include/ghc/filesystem.hpp"
#include "../include/cli11/CLI11.hpp"
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

size_t write_rad_header(mindex::reference_index& ri, std::ofstream& rad_file) {
    constexpr size_t bc_length = 16;
    constexpr size_t umi_length = 12;

    rad_writer bw;
    //  RADHeader
    rad_header rh;
    // right now, all formats are effectively single-end
    rh.is_paired(false);
    for (size_t i = 0; i < ri.num_refs(); ++i) { rh.add_refname(ri.ref_name(i)); }
    rh.dump_to_bin(bw);

    // where we will write the number of chunks when we know
    // how many there are
    size_t chunk_offset = bw.num_bytes() - sizeof(uint64_t);

    // ### start of tags

    // Tags we will have
    // write the tag meta-information section

    // File-level tag description
    uint16_t file_level_tags{2};
    bw << file_level_tags;

    // cblen
    uint8_t type_id{2};
    bw << std::string("cblen");
    bw << type_id;

    bw << std::string("ulen");
    bw << type_id;

    // read-level tag description
    uint16_t read_level_tags{2};
    bw << read_level_tags;

    // barcode
    bw << std::string("b");
    if (bc_length > 32) {
        type_id = 8;
    } else if (bc_length > 16) {
        // 17 - 32 bases
        type_id = 4;
    } else {
        // <= 16 bases
        type_id = 3;
    }
    bw << type_id;

    // umi
    bw << std::string("u");
    if (umi_length > 32) {
        type_id = 8;
    } else if (umi_length > 16) {
        // 17 - 32 bases
        type_id = 4;
    } else {
        // <= 16 bases
        type_id = 3;
    }
    bw << type_id;

    // alignment-level tag description
    uint16_t aln_level_tags{1};
    bw << aln_level_tags;
    // we maintain orientation
    // bw << std::string("orientation");
    // type_id = 1;
    // bw << type_id;

    // and reference id
    bw << std::string("compressed_ori_refid");
    type_id = 3;
    bw << type_id;

    // ### end of tag definitions

    // the actual file-level tags
    bw << static_cast<uint16_t>(16);
    bw << static_cast<uint16_t>(12);

    rad_file << bw;
    bw.clear();
    return chunk_offset;
}

// hardcode 10x v3 for now
class sc_protocol {
public:
    // We'd really like an std::optional<string&> here, but C++17
    // said no to that.
    std::string* extract_bc(std::string& r1, std::string& r2) {
        return (r1.length() >= bc_len) ? (bc.assign(r1, 0, bc_len), &bc) : nullptr;
    }

    // We'd really like an std::optional<string&> here, but C++17
    // said no to that.
    std::string* extract_umi(std::string& r1, std::string& r2) {
        return (r1.length() >= (bc_len + umi_len)) ? (umi.assign(r1, bc_len, umi_len), &umi)
                                                   : nullptr;
    }

    std::string* extract_mappable_read(std::string& r1, std::string& r2) { return &r2; }

private:
    std::string umi;
    std::string bc;
    const size_t bc_len = 16;
    const size_t umi_len = 12;
};

using umi_kmer_t = combinelib::kmers::Kmer<31, 2>;
using bc_kmer_t = combinelib::kmers::Kmer<31, 3>;

class output_info {
public:
    std::atomic<size_t> num_chunks{0};
    std::ofstream rad_file;
    std::mutex rad_mutex;

    std::ofstream unmapped_bc_file;
    std::mutex unmapped_bc_mutex;
};

enum class BarCodeRecovered : uint8_t { OK, RECOVERED, NOT_RECOVERED };

BarCodeRecovered recover_barcode(std::string& sequence) {
    size_t pos = sequence.find_first_not_of("ACTGactg");
    if (pos == std::string::npos) { return BarCodeRecovered::OK; }

    // Randomly assigning 'A' to first base with 'N'
    sequence[pos] = 'A';
    size_t invalid_pos = sequence.find_first_not_of("ACTGactg", pos);
    return (invalid_pos == std::string::npos) ? BarCodeRecovered::RECOVERED
                                              : BarCodeRecovered::NOT_RECOVERED;
}

inline void write_to_rad_stream(bc_kmer_t& bck, umi_kmer_t& umi,
                                mapping::util::MappingType map_type,
                                std::vector<mapping::util::simple_hit>& accepted_hits,
                                mindex::reference_index& ri,
                                phmap::flat_hash_map<uint64_t, uint32_t>& unmapped_bc_map,
                                uint32_t& num_reads_in_chunk, rad_writer& bw) {
    constexpr uint32_t barcode_len = 16;
    constexpr uint32_t umi_len = 12;

    if (map_type == mapping::util::MappingType::SINGLE_MAPPED) {
        // number of mappings
        bw << static_cast<uint32_t>(accepted_hits.size());

        // bc
        // if we can fit the barcode into an integer
        if (barcode_len <= 32) {
            if (barcode_len <= 16) {  // can use 32-bit int
                uint32_t shortbck = static_cast<uint32_t>(0x00000000FFFFFFFF & bck.word(0));
                bw << shortbck;
            } else {  // must use 64-bit int
                bw << bck.word(0);
            }
        } else {  // must use a string for the barcode
            std::cerr << "should not happen\n";
        }

        // umi
        if (umi_len <= 16) {  // if we can use 32-bit int
            uint64_t umiint = umi.word(0);
            uint32_t shortumi = static_cast<uint32_t>(0x00000000FFFFFFFF & umiint);
            bw << shortumi;
        } else if (umi_len <= 32) {  // if we can use 64-bit int
            uint64_t umiint = umi.word(0);
            bw << umiint;
        } else {  // must use string
            std::cerr << "should not happen\n";
        }

        for (auto& aln : accepted_hits) {
            uint32_t fw_mask = aln.is_fw ? 0x80000000 : 0x00000000;
            bw << (aln.tid | fw_mask);
        }
        ++num_reads_in_chunk;
    } else {
        unmapped_bc_map[bck.word(0)] += 1;
    }
}

inline void write_to_sam_stream(fastx_parser::ReadSeq& record,
                                std::vector<mapping::util::simple_hit>& accepted_hits,
                                mindex::reference_index& ri, std::string& workstr,
                                std::ostringstream& osstream) {
    constexpr uint16_t is_secondary = 256;
    constexpr uint16_t is_rc = 16;

    if (!accepted_hits.empty()) {
        bool secondary = false;
        size_t nh = accepted_hits.size();
        size_t hit_index = 0;
        for (auto& ah : accepted_hits) {
            uint16_t flag = secondary ? is_secondary : 0;
            // flag += 2;
            flag += ah.is_fw ? 0 : is_rc;
            // flag += first_seg;

            std::string* sptr = nullptr;
            if (is_rc) {
                combinelib::kmers::reverseComplement(record.seq, workstr);
                sptr = &workstr;
            } else {
                sptr = &record.seq;
            }
            osstream << record.name << "\t" << flag << "\t" << ri.ref_name(ah.tid) << "\t"
                     << ah.pos + 1 << "\t255\t*\t*\t0\t" << record.seq.length() << "\t" << *sptr
                     << "\t*\tNH:i:" << nh << "\tHI:i:" << hit_index << "\n";
            secondary = true;
            hit_index++;
        }
    } else {
        osstream << record.name << "\t" << 4 << "\t"
                 << "*\t0\t0\t*\t*\t0\t0\t" << record.seq << "\t*\tNH:i:0\n";
    }
}

void do_map(mindex::reference_index& ri, fastx_parser::FastxParser<fastx_parser::ReadPair>& parser,
            std::atomic<uint64_t>& global_nr, std::atomic<uint64_t>& global_nhits,
            output_info& out_info, std::mutex& iomut) {
    CanonicalKmer::k(ri.k());

    // put these in struct
    size_t num_short_umi{0};
    size_t num_ambig_umi{0};

    // map from reference id to hit info
    phmap::flat_hash_map<uint32_t, mapping::util::sketch_hit_info> hit_map;
    std::vector<mapping::util::simple_hit> accepted_hits;
    mapping::util::MappingType map_type = mapping::util::MappingType::UNMAPPED;

    // map to recall the number of unmapped reads we see
    // for each barcode
    phmap::flat_hash_map<uint64_t, uint32_t> unmapped_bc_map;

    size_t max_occ_default = 200;
    size_t max_occ_recover = 1000;
    const bool attempt_occ_recover = (max_occ_recover > max_occ_default);
    size_t alt_max_occ = 2500;

    // sshash::contig_info_query_canonical_parsing q(ri.get_dict());
    sshash::streaming_query_canonical_parsing q(ri.get_dict());
    mindex::hit_searcher hs(&ri);
    uint64_t read_num = 0;
    uint64_t processed = 0;
    uint64_t buff_size = 10000;
    size_t max_chunk_reads = 5000;
    std::string workstr;

    std::ostringstream osstream;

    int32_t k = static_cast<int32_t>(ri.k());
    // Get the read group by which this thread will
    // communicate with the parser (*once per-thread*)
    auto rg = parser.getReadGroup();

    sc_protocol protocol;
    rad_writer rad_w;
    // reserve the space to later write
    // down the number of reads in the
    // first chunk.
    uint32_t num_reads_in_chunk{0};
    rad_w << num_reads_in_chunk;
    rad_w << num_reads_in_chunk;

    while (parser.refill(rg)) {
        // Here, rg will contain a chunk of read pairs
        // we can process.
        for (auto& record : rg) {
            ++global_nr;
            ++read_num;
            auto rctr = global_nr.load();
            auto hctr = global_nhits.load();
            if (rctr % 500000 == 0) {
              iomut.lock();
                std::cerr << "\rreadnum : " << rctr << ", num_hits : " << hctr ;
                iomut.unlock();
            }

            // first extract the barcode
            std::string* bc = protocol.extract_bc(record.first.seq, record.second.seq);
            // if we couldn't get it, don't bother with
            // anything else for this read.
            if (bc == nullptr) { continue; }

            // correct up to one `N` in the barcode
            auto recovered = recover_barcode(*bc);
            // if we couldn't correct it with 1 `N`, then skip.
            if (recovered == BarCodeRecovered::NOT_RECOVERED) { continue; }

            // convert it to a k-mer type
            bc_kmer_t bc_kmer;
            bool bc_ok = bc_kmer.fromChars(*bc);
            if (!bc_ok) { continue; }

            std::string* umi = protocol.extract_umi(record.first.seq, record.second.seq);
            if (umi == nullptr) {
                num_short_umi++;
                continue;
            }

            // convert it to a k-mer type
            umi_kmer_t umi_kmer;
            bool umi_ok = umi_kmer.fromChars(*umi);
            if (!umi_ok) {
                num_ambig_umi++;
                continue;
            }

            // alt_max_occ = 0;
            std::string* read_seq =
                protocol.extract_mappable_read(record.first.seq, record.second.seq);
            map_type = mapping::util::MappingType::UNMAPPED;
            q.start();
            hs.clear();
            hit_map.clear();
            accepted_hits.clear();
            bool had_left_hit = hs.get_raw_hits_sketch(*read_seq, q, true, false);
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
                     &early_stop, signed_rl,
                     k](auto& raw_hits, auto& prev_read_pos, auto& max_allowed_occ,
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
                    early_stop = collect_mappings_from_hits(raw_hits, prev_read_pos,
                                                            max_allowed_occ, had_alt_max_occ);
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
                        best_alt_hits = simple_hit.num_hits > best_alt_hits ? simple_hit.num_hits
                                                                            : best_alt_hits;
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
            if (accepted_hits.size() > alt_max_occ) {
                accepted_hits.clear();
                map_type = mapping::util::MappingType::UNMAPPED;
                std::cerr << "boo\n";
            } else if (!accepted_hits.empty()) {
                map_type = mapping::util::MappingType::SINGLE_MAPPED;
            }

            global_nhits += accepted_hits.empty() ? 0 : 1;
            write_to_rad_stream(bc_kmer, umi_kmer, map_type, accepted_hits, ri, unmapped_bc_map,
                                num_reads_in_chunk, rad_w);

            // dump buffer
            if (num_reads_in_chunk > max_chunk_reads) {
                out_info.num_chunks++;
                uint32_t num_bytes = rad_w.num_bytes();
                rad_w.write_integer_at_offset(0, num_bytes);
                rad_w.write_integer_at_offset(sizeof(num_bytes), num_reads_in_chunk);
                out_info.rad_mutex.lock();
                out_info.rad_file << rad_w;
                out_info.rad_mutex.unlock();
                rad_w.clear();
                num_reads_in_chunk = 0;

                // reserve space for headers of next chunk
                rad_w << num_reads_in_chunk;
                rad_w << num_reads_in_chunk;
            }
        }
    }

    // dump any remaining output
    if (num_reads_in_chunk > 0) {
        out_info.num_chunks++;
        uint32_t num_bytes = rad_w.num_bytes();
        rad_w.write_integer_at_offset(0, num_bytes);
        rad_w.write_integer_at_offset(sizeof(num_bytes), num_reads_in_chunk);
        out_info.rad_mutex.lock();
        out_info.rad_file << rad_w;
        out_info.rad_mutex.unlock();
        rad_w.clear();
        num_reads_in_chunk = 0;
    }

    // unmapped barcode writer
    {  // make a scope and dump the unmapped barcode counts
        rad_writer ubcw;
        for (auto& kv : unmapped_bc_map) {
            ubcw << kv.first;
            ubcw << kv.second;
        }
        out_info.unmapped_bc_mutex.lock();
        out_info.unmapped_bc_file << ubcw;
        out_info.unmapped_bc_mutex.unlock();
        ubcw.clear();
    }
}

int main(int argc, char** argv) {
    /**
     * PESC : Pseudoalignment Enhanced with Structural Constraints
     **/
    std::ios_base::sync_with_stdio(false);

    std::string input_filename;
    std::vector<std::string> filenames1;
    std::vector<std::string> filenames2;
    std::string output_dirname;
    size_t nthread;
    CLI::App app{"PESC — single-cell RNA-seq mapper for alevin-fry"};
    app.add_option("-i,--index", input_filename, "input index prefix")->required();
    app.add_option("-1,--read1", filenames1, "path to list of read 1 files")->required()->delimiter(',');
    app.add_option("-2,--read2", filenames2, "path to list of read 2 files")->required()->delimiter(',');
    app.add_option("-o,--output", output_dirname, "path to output directory")->required();
    app.add_option("-t,--threads", nthread, "An integer that specifies the number of threads to use")->default_val(16);
    CLI11_PARSE(app, argc, argv);

    /*
    cmd_line_parser::parser parser(argc, argv);

    size_t nthread = 16;
    // mandatory arguments
    parser.add("input_filename", "input index prefix.");
    parser.add("read1", "read 1 filename.");
    parser.add("read2", "read 2 filename.");
    parser.add("output", "output directory.");
    parser.add("t", "A (integer) that specifies the number of threads to use.", "-t", false);
    // parser.add("reads", "read filename.");
    if (!parser.parse()) return 1;
    */

    umi_kmer_t::k(12);
    bc_kmer_t::k(16);

    // RAD file path
    ghc::filesystem::path output_path(output_dirname);
    ghc::filesystem::create_directories(output_path);

    ghc::filesystem::path rad_file_path = output_path / "map.rad";
    ghc::filesystem::path unmapped_bc_file_path = output_path / "unmapped_bc_count.bin";

    std::cerr << "thread: " << nthread << "\n";

    mindex::reference_index ri(input_filename);

    std::string cmdline;
    size_t narg = static_cast<size_t>(argc);
    for (size_t i = 0; i < narg; ++i) {
        cmdline += std::string(argv[i]);
        cmdline.push_back(' ');
    }
    cmdline.pop_back();
    // print_header(ri, cmdline);

    std::ofstream rad_file(rad_file_path.string());
    std::ofstream unmapped_bc_file(unmapped_bc_file_path.string());

    if (!rad_file.good()) {
        // alevinOpts.jointLog->error("Could not open {} for writing.", rad_file_path.string());
        throw std::runtime_error("error creating output file.");
    }

    if (!unmapped_bc_file.good()) {
        // alevinOpts.jointLog->error("Could not open {} for writing.",
        // unmapped_bc_file_path.string());
        throw std::runtime_error("error creating output file.");
    }
    output_info out_info;
    out_info.rad_file = std::move(rad_file);
    out_info.unmapped_bc_file = std::move(unmapped_bc_file);

    size_t chunk_offset = write_rad_header(ri, out_info.rad_file);

    std::atomic<uint64_t> num_chunks{0};
    std::mutex iomut;

    uint32_t np = 1;
    if ((filenames1.size() > 1) and (nthread >= 6)) {
        np += 1;
        nthread -= 1;
    }

    fastx_parser::FastxParser<fastx_parser::ReadPair> rparser(filenames1, filenames2, nthread, np);
    rparser.start();

    std::atomic<uint64_t> global_nr{0};
    std::atomic<uint64_t> global_nh{0};
    std::vector<std::thread> workers;
    for (size_t i = 0; i < nthread; ++i) {
        workers.push_back(std::thread([&ri, &rparser, &global_nr, &global_nh, &iomut, &out_info]() {
            do_map(ri, rparser, global_nr, global_nh, out_info, iomut);
        }));
    }

    for (auto& w : workers) { w.join(); }
    rparser.stop();
    std::cerr << "\n";

    out_info.rad_file.seekp(chunk_offset);
    uint64_t nc = out_info.num_chunks.load();
    out_info.rad_file.write(reinterpret_cast<char*>(&nc), sizeof(nc));

    out_info.rad_file.close();

    // We want to check if the RAD file stream was written to properly.
    // While we likely would have caught this earlier, it is possible the
    // badbit may not be set until the stream actually flushes (perhaps even
    // at close), so we check here one final time that the status of the
    // stream is as expected.
    // see :
    // https://stackoverflow.com/questions/28342660/error-handling-in-stdofstream-while-writing-data
    if (!out_info.rad_file) {
        std::cerr
            << "The RAD file stream had an invalid status after close; so some operation(s) may"
            << "have failed!\nA common cause for this is lack of output disk space.\n"
            << "Consequently, the output may be corrupted.\n\n";
        return 1;
    }

    out_info.unmapped_bc_file.close();

    return 0;
}