#include "../include/reference_index.hpp"
#include "../include/util.hpp"
#include "../include/mapping/utils.hpp"
#include "../include/spdlog_piscem/spdlog.h"
#include "../include/spdlog_piscem/sinks/stdout_color_sinks.h"
#include "../include/rad/rad_writer.hpp"
#include "../include/rad/rad_header.hpp"
#include "../include/rad/util.hpp"
#include "../include/ghc/filesystem.hpp"
#include "../include/cli11/CLI11.hpp"
#include "../include/meta_info.hpp"
#include "../include/FastxParser.hpp"
#include "../include/poison_table.hpp"
#include "zlib.h"
//#include "FastxParser.cpp"
//#include "../src/hit_searcher.cpp"

#include <atomic>
#include <iostream>
#include <vector>
#include <numeric>
#include <cstdio>
#include <thread>
#include <sstream>
#include <type_traits>

using namespace klibpp;
using mapping::util::mapping_cache_info;
using mapping::util::poison_state_t;

//using poison_map_t = phmap::flat_hash_set<uint64_t, sshash::RobinHoodHash>;

// utility class that wraps the information we will
// need access to when writing output within each thread
// as well as information we'll need to update for the
// caller.
class mapping_output_info {
public:
    // will keep track of the total number
    // of chunks written to rad_file
    std::atomic<size_t> num_chunks{0};
    // the output stream where the actual
    // RAD records are written
    std::ofstream rad_file;
    // the mutex for safely writing to
    // rad_file
    std::mutex rad_mutex;
    // will record the total number of observed fragments
    std::atomic<size_t> observed_fragments{0};
};

void print_header(mindex::reference_index& ri, std::string& cmdline) {
    std::cout << "@HD\tVN:1.0\tSO:unsorted\n";
    for (uint64_t i = 0; i < ri.num_refs(); ++i) {
        std::cout << "@SQ\tSN:" << ri.ref_name(i) << "\tLN:" << ri.ref_len(i) << "\n";
    }
    std::cout << "@PG\tID:mindex_map\tPN:mapper\tVN:0.0.1\t"
              << "CL:" << cmdline << "\n";
}

// single-end
bool map_fragment(fastx_parser::ReadSeq& record, poison_state_t& poison_state, mapping_cache_info& map_cache_left,
                  mapping_cache_info& map_cache_right, mapping_cache_info& map_cache_out) {
    (void)map_cache_left;
    (void)map_cache_right;
    /*
    if (poison_state.is_poisoned()) {
      map_cache_out.clear();
      return false;
    } else {
    }
    */
    return mapping::util::map_read(&record.seq, map_cache_out, poison_state, false);
}

// paried-end
bool map_fragment(fastx_parser::ReadPair& record, poison_state_t& poison_state, mapping_cache_info& map_cache_left,
                  mapping_cache_info& map_cache_right, mapping_cache_info& map_cache_out) {
    // don't map a poisned read pair
    /*
    if (poison_state.is_poisoned()) {
      map_cache_out.clear();
      return false;
    }
    */

    bool early_exit_left = mapping::util::map_read(&record.first.seq, map_cache_left, poison_state, false);
    if (poison_state.is_poisoned()) { return false; }
    bool early_exit_right = mapping::util::map_read(&record.second.seq, map_cache_right, poison_state, false);

    int32_t left_len = static_cast<int32_t>(record.first.seq.length());
    int32_t right_len = static_cast<int32_t>(record.second.seq.length());
    
    /*
    for (auto& lh : map_cache_left.accepted_hits) {
      std::cerr << "left: " << lh.tid << ", " << lh.pos << " (" << (lh.is_fw ? "fw" : "rc") << ")\n"; 
    }
    for (auto& lh : map_cache_right.accepted_hits) {
      std::cerr << "right: " << lh.tid << ", " << lh.pos << " (" << (lh.is_fw ? "fw" : "rc") << ")\n"; 
    }
    */

    mapping::util::merge_se_mappings(map_cache_left, map_cache_right, left_len, right_len,
                                     map_cache_out);

    return (early_exit_left or early_exit_right);
}

// single-end
inline void write_sam_mappings(mapping_cache_info& map_cache_out, fastx_parser::ReadSeq& record,
                               std::string& workstr_left, std::string& workstr_right,
                               std::atomic<uint64_t>& global_nhits, std::ostringstream& osstream) {
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
            osstream << record.name << "\t" << flag << "\t"
                     << map_cache_out.hs.get_index()->ref_name(ah.tid) << "\t" << ah.pos + 1
                     << "\t255\t*\t*\t0\t" << record.seq.length() << "\t" << *sptr << "\t*\n";
            secondary = true;
        }
    } else {
        osstream << record.name << "\t" << 4 << "\t"
                 << "*\t0\t0\t*\t*\t0\t0\t" << record.seq << "\t*\n";
    }
}

// paired-end
inline void write_sam_mappings(mapping_cache_info& map_cache_out, fastx_parser::ReadPair& record,
                               std::string& workstr_left, std::string& workstr_right,
                               std::atomic<uint64_t>& global_nhits, std::ostringstream& osstream) {
    (void)workstr_right;
    constexpr uint16_t is_secondary = 256;
    constexpr uint16_t is_rc = 16;
    constexpr uint16_t mate_rc = 32;
    constexpr uint16_t unmapped = 4;
    constexpr uint16_t mate_unmapped = 8;

    auto map_type = map_cache_out.map_type;

    if (!map_cache_out.accepted_hits.empty()) {
        ++global_nhits;
        bool secondary = false;

        uint16_t base_flag_first = 0;
        uint16_t base_flag_second = 0;
        switch (map_type) {
            case mapping::util::MappingType::MAPPED_FIRST_ORPHAN:
                base_flag_first = 73;
                base_flag_second = 133;
                break;
            case mapping::util::MappingType::MAPPED_SECOND_ORPHAN:
                base_flag_first = 69;
                base_flag_second = 137;
                break;
            case mapping::util::MappingType::MAPPED_PAIR:
                base_flag_first = 67;
                base_flag_second = 131;
                break;
            default:
                break;
        }

        bool have_rc_first = false;
        bool have_rc_second = false;
        for (auto& ah : map_cache_out.accepted_hits) {
            uint16_t flag_first = secondary ? base_flag_first + is_secondary : base_flag_first;
            uint16_t flag_second = secondary ? base_flag_second + is_secondary : base_flag_second;

            std::string* sptr_first = nullptr;
            std::string* sptr_second = nullptr;
            int32_t pos_first = 0;
            int32_t pos_second = 0;

            // if both reads are mapped
            if (map_type == mapping::util::MappingType::MAPPED_PAIR) {
                pos_first = ah.pos + 1;
                pos_second = ah.mate_pos + 1;

                if (ah.is_fw) {
                    flag_first += mate_rc;
                    sptr_first = &record.first.seq;

                    flag_second += is_rc;
                    if (!have_rc_second) {
                        have_rc_second = true;
                        combinelib::kmers::reverseComplement(record.second.seq, workstr_right);
                    }
                    sptr_second = &workstr_right;
                } else {
                    flag_first += is_rc;
                    if (!have_rc_first) {
                        have_rc_first = true;
                        combinelib::kmers::reverseComplement(record.first.seq, workstr_left);
                    }
                    sptr_first = &workstr_left;

                    flag_second += mate_rc;
                    sptr_second = &record.second.seq;
                }
            } else if (map_type == mapping::util::MappingType::MAPPED_FIRST_ORPHAN) {
                pos_first = ah.pos;
                pos_second = 0;

                sptr_first = &record.first.seq;
                sptr_second = &record.second.seq;

                if (!ah.is_fw) {  // if the mapped read is rc
                    flag_first += is_rc;
                    if (!have_rc_first) {
                        have_rc_first = true;
                        combinelib::kmers::reverseComplement(record.first.seq, workstr_left);
                    }
                    sptr_first = &workstr_left;

                    flag_second += mate_rc;
                }

            } else if (map_type == mapping::util::MappingType::MAPPED_SECOND_ORPHAN) {
                pos_first = 0;
                pos_second = ah.pos + 1;

                sptr_first = &record.first.seq;
                sptr_second = &record.second.seq;
                if (!ah.is_fw) {
                    flag_first += mate_rc;
                    flag_second += is_rc;
                    if (!have_rc_second) {
                        have_rc_second = true;
                        combinelib::kmers::reverseComplement(record.second.seq, workstr_right);
                    }
                    sptr_second = &workstr_right;
                }
            }

            auto print_pos_mapq_cigar = [](bool mapped, int32_t pos, int32_t read_len,
                                           int32_t ref_len, std::ostream& os) {
                if (!mapped) {
                    os << "0\t255\t*\t";
                    return;
                }
                int32_t pad_start = 0;
                int32_t pad_end = 0;
                int32_t m_len = read_len;

                if (pos + read_len >= ref_len) {
                    pad_end = (pos + read_len) - ref_len + 1;
                    m_len -= pad_end;
                }

                if (pos <= 0) {
                    pad_start = (-pos) + 1;
                    m_len -= pad_start;
                    pos = 1;
                }

                os << pos << "\t255\t";

                if (pad_start > 0) { os << pad_start << "S"; }
                if (pad_end > 0) {
                    os << m_len << 'M' << pad_end << "S\t";
                } else {
                    os << m_len << "M\t";
                }
            };

            const auto ref_name = map_cache_out.hs.get_index()->ref_name(ah.tid);
            const int32_t ref_len =
                static_cast<int32_t>(map_cache_out.hs.get_index()->ref_len(ah.tid));

            osstream << record.first.name << "\t" << flag_first << "\t"
                     << ((flag_first & unmapped) ? "*" : ref_name)
                     << '\t';  // if mapped RNAME, else *
            print_pos_mapq_cigar(!(flag_first & unmapped), pos_first,
                                 static_cast<int32_t>(record.first.seq.length()), ref_len,
                                 osstream);
            osstream << ((flag_first & mate_unmapped) ? '*' : '=') << '\t'  // RNEXT
                     << ((flag_first & mate_unmapped) ? 0 : std::max(1, pos_second))
                     << '\t'  // PNEXT
                     << ah.frag_len() << '\t' << *sptr_first << "\t*\n";
            osstream << record.second.name << "\t" << flag_second << "\t"
                     << ((flag_second & unmapped) ? "*" : ref_name)
                     << '\t';  // if mapped RNAME, else *
            print_pos_mapq_cigar(!(flag_second & unmapped), pos_second,
                                 static_cast<int32_t>(record.second.seq.length()), ref_len,
                                 osstream);
            osstream << ((flag_second & mate_unmapped) ? '*' : '=') << '\t'  // RNEXT
                     << ((flag_second & mate_unmapped) ? 0 : std::max(1, pos_first))
                     << '\t'  // PNEXT
                     << -ah.frag_len() << '\t' << *sptr_second << "\t*\n";
            secondary = true;
        }
    } else {
        osstream << record.first.name << "\t" << 77 << "\t"
                 << "*\t0\t0\t*\t*\t0\t0\t" << record.first.seq << "\t*\n";
        osstream << record.second.name << "\t" << 141 << "\t"
                 << "*\t0\t0\t*\t*\t0\t0\t" << record.second.seq << "\t*\n";
    }
}

std::string& get_name(fastx_parser::ReadSeq& rs) { return rs.name; }

std::string& get_name(fastx_parser::ReadPair& rs) { return rs.first.name; }

template <typename FragT>
void do_map(mindex::reference_index& ri, fastx_parser::FastxParser<FragT>& parser,
            poison_table& poison_map,
            std::atomic<uint64_t>& global_npoisoned,
            std::atomic<uint64_t>& global_nr, std::atomic<uint64_t>& global_nhits,
            mapping_output_info& out_info, std::mutex& iomut) {
    auto log_level = spdlog_piscem::get_level();
    auto write_mapping_rate = false;
    switch (log_level) {
        case spdlog_piscem::level::level_enum::trace:
            write_mapping_rate = true;
            break;
        case spdlog_piscem::level::level_enum::debug:
            write_mapping_rate = true;
            break;
        case spdlog_piscem::level::level_enum::info:
            write_mapping_rate = true;
            break;
        case spdlog_piscem::level::level_enum::warn:
            write_mapping_rate = false;
            break;
        case spdlog_piscem::level::level_enum::err:
            write_mapping_rate = false;
            break;
        case spdlog_piscem::level::level_enum::critical:
            write_mapping_rate = false;
            break;
        case spdlog_piscem::level::level_enum::off:
            write_mapping_rate = false;
            break;
        default:
            write_mapping_rate = false;
    }
    
    bool use_poison = !poison_map.empty();
        
    pufferfish::CanonicalKmerIterator kit_end;
    auto pmap_end = poison_map.kmer_end();
    
    // checks if a read is "poisned", returns true if it is
    // and false otherwise.
    auto is_poisoned = [&](const std::string& seq) -> bool {
      pufferfish::CanonicalKmerIterator kit(seq);
      while (kit != kit_end) {
        // current canonical k-mer
        if (kit->first.is_homopolymer()) { 
          kit++;
          continue; 
        }
        

        auto pmap_it = poison_map.find_kmer(kit->first.getCanonicalWord());
        if (pmap_it != pmap_end) {
          return true;
        }
        ++kit;
      }
      return false;
    };

    mapping_cache_info map_cache_left(ri);
    mapping_cache_info map_cache_right(ri);
    mapping_cache_info map_cache_out(ri);
    poison_state_t poison_state;
    if (use_poison) {
      poison_state.ptab = &poison_map;
    }
    
    // the reads are paired
    if constexpr(std::is_same_v<fastx_parser::ReadPair, FragT>) {
      poison_state.paired_for_mapping = true;
    }

    rad_writer rad_w;
    size_t max_chunk_reads = 5000;
    // reserve the space to later write
    // down the number of reads in the
    // first chunk.
    uint32_t num_reads_in_chunk{0};
    rad_w << num_reads_in_chunk;
    rad_w << num_reads_in_chunk;

    sshash::streaming_query_canonical_parsing q(ri.get_dict());
    mindex::hit_searcher hs(&ri);
    uint64_t read_num = 0;
    // SAM output
    // uint64_t processed = 0;
    // uint64_t buff_size = 10000;

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

            if (write_mapping_rate and (rctr % 500000 == 0)) {
                iomut.lock();
                std::cerr << "\rprocessed (" << rctr << ") reads; (" << hctr << ") had mappings.";
                iomut.unlock();
            }
            
            poison_state.clear();
            /*
            if (use_poison) {
              if constexpr(std::is_same_v<fastx_parser::ReadSeq, FragT>) {
                  poison_state.poisoned_left = is_poisoned(record.seq);
              } else if constexpr(std::is_same_v<fastx_parser::ReadPair, FragT>) {
                   poison_state.poisoned_left = is_poisoned(record.first.seq);
                   poison_state.poisoned_right = is_poisoned(record.second.seq);
              }
            }
            */

            // this *overloaded* function will just do the right thing.
            // If record is single-end, just map that read, otherwise, map both and look
            // for proper pairs.
            bool had_early_stop =
                map_fragment(record, poison_state, map_cache_left, map_cache_right, map_cache_out);
            (void)had_early_stop;
            if (poison_state.is_poisoned()) {
              global_npoisoned++;
              //if (global_npoisned > 0 and global_npoisned % 100000 == 0) {
              // spdlog_piscem::info("number of mapped poisned reads {}", global_npoisned);
              //}
            }
            // to write unmapped names
            /*
            if (map_cache_out.accepted_hits.empty()) {
                iomut.lock();
                std::cout << get_name(record) << "\n";
                iomut.unlock();
            }
            */
            /*
            if constexpr( std::is_same_v<FragT, fastx_parser::ReadSeq> ) {
              if (map_cache_out.accepted_hits.empty()) {
                std::cout << ">" << record.name << "\n";
                std::cout << record.seq << "\n";
              }
            }
            */
            // RAD output
            global_nhits += map_cache_out.accepted_hits.empty() ? 0 : 1;
            rad::util::write_to_rad_stream_bulk(map_cache_out.map_type, map_cache_out.accepted_hits,
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

            // SAM output
            /*
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
            */
            
        }
    }

    // RAD output: dump any remaining output
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

    // SAM output
    // dump any remaining output
    /*
      std::string o = osstream.str();
      iomut.lock();
      std::cout << o;
      iomut.unlock();
      osstream.clear();
    */
    // don't need this here because osstream goes away at end of scope
    // osstream.str("");
    
}

#ifdef __cplusplus
extern "C" {
#endif
  int run_pesc_bulk(int argc, char** argv);
#ifdef __cplusplus
}
#endif

int run_pesc_bulk(int argc, char** argv) {
    /**
     * Mapper
     **/
    std::ios_base::sync_with_stdio(false);

    std::string index_basename;
    std::vector<std::string> left_read_filenames;
    std::vector<std::string> right_read_filenames;
    std::vector<std::string> single_read_filenames;
    std::string output_stem;
    size_t nthread{16};
    bool quiet{false};
    bool no_poison{false};

    CLI::App app{"Bulk mapper"};
    app.add_option("-i,--index", index_basename, "Input index prefix")->required();

    auto ogroup = app.add_option_group("input reads", "provide input reads");

    CLI::Option* read_opt =
        ogroup->add_option("-r,--reads", single_read_filenames, "Path to list (comma separated) of single-end files")
            ->delimiter(',');
    CLI::Option* paired_left_opt =
        ogroup->add_option("-1,--read1", left_read_filenames, "Path to list (comma separated) of read 1 files")
            ->delimiter(',');
    CLI::Option* paired_right_opt =
        ogroup->add_option("-2,--read2", right_read_filenames, "Path to list (comma separated) of read 2 files")
            ->delimiter(',');

    paired_left_opt->excludes(read_opt);
    paired_right_opt->excludes(read_opt);
    read_opt->excludes(paired_left_opt, paired_right_opt);
    paired_left_opt->needs(paired_right_opt);
    paired_right_opt->needs(paired_left_opt);

    ogroup->require_option(1, 2);

    app.add_option("-o,--output", output_stem, "The file stem where output should be written.")
        ->required();
    app.add_option("-t,--threads", nthread,
                   "An integer that specifies the number of threads to use.")
        ->default_val(16);
    app.add_flag("--no-poison", no_poison, "Do not filter reads for poison k-mers, even if a poison table exists for the index");
    app.add_flag("--quiet", quiet, "Try to be quiet in terms of console output");

    CLI11_PARSE(app, argc, argv);

    auto input_filename = index_basename;
    auto read_filename = single_read_filenames;

    spdlog_piscem::drop_all();
    auto logger = spdlog_piscem::create<spdlog_piscem::sinks::stderr_color_sink_mt>("");
    logger->set_pattern("%+");
    spdlog_piscem::set_default_logger(logger);

    if (quiet) {
      spdlog_piscem::set_level(spdlog_piscem::level::warn);
    }

    // start the timer
    auto start_t = std::chrono::high_resolution_clock::now();

    mindex::reference_index ri(input_filename);

    bool is_paired = read_opt->empty();

    std::string cmdline;
    size_t narg = static_cast<size_t>(argc);
    for (size_t i = 0; i < narg; ++i) {
        cmdline += std::string(argv[i]);
        cmdline.push_back(' ');
    }
    cmdline.pop_back();
    // print_header(ri, cmdline);

    ghc::filesystem::path rad_file_path = output_stem + ".rad";

    std::ofstream rad_file(rad_file_path.string());
    if (!rad_file.good()) {
        spdlog_piscem::critical("Could not open {} for writing.", rad_file_path.string());
        throw std::runtime_error("error creating output file.");
    }

    mapping_output_info out_info;
    out_info.rad_file = std::move(rad_file);
    size_t chunk_offset = rad::util::write_rad_header_bulk(ri, is_paired, out_info.rad_file);

    std::mutex iomut;

    // set the canonical k-mer size globally
    CanonicalKmer::k(ri.k());
    std::atomic<uint64_t> global_nr{0};
    std::atomic<uint64_t> global_nh{0};
    std::atomic<uint64_t> global_np{0};
    uint32_t np = 1;

    // load a poison map if we had one
    poison_table ptab;
    if (!no_poison and poison_table::exists(input_filename)) {
      /*
      spdlog_piscem::info("Loading poison k-mer map...");
      phmap::BinaryInputArchive ar_in(pmap_file.c_str());
      poison_map.phmap_load(ar_in);
      spdlog_piscem::info("done");
      */
      poison_table ptab_tmp(input_filename);
      ptab = std::move(ptab_tmp);
    } else {
      spdlog_piscem::info("No poison k-mer map exists, or it was requested not to be used");
    }

    // if we have paired-end data
    if (read_opt->empty()) {
        std::vector<std::thread> workers;
        if ((left_read_filenames.size() > 1) and (nthread >= 6)) {
            np += 1;
            nthread -= 1;
        }

        fastx_parser::FastxParser<fastx_parser::ReadPair> rparser(
            left_read_filenames, right_read_filenames, nthread, np);
        rparser.start();

        for (size_t i = 0; i < nthread; ++i) {
            workers.push_back(
                std::thread([&ri, &rparser, &ptab, &global_np, &global_nr, &global_nh, &out_info, &iomut]() {
                    do_map(ri, rparser, ptab, global_np, global_nr, global_nh, out_info, iomut);
                }));
        }

        for (auto& w : workers) { w.join(); }
        rparser.stop();
    } else {  // single-end
        std::vector<std::thread> workers;
        if ((single_read_filenames.size() > 1) and (nthread >= 6)) {
            np += 1;
            nthread -= 1;
        }

        fastx_parser::FastxParser<fastx_parser::ReadSeq> rparser(single_read_filenames, nthread,
                                                                 np);
        rparser.start();

        for (size_t i = 0; i < nthread; ++i) {
            workers.push_back(
                std::thread([&ri, &rparser, &ptab, &global_np, &global_nr, &global_nh, &out_info, &iomut]() {
                    do_map(ri, rparser, ptab, global_np, global_nr, global_nh, out_info, iomut);
                }));
        }

        for (auto& w : workers) { w.join(); }
        rparser.stop();
    }

    spdlog_piscem::info("finished mapping.");
    if (!ptab.empty()) {
      spdlog_piscem::info("number of reads discarded because of poison k-mers: {}", global_np);
    }
    // rewind to the start of the file and write the number of
    // chunks that we actually produced.
    out_info.rad_file.seekp(chunk_offset);
    uint64_t nc = out_info.num_chunks.load();
    out_info.rad_file.write(reinterpret_cast<char*>(&nc), sizeof(nc));

    out_info.rad_file.close();

    // We want to check if the RAD file stream was written to
    // properly. While we likely would have caught this earlier,
    // it is possible the badbit may not be set until the stream
    // actually flushes (perhaps even at close), so we check here
    // one final time that the status of the stream is as
    // expected. see :
    // https://stackoverflow.com/questions/28342660/error-handling-in-stdofstream-while-writing-data
    if (!out_info.rad_file) {
        spdlog_piscem::critical(
            "The RAD file stream had an invalid status after "
            "close; so some operation(s) may"
            "have failed!\nA common cause for this is lack "
            "of output disk space.\n"
            "Consequently, the output may be corrupted.\n\n");
        return 1;
    }

    auto end_t = std::chrono::high_resolution_clock::now();
    auto num_sec = std::chrono::duration_cast<std::chrono::seconds>(end_t - start_t);
    piscem::meta_info::run_stats rs;
    rs.cmd_line(cmdline);
    rs.num_reads(global_nr.load());
    rs.num_hits(global_nh.load());
    rs.num_poisoned(global_np.load());
    rs.num_seconds(num_sec.count());

    ghc::filesystem::path map_info_file_path = output_stem + ".map_info.json";
    bool info_ok = piscem::meta_info::write_map_info(rs, map_info_file_path);
    if (!info_ok) { spdlog_piscem::critical("failed to write map_info.json file"); }

    return 0;
}
