#include "../include/reference_index.hpp"
#include "../include/util.hpp"
#include "../include/mapping/utils.hpp"
#include "../include/spdlog/spdlog.h"
#include "../include/spdlog/sinks/stdout_color_sinks.h"
#include "../include/cli11/CLI11.hpp"
#include "../include/FastxParser.hpp"
#include "FastxParser.cpp"
#include "zlib.h"
//#include "../src/hit_searcher.cpp"

#include <atomic>
#include <iostream>
#include <vector>
#include <numeric>
#include <cstdio>
#include <thread>
#include <sstream>

using namespace klibpp;
using mapping::util::mapping_cache_info;

void print_header(mindex::reference_index& ri, std::string& cmdline) {
    std::cout << "@HD\tVN:1.0\tSO:unsorted\n";
    for (uint64_t i = 0; i < ri.num_refs(); ++i) {
        std::cout << "@SQ\tNS:" << ri.ref_name(i) << "\tLN:" << ri.ref_len(i) << "\n";
    }
    std::cout << "@PG\tID:mindex_map\tPN:mapper\tVN:0.0.1\t"
              << "CL:" << cmdline << "\n";
}

// single-end
bool map_fragment(fastx_parser::ReadSeq& record, mapping_cache_info& map_cache_left,
                  mapping_cache_info& map_cache_right, mapping_cache_info& map_cache_out) {
    (void)map_cache_left;
    (void)map_cache_right;
    return mapping::util::map_read(&record.seq, map_cache_out);
}

// paried-end
bool map_fragment(fastx_parser::ReadPair& record, mapping_cache_info& map_cache_left,
                  mapping_cache_info& map_cache_right, mapping_cache_info& map_cache_out) {
    bool early_exit_left = mapping::util::map_read(&record.first.seq, map_cache_left);
    bool early_exit_right = mapping::util::map_read(&record.second.seq, map_cache_right);

    int32_t left_len = static_cast<int32_t>(record.first.seq.length());
    int32_t right_len = static_cast<int32_t>(record.second.seq.length());
    mapping::util::merge_se_mappings(map_cache_left, map_cache_right, left_len, right_len,
                                     map_cache_out);

    return (early_exit_left or early_exit_right);
}

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
              pos_first = ah.pos;
              pos_second = ah.mate_pos;

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

              if (!ah.is_fw) { // if the mapped read is rc
                flag_first += is_rc;
                if (!have_rc_first){
                  have_rc_first = true;
                  combinelib::kmers::reverseComplement(record.first.seq, workstr_left);
                }
                sptr_first = &workstr_left;

                flag_second += mate_rc;
              } 

            } else if (map_type == mapping::util::MappingType::MAPPED_SECOND_ORPHAN) {

              pos_first = 0;
              pos_second = ah.pos; 

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

            const auto ref_name = map_cache_out.hs.get_index()->ref_name(ah.tid);
            osstream << record.first.name << "\t" << flag_first << "\t"
                     << ((flag_first & unmapped) ? "*" : ref_name)  << '\t' // if mapped RNAME, else *
                     << pos_first << '\t' // POS 
                     << "255\t*\t" // MAPQ & CIGAR
                     << ((flag_first & mate_unmapped) ? '*' : '=') << '\t' // RNEXT
                     << pos_second << '\t' // PNEXT
                     << ah.frag_len() << '\t' 
                     << *sptr_first << "\t*\n";
            osstream << record.second.name << "\t" << flag_second << "\t"
                     << ((flag_second & unmapped) ? "*" : ref_name)  << '\t' // if mapped RNAME, else *
                     << pos_second << '\t' // POS 
                     << "255\t*\t" // MAPQ & CIGAR
                     << ((flag_second & mate_unmapped) ? '*' : '=') << '\t' // RNEXT
                     << pos_first << '\t' // PNEXT
                     << -ah.frag_len() << '\t' 
                     << *sptr_second << "\t*\n";
            secondary = true;
        }
    } else {
        osstream << record.first.name << "\t" << 77 << "\t"
                 << "*\t0\t0\t*\t*\t0\t0\t" << record.first.seq << "\t*\n";
        osstream << record.second.name << "\t" << 141 << "\t"
                 << "*\t0\t0\t*\t*\t0\t0\t" << record.second.seq << "\t*\n";
    }
}

template <typename FragT>
void do_map(mindex::reference_index& ri, fastx_parser::FastxParser<FragT>& parser,
            std::atomic<uint64_t>& global_nr, std::atomic<uint64_t>& global_nhits,
            std::mutex& iomut) {
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

    // set the canonical k-mer size globally
    CanonicalKmer::k(ri.k());

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
