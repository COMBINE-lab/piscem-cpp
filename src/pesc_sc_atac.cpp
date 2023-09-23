#include "../external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "../include/reference_index.hpp"
#include "../include/CanonicalKmerIterator.hpp"
#include "../include/Kmer.hpp"
#include "../include/query/streaming_query_canonical_parsing.hpp"
#include "../include/projected_hits.hpp"
#include "../include/util.hpp"
#include "../include/mapping/utils.hpp"
#include "../include/parallel_hashmap/phmap.h"
#include "../include/FastxParser.hpp"
#include "../include/rad/rad_writer.hpp"
#include "../include/rad/rad_header.hpp"
#include "../include/rad/util.hpp"
#include "../include/ghc/filesystem.hpp"
#include "../include/cli11/CLI11.hpp"
#include "../include/sc/util.hpp"
#include "../include/spdlog/spdlog.h"
#include "../include/spdlog/sinks/stdout_color_sinks.h"
#include "../include/meta_info.hpp"
#include "../include/mapping/utils.hpp"
#include "check_overlap.cpp"
//#include "FastxParser.cpp"
//#include "hit_searcher.cpp"
#include "zlib.h"

#include <atomic>
#include <chrono>
#include <iostream>
#include <vector>
#include <memory>
#include <numeric>
#include <optional>
#include <cstdio>
#include <thread>
#include <sstream>

using namespace klibpp;
using BarCodeRecovered = single_cell::util::BarCodeRecovered;
using bc_kmer_t = rad::util::bc_kmer_t;
using mapping::util::mapping_cache_info;

struct pesc_atac_options {
    std::string index_basename;
    std::vector<std::string> left_read_filenames;
    std::vector<std::string> right_read_filenames;
    std::vector<std::string> barcode_filenames;
    std::string output_dirname;
    bool psc_off{false};
    bool ps_skip{true};
    bool quiet{false};
    bool check_ambig_hits{false};
    uint32_t max_ec_card{256};
    size_t nthread{16};
};

// bool map_fragment(fastx_parser::ReadSeq& record, mapping_cache_info& map_cache_left,
//                   mapping_cache_info& map_cache_right, mapping_cache_info& map_cache_out) {
//     (void)map_cache_left;
//     (void)map_cache_right;
//     return mapping::util::map_read(&record.seq, map_cache_out, false);
// }

// paried-end
bool map_fragment(fastx_parser::ReadTrip& record, mapping_cache_info& map_cache_left,
                  mapping_cache_info& map_cache_right, mapping_cache_info& map_cache_out, 
                  std::atomic<uint64_t>& k_match, bool psc_off, bool ps_skip) {
// bool map_fragment(fastx_parser::ReadTrip& record, mapping_cache_info& map_cache_left,
//                   mapping_cache_info& map_cache_right, mapping_cache_info& map_cache_out) {
    check_overlap::MateOverlap mov;
    check_overlap::findOverlapBetweenPairedEndReads(record.first.seq, record.second.seq, mov, 30);
    bool km = false; //kmatch checker
    // std::cout << mov.frag << std::endl;
    // std::cout << "aaa\n" ;
    if(mov.frag != "") {
        (void)map_cache_left;
        (void)map_cache_right;
        bool read_map = mapping::util::map_atac_read(&mov.frag, map_cache_out, false, km, psc_off, ps_skip);
        if (km) {
            ++k_match;
        }
        return read_map;
    }
    return true;
    // bool early_exit_left = mapping::util::map_atac_read(&record.first.seq, map_cache_left, false, km, psc_off, ps_skip);
    
    // bool right_km = false;
    // bool early_exit_right = mapping::util::map_atac_read(&record.second.seq, map_cache_right, false, right_km, psc_off, ps_skip);

    // if(km | right_km) {
    //     ++k_match;
    // }

    // int32_t left_len = static_cast<int32_t>(record.first.seq.length());
    // int32_t right_len = static_cast<int32_t>(record.second.seq.length());
    
    // mapping::util::merge_se_mappings(map_cache_left, map_cache_right, left_len, right_len,
    //                                  map_cache_out);

    // return (early_exit_left or early_exit_right);
}

// utility class that wraps the information we will
// need access to when writing output within each thread
// as well as information we'll need to update for the
// caller.
class pesc_output_info {
public:
    // will keep track of the total number
    // of chunks written to be_file
    std::atomic<size_t> num_chunks{0};
    // the output stream where the actual
    // RAD records are written
    std::ofstream bed_file;
    // the mutex for safely writing to
    // rad_file
    std::mutex bed_mutex;

    // the output stream where the counts
    // of observed barcodes for unmapped
    // reads will be written
    std::ofstream unmapped_bc_file;
    // the mutex for safely writing to
    // unmapped_bc_file
    std::mutex unmapped_bc_mutex;
};

template <typename FragT>
void do_map(mindex::reference_index& ri,
                    fastx_parser::FastxParser<FragT>& parser, 
                    std::atomic<uint64_t>& global_nr, 
                    std::atomic<uint64_t>& global_nhits,
                    pesc_output_info& out_info,
                    std::mutex& iomut,
                    std::atomic<uint64_t>& k_match,
                    bool psc_off,
                    bool ps_skip) {

    auto log_level = spdlog::get_level();
    auto write_mapping_rate = false;
    switch (log_level) {
        case spdlog::level::level_enum::trace:
            write_mapping_rate = true;
            break;
        case spdlog::level::level_enum::debug:
            write_mapping_rate = true;
            break;
        case spdlog::level::level_enum::info:
            write_mapping_rate = true;
            break;
        case spdlog::level::level_enum::warn:
            write_mapping_rate = false;
            break;
        case spdlog::level::level_enum::err:
            write_mapping_rate = false;
            break;
        case spdlog::level::level_enum::critical:
            write_mapping_rate = false;
            break;
        case spdlog::level::level_enum::off:
            write_mapping_rate = false;
            break;
        default:
            write_mapping_rate = false;
    }

    mapping_cache_info map_cache_left(ri);
    mapping_cache_info map_cache_right(ri);
    mapping_cache_info map_cache_out(ri);

    size_t max_chunk_reads = 5000;

    auto rg = parser.getReadGroup();
    uint32_t num_reads_in_chunk{0};

    sshash::streaming_query_canonical_parsing q(ri.get_dict());
    mindex::hit_searcher hs(&ri);
    uint64_t read_num = 0;

    std::string temp_buff="";
    temp_buff.reserve(max_chunk_reads*(9*2+16+5+12+5));
    while (parser.refill(rg)) {
        for (auto& record : rg) {
            ++global_nr;
            ++read_num;
            auto rctr = global_nr.load();
            auto hctr = global_nhits.load();
            // if(read_num > 100) {
            //     break;
            // }

            if (write_mapping_rate and (rctr % 500000 == 0)) {
                iomut.lock();
                std::cerr << "\rprocessed (" << rctr << ") reads; (" << hctr << ") had mappings.";
                iomut.unlock();
            }
            std::string* bc = &record.third.seq; // need to modify this
            bc_kmer_t bc_kmer;

            auto recovered = single_cell::util::recover_barcode(*bc);
            // if we couldn't correct it with 1 `N`, then skip.
            if (recovered == BarCodeRecovered::NOT_RECOVERED) { continue; }
            bool had_early_stop = 
               map_fragment(record, map_cache_left, map_cache_right, map_cache_out, k_match,
                   psc_off, ps_skip);
            (void)had_early_stop;
            
            global_nhits += map_cache_out.accepted_hits.empty() ? 0 : 1;
            rad::util::write_to_rad_stream_atac(bc_kmer, map_cache_out.map_type, map_cache_out.accepted_hits,
                                                map_cache_out.unmapped_bc_map, num_reads_in_chunk, 
                                                temp_buff, *bc, ri);

            // dump buffer
            if (num_reads_in_chunk > max_chunk_reads) {
                out_info.num_chunks++;
                out_info.bed_mutex.lock();
                out_info.bed_file << temp_buff;
                out_info.bed_mutex.unlock();
                temp_buff = "";
                num_reads_in_chunk = 0;
            }
        }
        // if(read_num > 100) {
        //         break;
        // }
    }
    
    if (num_reads_in_chunk > 0) {
        out_info.num_chunks++;
        
        out_info.bed_mutex.lock();
        out_info.bed_file << temp_buff;
        out_info.bed_mutex.unlock();
        temp_buff = "";
        num_reads_in_chunk = 0;
    }

    // unmapped barcode writer
    {  // make a scope and dump the unmapped barcode counts
        std::string ubcw="";
        for (auto& kv : map_cache_out.unmapped_bc_map) {
            ubcw += kv.first;
            ubcw += kv.second;
        }
        out_info.unmapped_bc_mutex.lock();
        out_info.unmapped_bc_file << ubcw;
        out_info.unmapped_bc_mutex.unlock();
        ubcw.clear();
    }
}

#ifdef __cplusplus
extern "C" {
#endif
  int run_pesc_sc_atac(int argc, char** argv);
#ifdef __cplusplus
}
#endif

int run_pesc_sc_atac(int argc, char** argv) {
    // std::vector<std::string> read_file1 = {"/fs/cbcb-lab/rob/students/noor/Atacseq/l1.fastq"};
    // std::vector<std::string> read_file2 = {"/fs/cbcb-lab/rob/students/noor/Atacseq/l2.fastq"};
    // std::vector<std::string> barcode_file = {"/fs/cbcb-lab/rob/students/noor/Atacseq/b.fastq"};

    std::ios_base::sync_with_stdio(false);
    pesc_atac_options po;
    // // std::string index_basename;
    // // std::vector<std::string> left_read_filenames;
    // // std::vector<std::string> right_read_filenames;
    // // std::vector<std::string> single_read_filenames;
    // std::string output_stem;
    size_t nthread{16};
    bool quiet{false};
    CLI::App app{"Single cell Atac Seq mapper"};
    app.add_option("-i,--index", po.index_basename, "input index prefix")->required();
    app.add_option("-1,--read1", po.left_read_filenames, "path to list of read 1 files")
        ->required()
        ->delimiter(',');
    app.add_option("-2,--read2", po.right_read_filenames, "path to list of read 2 files")
        ->required()
        ->delimiter(',');
    app.add_option("-b,--barcode", po.barcode_filenames, "path to list of barcodes")
        ->required()
        ->delimiter(',');
    app.add_option("-o,--output", po.output_dirname, "path to output directory")->required();
    // app.add_option("-g,--geometry", po.library_geometry, "geometry of barcode, umi and read")
    //     ->required();
    app.add_option("-t,--threads", po.nthread,
                   "An integer that specifies the number of threads to use")
        ->default_val(16);
    app.add_option("--psc_off", po.psc_off,
                   "whether to switch structural constraints off")
        ->default_val(false);
    app.add_option("--ps_skip", po.ps_skip,
                   "whether to implement pseudoalignment with skipping")
        ->default_val(true);
    app.add_flag("--quiet", po.quiet, "try to be quiet in terms of console output");
    auto check_ambig =
        app.add_flag("--check-ambig-hits", po.check_ambig_hits,
                     "check the existence of highly-frequent hits in mapped targets, rather than "
                     "ignoring them.");
    
    CLI11_PARSE(app, argc, argv);
    
    spdlog::drop_all();
    auto logger = spdlog::create<spdlog::sinks::stdout_color_sink_mt>("");
    logger->set_pattern("%+");

    if (po.quiet) { logger->set_level(spdlog::level::warn); }
    spdlog::set_default_logger(logger);

    auto start_t = std::chrono::high_resolution_clock::now();

    ghc::filesystem::path output_path(po.output_dirname);
    ghc::filesystem::create_directories(output_path);

    ghc::filesystem::path bed_file_path = output_path / "map.bed";
    ghc::filesystem::path unmapped_bc_file_path = output_path / "unmapped_bc_count.bin";

    std::ofstream bed_file(bed_file_path.string());
    if (!bed_file.good()) {
        spdlog::critical("Could not open {} for writing.", bed_file_path.string());
        throw std::runtime_error("error creating output file.");
    }

    std::ofstream unmapped_bc_file(unmapped_bc_file_path.string());
    if (!unmapped_bc_file.good()) {
        spdlog::critical("Could not open {} for writing.", unmapped_bc_file_path.string());
        throw std::runtime_error("error creating output file.");
    }

    pesc_output_info out_info;
    out_info.bed_file = std::move(bed_file);
    out_info.unmapped_bc_file = std::move(unmapped_bc_file);

    mindex::reference_index ri(po.index_basename);
    std::string cmdline;
    size_t narg = static_cast<size_t>(argc);
    for (size_t i = 0; i < narg; ++i) {
        cmdline += std::string(argv[i]);
        cmdline.push_back(' ');
    }
    cmdline.pop_back();
    CanonicalKmer::k(ri.k());

    uint32_t np = 1;
    
    std::atomic<uint64_t> global_nr{0};
    std::atomic<uint64_t> global_nh{0};
    std::atomic<uint64_t> k_match{0}; //whether the kmer exists in the unitig table
    std::mutex iomut;

    bool psc_off=po.psc_off;
    bool ps_skip=po.ps_skip;

    fastx_parser::FastxParser<fastx_parser::ReadTrip> rparser(
    po.left_read_filenames, po.right_read_filenames, po.barcode_filenames, nthread, np);
    rparser.start();
    if (nthread >= 6) {
            np += 1;
            nthread -= 1;
    }
    std::vector<std::thread> workers;
    for (size_t i = 0; i < nthread; ++i) {
        workers.push_back(std::thread(
            [&ri, &rparser, &global_nr, &global_nh, &out_info, &iomut, &k_match, &psc_off, &ps_skip]() {
                do_map(ri, rparser, global_nr, global_nh, out_info, iomut, k_match, psc_off, ps_skip);
            }));
    }
    for (auto& w : workers) { w.join(); }
    rparser.stop();
    spdlog::info("finished mapping.");

    out_info.bed_file.close();
    if (!out_info.bed_file) {
        spdlog::critical(
            "The RAD file stream had an invalid status after "
            "close; so some operation(s) may"
            "have failed!\nA common cause for this is lack "
            "of output disk space.\n"
            "Consequently, the output may be corrupted.\n\n");
        return 1;
    }
    out_info.unmapped_bc_file.close();
        auto end_t = std::chrono::high_resolution_clock::now();
    auto num_sec = std::chrono::duration_cast<std::chrono::seconds>(end_t - start_t);
    piscem::meta_info::run_stats rs;
    rs.cmd_line(cmdline);
    rs.num_reads(global_nr.load());
    rs.num_hits(global_nh.load());
    rs.num_seconds(num_sec.count());
    rs.num_kmatch(k_match.load());

    ghc::filesystem::path map_info_file_path = output_path / "map_info.json";
    bool info_ok = piscem::meta_info::write_map_info(rs, map_info_file_path);
    if (!info_ok) { spdlog::critical("failed to write map_info.json file"); }

    return 0;
}