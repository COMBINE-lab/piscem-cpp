#include "../external/sshash/external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "../include/reference_index.hpp"
#include "../include/CanonicalKmerIterator.hpp"
#include "../include/Kmer.hpp"
#include "../include/streaming_query.hpp"
#include "../include/projected_hits.hpp"
#include "../include/util_piscem.hpp"
#include "../external/sshash/include/util.hpp"
#include "../include/mapping/utils.hpp"
#include "../include/mapping/utils_bin.hpp"
#include "../include/parallel_hashmap/phmap.h"
#include "../include/FastxParser.hpp"
#include "../include/rad/rad_writer.hpp"
#include "../include/rad/rad_header.hpp"
#include "../include/rad/util.hpp"
#include "../include/ghc/filesystem.hpp"
#include "../include/cli11/CLI11.hpp"
#include "../include/sc/util.hpp"
#include "../include/spdlog_piscem/spdlog.h"
#include "../include/spdlog_piscem/sinks/stdout_color_sinks.h"
#include "../include/meta_info.hpp"
#include "../include/mapping/utils.hpp"
// #include "check_overlap.cpp"
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
using mapping::util::poison_state_t;
using mapping::util::bin_pos;

struct pesc_atac_options {
    std::string index_basename;
    std::vector<std::string> left_read_filenames;
    std::vector<std::string> right_read_filenames;
    std::vector<std::string> barcode_filenames;
    std::string output_dirname;
    bool no_poison{true};
    bool enable_structural_constraints{false};
    float thr{0.7};
    bool quiet{false};
    bool check_ambig_hits{false};
    uint32_t max_ec_card{256};
    uint64_t bin_size{20000};
    uint64_t overlap{300};
    size_t nthread{16};
};

template <typename mapping_cache_info_t>
bool map_fragment(fastx_parser::ReadTrip& record, 
                  poison_state_t& poison_state,
                  mapping_cache_info_t& map_cache_left,
                  mapping_cache_info_t& map_cache_right,
                  mapping_cache_info_t& map_cache_out, 
                  std::atomic<uint64_t>& k_match, 
                  std::atomic<uint64_t> &l_match,
                  std::atomic<uint64_t> &r_match,
                  mapping::util::bin_pos& binning) {

    bool km = false; //kmatch checker
    map_cache_out.clear();
    poison_state.clear();
    poison_state.set_fragment_end(mapping::util::fragment_end::LEFT);
    // std::cout << "left\n";
    bool early_exit_left = mapping::util::map_read(&record.first.seq, map_cache_left, poison_state, binning, km);
    if (poison_state.is_poisoned()) {
        return false;
    }

    bool right_km = false;
    poison_state.set_fragment_end(mapping::util::fragment_end::RIGHT);
    // std::cout << "right\n";
    bool early_exit_right = mapping::util::map_read(&record.second.seq, map_cache_right, poison_state, binning, right_km);
    if (poison_state.is_poisoned()) {
        return false;
    }  
    if(km or right_km) {
        ++k_match;
    }

    int32_t left_len = static_cast<int32_t>(record.first.seq.length());
    int32_t right_len = static_cast<int32_t>(record.second.seq.length());

    l_match += map_cache_left.accepted_hits.empty() ? 0 : 1;
    r_match += map_cache_right.accepted_hits.empty() ? 0 : 1;

    mapping::util_bin::merge_se_mappings(map_cache_left, map_cache_right, left_len, right_len,
                                     map_cache_out);
    // std::cout << "merged size " << map_cache_out.accepted_hits.size() << std::endl;
    return (early_exit_left or early_exit_right);
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
    // BED records are written
    std::ofstream bed_file;
    // the mutex for safely writing to
    // bed_file
    std::mutex bed_mutex;

    // the output stream where the actual
    // RAD records are written
    std::ofstream rad_file;
    // the mutex for safely writing to
    // rad_file
    std::mutex rad_mutex;

    // the output stream where the counts
    // of observed barcodes for unmapped
    // reads will be written
    std::ofstream unmapped_bc_file;
    // the mutex for safely writing to
    // unmapped_bc_file
    std::mutex unmapped_bc_mutex;

    // the file where the sequence from first fastq and string mapping will be stored
};

template <typename FragT, typename SketchHitT>
void do_map(mindex::reference_index& ri,
                    fastx_parser::FastxParser<FragT>& parser,
                    mapping::util::bin_pos& binning, 
                    poison_table &poison_map,  
                    std::atomic<uint64_t>& global_nr,
                    std::atomic<uint64_t>& global_nhits,
                    std::atomic<uint64_t>& global_nmult,
                    std::atomic<uint64_t>& k_match,
                    std::atomic<uint64_t>& l_match,
                    std::atomic<uint64_t>& r_match,
                    std::atomic<uint64_t>& global_npoisoned,
                    pesc_output_info& out_info,
                    std::mutex& iomut,
                    RAD::RAD_Writer& rw,
                    RAD::Token token) {

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
    poison_state_t poison_state;
    if (use_poison) {
        poison_state.ptab = &poison_map;
    }
    // the reads are paired
    if constexpr (std::is_same_v<fastx_parser::ReadTrip, FragT>) {
        poison_state.paired_for_mapping = true;
    }

    mapping_cache_info<SketchHitT> map_cache_left(ri);
    mapping_cache_info<SketchHitT> map_cache_right(ri);
    mapping_cache_info<SketchHitT> map_cache_out(ri);

    size_t max_chunk_reads = 5000;

    auto rg = parser.getReadGroup();
    uint32_t num_reads_in_chunk{0};

    piscem::streaming_query q(ri.get_dict());
    // for(int32_t i = 0; i < ri.num_refs(); i++) {
    //     std::cout << "i is " <<  i << " ref len is " << ri.ref_len(i) << " name is " << ri.ref_name(i) << std::endl;
    // }

    mindex::hit_searcher hs(&ri);
    uint64_t read_num = 0;

    std::string temp_buff="";
    while (parser.refill(rg)) {
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
            std::string* bc = &record.third.seq; // need to modify this
            bc_kmer_t bc_kmer;
        
            auto recovered = single_cell::util::recover_barcode(*bc);
            if (recovered == BarCodeRecovered::NOT_RECOVERED) { continue; }
            
            bool bc_ok = bc_kmer.fromChars(*bc);
            if (!bc_ok) { continue; }
            
            // if we couldn't correct it with 1 `N`, then skip.
            
            bool had_early_stop = 
               map_fragment(record, poison_state, map_cache_left, map_cache_right, map_cache_out, k_match,
                   l_match, r_match, binning);
            (void)had_early_stop;
            if (poison_state.is_poisoned()) {
                global_npoisoned++;
            }    
            
            global_nhits += map_cache_out.accepted_hits.empty() ? 0 : 1;
            global_nmult += map_cache_out.accepted_hits.size() > 1 ? 1 : 0;
            rad::util::write_to_rad_stream_atac(bc_kmer, map_cache_out.map_type, map_cache_out.accepted_hits,
                                                map_cache_out.unmapped_bc_map, num_reads_in_chunk, 
                                                temp_buff, *bc, ri, rw, token);
            
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
    }
    
    if (num_reads_in_chunk > 0) {
        out_info.num_chunks++;
        
        out_info.bed_mutex.lock();
        out_info.bed_file << temp_buff;
        out_info.bed_mutex.unlock();
        temp_buff = "";
        num_reads_in_chunk = 0;
    }

    // // unmapped barcode writer
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
	
    // std::cerr << "left new state: " << map_cache_left.hs.new_state_cnt << "\n";
    // std::cerr << "left matches : " << map_cache_left.hs.matches_cnt << "\n";
    // std::cerr << "left non matches : " << map_cache_left.hs.non_matches_cnt << "\n";
}

#ifdef __cplusplus
extern "C" {
#endif
  int run_pesc_sc_atac(int argc, char** argv);
#ifdef __cplusplus
}
#endif

int run_pesc_sc_atac(int argc, char** argv) {
    std::ios_base::sync_with_stdio(false);
    pesc_atac_options po;
    std::string skipping_rule;

    size_t nthread{16};
    CLI::App app{"Single cell Atac Seq mapper"};
    app.add_option("-i,--index", po.index_basename, "input index prefix")->required();
    app.add_option("-1,--read1", po.left_read_filenames, "path to list of read 1 files")
        ->required()
        ->delimiter(',');
    app.add_option("-2,--read2", po.right_read_filenames, "path to list of read 2 files")
        ->delimiter(',');
    app.add_option("-b,--barcode", po.barcode_filenames, "path to list of barcodes")
        ->required()
        ->delimiter(',');
    app.add_option("-o,--output", po.output_dirname, "path to output directory")->required();
    app.add_option("-t,--threads", po.nthread,
                   "An integer that specifies the number of threads to use")
        ->default_val(16);
    app.add_option("--no-poison", po.no_poison,
                "Do not filter reads for poison k-mers, even if a poison table "
                "exists for the index")
        ->default_val(true);
    app.add_flag("-c,--struct-constraints", po.enable_structural_constraints,
                "Apply structural constraints when performing mapping");
    app
        .add_option(
        "--skipping-strategy", skipping_rule,
        "Which skipping rule to use for pseudoalignment ({strict, permissive, strict})")
        ->default_val("permissive");
    app.add_flag("--quiet", po.quiet,
                "Try to be quiet in terms of console output");
    app.add_option("--thr", po.thr,
                "threshold for psa")
        ->default_val(0.7);
    app.add_option("--bin_size", po.bin_size,
                "size for binning")
        ->default_val(20000);
    app.add_option("--overlap", po.overlap,
                "size for overlap")
        ->default_val(300);
    auto check_ambig =
        app.add_flag("--check-ambig-hits", po.check_ambig_hits,
                    "check the existence of highly-frequent hits in mapped targets, rather than "
                    "ignoring them.");
    (void) check_ambig; // currently unused in atacseq mode
    
    CLI11_PARSE(app, argc, argv);
    bool paired_end = !po.right_read_filenames.empty();
    
    spdlog_piscem::drop_all();
    auto logger = spdlog_piscem::create<spdlog_piscem::sinks::stdout_color_sink_mt>("");
    logger->set_pattern("%+");

    if (po.quiet) { spdlog_piscem::set_level(spdlog_piscem::level::warn); }
    spdlog_piscem::set_default_logger(logger);

    spdlog_piscem::info("enable structural constraints : {}",
                      po.enable_structural_constraints);

    std::optional<mindex::SkippingStrategy> skip_strat_opt =
        mindex::SkippingStrategy::from_string(skipping_rule);
    if (!skip_strat_opt) {
        spdlog_piscem::critical("The skipping strategy must be one of \"strict\", "
                                " \"every\" or \"permissive\", but \"{}\" was passed in",
                                skipping_rule);
        return 1;
    }

    nthread = po.nthread;

    auto start_t = std::chrono::high_resolution_clock::now();

    ghc::filesystem::path output_path(po.output_dirname);
    ghc::filesystem::create_directories(output_path);

    ghc::filesystem::path bed_file_path = output_path / "map.bed";
    ghc::filesystem::path unmapped_bc_file_path = output_path / "unmapped_bc_count.bin";
    ghc::filesystem::path mapped_bc_file_path = output_path / "mapped_bc.txt";

    std::ofstream bed_file(bed_file_path.string());

    if (!bed_file.good()) {
        spdlog_piscem::critical("Could not open {} for writing.", bed_file_path.string());
        throw std::runtime_error("error creating bed file.");
    }

    std::ofstream unmapped_bc_file(unmapped_bc_file_path.string());
    if (!unmapped_bc_file.good()) {
        spdlog_piscem::critical("Could not open {} for writing.", unmapped_bc_file_path.string());
        throw std::runtime_error("error creating unmapped barcode file.");
    }

    poison_table ptab;
    if (!po.no_poison and poison_table::exists(po.index_basename)) {
        /*
        spdlog_piscem::info("Loading poison k-mer map...");
        phmap::BinaryInputArchive ar_in(pmap_file.c_str());
        poison_map.phmap_load(ar_in);
        spdlog_piscem::info("done");
        */
        poison_table ptab_tmp(po.index_basename);
        ptab = std::move(ptab_tmp);
    } else {
        spdlog_piscem::info(
        "No poison k-mer map exists, or it was requested not to be used");
    }

    pesc_output_info out_info;
    out_info.bed_file = std::move(bed_file);
    out_info.unmapped_bc_file = std::move(unmapped_bc_file);
    
    mindex::reference_index ri(po.index_basename, po.check_ambig_hits);
    std::string rad_file_path = output_path;
    rad_file_path.append("/map.rad");
    RAD::Tag_Defn tag_defn;
    RAD::Tag_List file_tag_vals;
    file_tag_vals.add(RAD::Type::u16(16));

    std::vector<std::string> refs;
    bc_kmer_t::k(16);
    rad::util::write_rad_header_atac(ri, refs, tag_defn);
    mapping::util::bin_pos binning(&ri, po.thr, po.bin_size, po.overlap);
    
    const RAD::Header header(static_cast<uint8_t>(paired_end), refs.size(), refs);
    
    RAD::RAD_Writer rw(header, tag_defn, file_tag_vals, rad_file_path, nthread);
    
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
    std::atomic<uint64_t> global_nmult{0}; // number of multimapping
    std::atomic<uint64_t> k_match{0}; //whether the kmer exists in the unitig table
    std::atomic<uint64_t> l_match{0};
    std::atomic<uint64_t> r_match{0};
    std::atomic<uint64_t> global_np{0}; //whether the kmer exists in the unitig table
    std::mutex iomut;

    if (paired_end) {
        using FragmentT = fastx_parser::ReadTrip;
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
                [&ri, &po, &rparser, &binning, &ptab, &global_nr, &global_nh, &global_nmult, &k_match, &global_np,
                 &out_info, &iomut, &rw, &l_match, &r_match]() {
                    const auto token = rw.get_token();
                    if (!po.enable_structural_constraints) {
                        using SketchHitT = mapping::util::sketch_hit_info_no_struct_constraint;
                        do_map<FragmentT, SketchHitT>(ri, rparser, binning, ptab, global_nr, global_nh, global_nmult,
                            k_match, l_match, r_match, global_np, out_info, iomut, rw, token);
                    } else {
                        using SketchHitT = mapping::util::sketch_hit_info;   
                        do_map<FragmentT, SketchHitT>(ri, rparser, binning, ptab, global_nr, global_nh, global_nmult,
                            k_match, l_match, r_match, global_np, out_info, iomut, rw, token);
                    }
                }));
        }
    
        for (auto& w : workers) { w.join(); }
        rparser.stop();
    }

    spdlog_piscem::info("finished mapping.");
    rw.close();
    out_info.bed_file.close();
    
    if (!out_info.bed_file) {
        spdlog_piscem::critical(
            "The BED file stream had an invalid status after "
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
    rs.num_lmatch(l_match.load());
    rs.num_rmatch(r_match.load());

    ghc::filesystem::path map_info_file_path = output_path / "map_info.json";
    bool info_ok = piscem::meta_info::write_map_info(rs, map_info_file_path);
    if (!info_ok) { spdlog_piscem::critical("failed to write map_info.json file"); }

    return 0;
}
