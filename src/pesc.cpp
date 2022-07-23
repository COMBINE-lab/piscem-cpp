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
#include "../include/json.hpp"
#include "../include/mapping/utils.hpp"
#include "FastxParser.cpp"
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
using umi_kmer_t = rad::util::umi_kmer_t;
using bc_kmer_t = rad::util::bc_kmer_t;

// utility class that wraps the information we will
// need access to when writing output within each thread
// as well as information we'll need to update for the
// caller.
class pesc_output_info {
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

    // the output stream where the counts
    // of observed barcodes for unmapped
    // reads will be written
    std::ofstream unmapped_bc_file;
    // the mutex for safely writing to
    // unmapped_bc_file
    std::mutex unmapped_bc_mutex;
};

template <typename Protocol>
void do_map(mindex::reference_index& ri, fastx_parser::FastxParser<fastx_parser::ReadPair>& parser,
            const Protocol& p, std::atomic<uint64_t>& global_nr,
            std::atomic<uint64_t>& global_nhits, pesc_output_info& out_info, std::mutex& iomut) {
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

    // put these in struct
    size_t num_short_umi{0};
    size_t num_ambig_umi{0};

    mapping::util::mapping_cache_info map_cache(ri);

    size_t max_chunk_reads = 5000;
    // Get the read group by which this thread will
    // communicate with the parser (*once per-thread*)
    auto rg = parser.getReadGroup();

    Protocol protocol(p);
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
            auto rctr = global_nr.load();
            auto hctr = global_nhits.load();

            if (write_mapping_rate and (rctr % 500000 == 0)) {
                iomut.lock();
                std::cerr << "\rprocessed (" << rctr << ") reads; (" << hctr << ") had mappings.";
                iomut.unlock();
            }

            // first extract the barcode
            std::string* bc = protocol.extract_bc(record.first.seq, record.second.seq);
            // if we couldn't get it, don't bother with
            // anything else for this read.
            if (bc == nullptr) { continue; }

            // correct up to one `N` in the barcode
            auto recovered = single_cell::util::recover_barcode(*bc);
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

            bool had_early_stop = mapping::util::map_read(read_seq, map_cache);
            (void)had_early_stop;

            global_nhits += map_cache.accepted_hits.empty() ? 0 : 1;
            rad::util::write_to_rad_stream(bc_kmer, umi_kmer, map_cache.map_type,
                                           map_cache.accepted_hits, map_cache.unmapped_bc_map,
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
        for (auto& kv : map_cache.unmapped_bc_map) {
            ubcw << kv.first;
            ubcw << kv.second;
        }
        out_info.unmapped_bc_mutex.lock();
        out_info.unmapped_bc_file << ubcw;
        out_info.unmapped_bc_mutex.unlock();
        ubcw.clear();
    }
}

enum class protocol_t : uint8_t { CHROM_V2, CHROM_V3, CUSTOM };

bool set_geometry(std::string& library_geometry, protocol_t& pt,
                  std::unique_ptr<custom_protocol>& p) {
    if (library_geometry == "chromium_v2") {
        umi_kmer_t::k(10);
        bc_kmer_t::k(16);
        pt = protocol_t::CHROM_V2;
    } else if (library_geometry == "chromium_v3") {
        umi_kmer_t::k(12);
        bc_kmer_t::k(16);
        pt = protocol_t::CHROM_V3;
    } else {
        std::unique_ptr<custom_protocol> opt_cp =
            single_cell::util::parse_custom_geometry(library_geometry);
        if (opt_cp) {
            p.swap(opt_cp);
            umi_kmer_t::k(p->get_umi_len());
            bc_kmer_t::k(p->get_bc_len());
            pt = protocol_t::CUSTOM;
        } else {
            spdlog::critical("could not parse custom geometry description [{}]", library_geometry);
            return false;
        }
    }
    return true;
}

class run_stats {
public:
    run_stats() = default;
    void cmd_line(std::string& cmd_line_in) { cmd_line_ = cmd_line_in; }
    void num_reads(uint64_t num_reads_in) { num_reads_ = num_reads_in; }
    void num_hits(uint64_t num_hits_in) { num_hits_ = num_hits_in; }
    void num_seconds(double num_sec) { num_seconds_ = num_sec; }

    std::string cmd_line() const { return cmd_line_; }
    uint64_t num_reads() const { return num_reads_; }
    uint64_t num_hits() const { return num_hits_; }
    double num_seconds() const { return num_seconds_; }

private:
    std::string cmd_line_{""};
    uint64_t num_reads_{0};
    uint64_t num_hits_{0};
    double num_seconds_{0};
};

bool write_map_info(run_stats& rs, ghc::filesystem::path& map_info_file_path) {
    using json = nlohmann::json;

    json j;
    j["cmdline"] = rs.cmd_line();
    j["num_reads"] = rs.num_reads();
    j["num_mapped"] = rs.num_hits();
    double percent_mapped = (100.0 * static_cast<double>(rs.num_hits())) / rs.num_reads();
    j["percent_mapped"] = percent_mapped;
    j["runtime_seconds"] = rs.num_seconds();
    // write prettified JSON to another file
    std::ofstream o(map_info_file_path.string());

    if (!o.good()) { return false; }

    o << std::setw(4) << j << std::endl;

    if (!o) { return false; }
    return true;
}

struct pesc_options {
  std::string index_basename;
  std::vector<std::string> left_read_filenames;
  std::vector<std::string> right_read_filenames;
  std::string output_dirname;
  std::string library_geometry;
  protocol_t pt{protocol_t::CUSTOM};
  std::unique_ptr<custom_protocol> p{nullptr};
  bool quiet{false};
  size_t nthread{16};
};

#ifdef __cplusplus
extern "C" {
#endif
  int run_pesc(int argc, char** argv);
#ifdef __cplusplus
}
#endif

int run_pesc(int argc, char** argv) {
    /**
     * PESC : Pseudoalignment Enhanced with Structural Constraints
     **/
    std::ios_base::sync_with_stdio(false);

    pesc_options po;
    CLI::App app{"PESC — single-cell RNA-seq mapper for alevin-fry"};
    app.add_option("-i,--index", po.index_basename, "input index prefix")->required();
    app.add_option("-1,--read1", po.left_read_filenames, "path to list of read 1 files")
        ->required()
        ->delimiter(',');
    app.add_option("-2,--read2", po.right_read_filenames, "path to list of read 2 files")
        ->required()
        ->delimiter(',');
    app.add_option("-o,--output", po.output_dirname, "path to output directory")->required();
    app.add_option("-g,--geometry", po.library_geometry, "geometry of barcode, umi and read")
        ->required();
    app.add_option("-t,--threads", po.nthread,
                   "An integer that specifies the number of threads to use")
        ->default_val(16);
    app.add_flag("--quiet", po.quiet, "try to be quiet in terms of console output");
    CLI11_PARSE(app, argc, argv);

    spdlog::drop_all();
    auto logger = spdlog::create<spdlog::sinks::stdout_color_sink_mt>("");
    logger->set_pattern("%+");

    if (po.quiet) {
      logger->set_level(spdlog::level::warn);
    }
    spdlog::set_default_logger(logger);

    // start the timer
    auto start_t = std::chrono::high_resolution_clock::now();

    bool geom_ok = set_geometry(po.library_geometry, po.pt, po.p);
    if (!geom_ok) {
        spdlog::critical("could not set the library geometry properly.");
        return 1;
    }

    // RAD file path
    ghc::filesystem::path output_path(po.output_dirname);
    ghc::filesystem::create_directories(output_path);

    ghc::filesystem::path rad_file_path = output_path / "map.rad";
    ghc::filesystem::path unmapped_bc_file_path = output_path / "unmapped_bc_count.bin";

    mindex::reference_index ri(po.index_basename);

    std::string cmdline;
    size_t narg = static_cast<size_t>(argc);
    for (size_t i = 0; i < narg; ++i) {
        cmdline += std::string(argv[i]);
        cmdline.push_back(' ');
    }
    cmdline.pop_back();

    std::ofstream rad_file(rad_file_path.string());
    std::ofstream unmapped_bc_file(unmapped_bc_file_path.string());

    if (!rad_file.good()) {
        spdlog::critical("Could not open {} for writing.", rad_file_path.string());
        throw std::runtime_error("error creating output file.");
    }

    if (!unmapped_bc_file.good()) {
        spdlog::critical("Could not open {} for writing.", unmapped_bc_file_path.string());
        throw std::runtime_error("error creating output file.");
    }

    pesc_output_info out_info;
    out_info.rad_file = std::move(rad_file);
    out_info.unmapped_bc_file = std::move(unmapped_bc_file);

    size_t bc_length = bc_kmer_t::k();
    size_t umi_length = umi_kmer_t::k();
    size_t chunk_offset = rad::util::write_rad_header(ri, bc_length, umi_length, out_info.rad_file);

    std::atomic<uint64_t> num_chunks{0};
    std::mutex iomut;

    uint32_t np = 1;
    if ((po.left_read_filenames.size() > 1) and (po.nthread >= 6)) {
        np += 1;
        po.nthread -= 1;
    }

    fastx_parser::FastxParser<fastx_parser::ReadPair> rparser(po.left_read_filenames, po.right_read_filenames, po.nthread, np);
    rparser.start();

    // set the k-mer size for the
    // CanonicalKmer type.
    CanonicalKmer::k(ri.k());

    std::atomic<uint64_t> global_nr{0};
    std::atomic<uint64_t> global_nh{0};
    std::vector<std::thread> workers;
    for (size_t i = 0; i < po.nthread; ++i) {
        switch (po.pt) {
            case protocol_t::CHROM_V2: {
                chromium_v2 prot;
                workers.push_back(std::thread(
                    [&ri, &rparser, &prot, &global_nr, &global_nh, &iomut, &out_info]() {
                        do_map(ri, rparser, prot, global_nr, global_nh, out_info, iomut);
                    }));
                break;
            }
            case protocol_t::CHROM_V3: {
                chromium_v3 prot;
                workers.push_back(std::thread(
                    [&ri, &rparser, &prot, &global_nr, &global_nh, &iomut, &out_info]() {
                        do_map(ri, rparser, prot, global_nr, global_nh, out_info, iomut);
                    }));
                break;
            }
            case protocol_t::CUSTOM: {
                workers.push_back(
                    std::thread([&ri, &rparser, &po, &global_nr, &global_nh, &iomut, &out_info]() {
                        do_map(ri, rparser, *(po.p), global_nr, global_nh, out_info, iomut);
                    }));
                break;
            }
        }
    }

    for (auto& w : workers) { w.join(); }
    rparser.stop();

    spdlog::info("finished mapping.");
    
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
    run_stats rs;
    rs.cmd_line(cmdline);
    rs.num_reads(global_nr.load());
    rs.num_hits(global_nh.load());
    rs.num_seconds(num_sec.count());

    ghc::filesystem::path map_info_file_path = output_path / "map_info.json";
    bool info_ok = write_map_info(rs, map_info_file_path);
    if (!info_ok) { spdlog::critical("failed to write map_info.json file"); }

    return 0;
}
