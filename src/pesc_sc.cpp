#include "../external/sshash/include/util.hpp"
#include "../include/CanonicalKmerIterator.hpp"
#include "../include/FastxParser.hpp"
#include "../include/Kmer.hpp"
#include "../include/cli11/CLI11.hpp"
#include "../include/defaults.hpp"
#include "../include/ghc/filesystem.hpp"
#include "../include/mapping/utils.hpp"
#include "../include/meta_info.hpp"
#include "../include/parallel_hashmap/phmap.h"
#include "../include/parallel_hashmap/phmap_dump.h"
#include "../include/poison_table.hpp"
#include "../include/projected_hits.hpp"
#include "../include/rad/rad_header.hpp"
#include "../include/rad/rad_writer.hpp"
#include "../include/rad/util.hpp"
#include "../include/reference_index.hpp"
#include "../include/sc/util.hpp"
#include "../include/spdlog_piscem/sinks/stdout_color_sinks.h"
#include "../include/spdlog_piscem/spdlog.h"
#include "../include/util_piscem.hpp"
#include "zlib.h"

#include <atomic>
#include <chrono>
#include <cstdio>
#include <iostream>
#include <memory>
#include <numeric>
#include <optional>
#include <sstream>
#include <thread>
#include <type_traits>
#include <vector>

using namespace klibpp;
using BarCodeRecovered = single_cell::util::BarCodeRecovered;
using umi_kmer_t = rad::util::umi_kmer_t;
using bc_kmer_t = rad::util::bc_kmer_t;
using mapping::util::mapping_cache_info;
using mapping::util::poison_state_t;

struct SingleEndBioSeq;
struct PairedEndBioSeq;

enum class protocol_t : uint8_t {
  CHROM_V2,
  CHROM_V2_5P,
  CHROM_V3,
  CHROM_V3_5P,
  CHROM_V4_3P,
  CUSTOM
};

struct pesc_sc_options {
  std::string index_basename;
  std::vector<std::string> left_read_filenames;
  std::vector<std::string> right_read_filenames;
  std::string output_dirname;
  std::string library_geometry;
  protocol_t pt{protocol_t::CUSTOM};
  std::unique_ptr<custom_protocol> p{nullptr};
  bool no_poison{false};
  bool quiet{false};
  bool enable_structural_constraints{false};
  bool ignore_ambig_hits{false};
  uint32_t max_ec_card{piscem::defaults::max_ec_card};
  uint32_t max_hit_occ{piscem::defaults::max_hit_occ};
  uint32_t max_hit_occ_recover{piscem::defaults::max_hit_occ_recover};
  uint32_t max_read_occ{piscem::defaults::max_read_occ};
  bool attempt_occ_recover{true};
  size_t nthread{16};
  mindex::SkippingStrategy skip_strat{mindex::SkippingStrategy::PERMISSIVE};
};

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

// single-end
template <typename mapping_cache_info_t>
bool map_se_fragment(AlignableReadSeqs &records, poison_state_t &poison_state,
                     mindex::SkippingStrategy skip_strat,
                     mapping_cache_info_t &map_cache_left,
                     mapping_cache_info_t &map_cache_right,
                     mapping_cache_info_t &map_cache_out) {
  (void)map_cache_left;
  (void)map_cache_right;
  poison_state.clear();
  return mapping::util::map_read(records.get_alignable_seq(), map_cache_out,
                                 poison_state, skip_strat);

  /*
  std::cerr << "\n\nnew_read: " << *records.get_alignable_seq() << "\n";
  for (auto& lh : map_cache_out.accepted_hits) {
    std::cerr << "single: " << lh.tid << ", " << lh.pos << " (" << (lh.is_fw ?
  "fw" : "rc") <<
  ")\n";
  }
  return r;
  */
}

// paried-end
template <typename mapping_cache_info_t>
bool map_pe_fragment(AlignableReadSeqs &records, poison_state_t &poison_state,
                     mindex::SkippingStrategy skip_strat,
                     mapping_cache_info_t &map_cache_left,
                     mapping_cache_info_t &map_cache_right,
                     mapping_cache_info_t &map_cache_out) {

  poison_state.clear();

  // we have to clear map_cache_out here in case we
  // exit early.
  map_cache_out.clear();

  // don't map a poisned read pair
  poison_state.set_fragment_end(mapping::util::fragment_end::LEFT);
  // if we have a first read for this record, map it
  bool early_exit_left =
    (records.seq1 == nullptr)
      ? (map_cache_left.clear(), false)
      : mapping::util::map_read(records.seq1, map_cache_left, poison_state,
                                skip_strat);
  if (poison_state.is_poisoned()) {
    return false;
  }

  poison_state.set_fragment_end(mapping::util::fragment_end::RIGHT);
  // if we have a second read for this record, map it
  bool early_exit_right =
    (records.seq2 == nullptr)
      ? (map_cache_right.clear(), false)
      : mapping::util::map_read(records.seq2, map_cache_right, poison_state,
                                skip_strat);

  if (poison_state.is_poisoned()) {
    return false;
  }

  int32_t left_len = (records.seq1 == nullptr)
                       ? 0
                       : static_cast<int32_t>(records.seq1->length());
  int32_t right_len = (records.seq2 == nullptr)
                        ? 0
                        : static_cast<int32_t>(records.seq2->length());

  /*
  std::cerr << "\n\nnew_read: " << *records.seq2 << "\n";
  for (auto& lh : map_cache_left.accepted_hits) {
    std::cerr << "left: " << lh.tid << ", " << lh.pos << " (" << (lh.is_fw ?
  "fw" : "rc") <<
  ")\n";
  }
  for (auto& lh : map_cache_right.accepted_hits) {
    std::cerr << "right: " << lh.tid << ", " << lh.pos << " (" << (lh.is_fw ?
  "fw" : "rc") <<
  ")\n";
  }
  */
  mapping::util::merge_se_mappings(map_cache_left, map_cache_right, left_len,
                                   right_len, map_cache_out);

  return (early_exit_left or early_exit_right);
}

template <typename Protocol, typename SketchHitT, typename PairedEndT>
void do_map(mindex::reference_index &ri,
            fastx_parser::FastxParser<fastx_parser::ReadPair> &parser,
            poison_table &poison_map, const Protocol &p,
            const pesc_sc_options &po, std::atomic<uint64_t> &global_nr,
            std::atomic<uint64_t> &global_nhits,
            std::atomic<uint64_t> &global_npoisoned, pesc_output_info &out_info,
            std::mutex &iomut) {

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

  pufferfish::CanonicalKmerIterator kit_end;

  // set up the poison table
  bool use_poison = !poison_map.empty();
  mapping::util::poison_state_t poison_state;
  if (use_poison) {
    poison_state.ptab = &poison_map;
  }

  poison_state.paired_for_mapping = std::is_same_v<PairedEndT, PairedEndBioSeq>;

  // put these in struct
  size_t num_short_umi{0};
  size_t num_ambig_umi{0};
  (void)num_short_umi;
  (void)num_ambig_umi;

  mapping_cache_info<SketchHitT, piscem::streaming_query<false>> map_cache_left(ri);
  map_cache_left.max_ec_card = po.max_ec_card;
  map_cache_left.max_hit_occ = po.max_hit_occ;
  map_cache_left.max_hit_occ_recover = po.max_hit_occ_recover;
  map_cache_left.max_read_occ = po.max_read_occ;
  map_cache_left.attempt_occ_recover = po.attempt_occ_recover;

  mapping_cache_info<SketchHitT, piscem::streaming_query<false>> map_cache_right(ri);
  map_cache_right.max_ec_card = po.max_ec_card;
  map_cache_right.max_hit_occ = po.max_hit_occ;
  map_cache_right.max_hit_occ_recover = po.max_hit_occ_recover;
  map_cache_right.max_read_occ = po.max_read_occ;
  map_cache_right.attempt_occ_recover = po.attempt_occ_recover;

  mapping_cache_info<SketchHitT, piscem::streaming_query<false>> map_cache_out(ri);
  map_cache_out.max_ec_card = po.max_ec_card;
  map_cache_out.max_hit_occ = po.max_hit_occ;
  map_cache_out.max_hit_occ_recover = po.max_hit_occ_recover;
  map_cache_out.max_read_occ = po.max_read_occ;
  map_cache_out.attempt_occ_recover = po.attempt_occ_recover;

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
    for (auto &record : rg) {
      ++global_nr;
      auto rctr = global_nr.load();
      auto hctr = global_nhits.load();

      if (write_mapping_rate and (rctr % 500000 == 0)) {
        iomut.lock();
        std::cerr << "\rprocessed (" << rctr << ") reads; (" << hctr
                  << ") had mappings.";
        iomut.unlock();
      }

      // first extract the barcode
      std::string *bc =
        protocol.extract_bc(record.first.seq, record.second.seq);
      // if we couldn't get it, don't bother with
      // anything else for this read.
      if (bc == nullptr) {
        continue;
      }

      // correct up to one `N` in the barcode
      auto recovered = single_cell::util::recover_barcode(*bc);
      // if we couldn't correct it with 1 `N`, then skip.
      if (recovered == BarCodeRecovered::NOT_RECOVERED) {
        continue;
      }

      // convert it to a k-mer type
      bc_kmer_t bc_kmer;
      bool bc_ok = bc_kmer.fromChars(*bc);
      if (!bc_ok) {
        continue;
      }

      std::string *umi =
        protocol.extract_umi(record.first.seq, record.second.seq);
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
      AlignableReadSeqs read_seqs = protocol.get_mappable_read_sequences(
        record.first.seq, record.second.seq);

      bool had_early_stop = false;
      // dispatch on the *compile-time determined* paired-endness of this
      // protocol.
      if constexpr (std::is_same_v<PairedEndT, SingleEndBioSeq>) {
        had_early_stop =
          map_se_fragment(read_seqs, poison_state, po.skip_strat,
                          map_cache_left, map_cache_right, map_cache_out);
      } else {
        had_early_stop =
          map_pe_fragment(read_seqs, poison_state, po.skip_strat,
                          map_cache_left, map_cache_right, map_cache_out);
      }
      (void)had_early_stop;
      if (poison_state.is_poisoned()) {
        global_npoisoned++;
        // should we do this here?
        // map_cache_left.clear()
        // map_cache_right.clear()
        // map_cache_out.clear()
      }

      global_nhits += map_cache_out.accepted_hits.empty() ? 0 : 1;
      rad::util::write_to_rad_stream(
        bc_kmer, umi_kmer, map_cache_out.map_type, map_cache_out.accepted_hits,
        map_cache_out.unmapped_bc_map, num_reads_in_chunk, rad_w);

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
  { // make a scope and dump the unmapped barcode counts
    rad_writer ubcw;
    for (auto &kv : map_cache_out.unmapped_bc_map) {
      ubcw << kv.first;
      ubcw << kv.second;
    }
    out_info.unmapped_bc_mutex.lock();
    out_info.unmapped_bc_file << ubcw;
    out_info.unmapped_bc_mutex.unlock();
    ubcw.clear();
  }
}

template <typename Protocol>
void do_map_dispatch(mindex::reference_index &ri,
                     fastx_parser::FastxParser<fastx_parser::ReadPair> &parser,
                     poison_table &poison_map, const Protocol &p,
                     const pesc_sc_options &po,
                     std::atomic<uint64_t> &global_nr,
                     std::atomic<uint64_t> &global_nhits,
                     std::atomic<uint64_t> &global_npoisoned,
                     pesc_output_info &out_info, std::mutex &iomut) {

  const bool bio_paired_end = p.is_bio_paired_end();
  spdlog_piscem::debug("mapping configuration: protocol = {}, structural "
                       "constraints = {}, bio. paired end = {}",
                       p.get_name(), po.enable_structural_constraints,
                       bio_paired_end);
  if (!po.enable_structural_constraints) {
    using SketchHitT = mapping::util::sketch_hit_info_no_struct_constraint;
    if (bio_paired_end) {
      do_map<Protocol, SketchHitT, PairedEndBioSeq>(
        ri, parser, poison_map, p, po, global_nr, global_nhits,
        global_npoisoned, out_info, iomut);
    } else {
      do_map<Protocol, SketchHitT, SingleEndBioSeq>(
        ri, parser, poison_map, p, po, global_nr, global_nhits,
        global_npoisoned, out_info, iomut);
    }
  } else {
    using SketchHitT = mapping::util::sketch_hit_info;
    if (bio_paired_end) {
      do_map<Protocol, SketchHitT, PairedEndBioSeq>(
        ri, parser, poison_map, p, po, global_nr, global_nhits,
        global_npoisoned, out_info, iomut);
    } else {
      do_map<Protocol, SketchHitT, SingleEndBioSeq>(
        ri, parser, poison_map, p, po, global_nr, global_nhits,
        global_npoisoned, out_info, iomut);
    }
  }
}

bool set_geometry(std::string &library_geometry, protocol_t &pt,
                  std::unique_ptr<custom_protocol> &p) {
  if (library_geometry == "chromium_v2") {
    umi_kmer_t::k(10);
    bc_kmer_t::k(16);
    pt = protocol_t::CHROM_V2;
  } else if (library_geometry == "chromium_v2_5p") {
    umi_kmer_t::k(10);
    bc_kmer_t::k(16);
    pt = protocol_t::CHROM_V2_5P;
  } else if (library_geometry == "chromium_v3") {
    umi_kmer_t::k(12);
    bc_kmer_t::k(16);
    pt = protocol_t::CHROM_V3;
  } else if (library_geometry == "chromium_v3_5p") {
    umi_kmer_t::k(12);
    bc_kmer_t::k(16);
    pt = protocol_t::CHROM_V3_5P;
  } else if (library_geometry == "chromium_v4_3p") {
    umi_kmer_t::k(12);
    bc_kmer_t::k(16);
    pt = protocol_t::CHROM_V4_3P;
  } else {
    std::unique_ptr<custom_protocol> opt_cp =
      single_cell::util::parse_custom_geometry(library_geometry);
    if (opt_cp) {
      p.swap(opt_cp);
      umi_kmer_t::k(p->get_umi_len());
      bc_kmer_t::k(p->get_bc_len());
      pt = protocol_t::CUSTOM;
    } else {
      spdlog_piscem::critical(
        "could not parse custom geometry description [{}]", library_geometry);
      return false;
    }
  }
  return true;
}

#ifdef __cplusplus
extern "C" {
#endif
int run_pesc_sc(int argc, char **argv);
#ifdef __cplusplus
}
#endif

int run_pesc_sc(int argc, char **argv) {
  /**
   * piscem single-cell mapper
   **/
  std::ios_base::sync_with_stdio(false);

  std::string skipping_rule;
  pesc_sc_options po;
  { // we need to ensure that no resource within po is collected until the end
    // of this scope.

    // if the user requested to list the geometries, then do that and
    // exit.
    auto list_geometries = [](std::size_t s) {
      (void)s;
      std::cout << "{\n \"supported_geometries\": [\n";
      for (size_t geo_idx = 0; geo_idx < builtin_geometries.size(); ++geo_idx) {
        auto &geo = builtin_geometries[geo_idx];
        std::cout << "  \"" << geo << "\"";
        if (geo_idx < builtin_geometries.size() - 1) {
          std::cout << ", ";
        }
        std::cout << '\n';
      }
      std::cout << "  ]\n}\n";
      std::exit(0);
    };

    CLI::App app{"PESC — single-cell RNA-seq mapper for alevin-fry"};
    app
      .add_flag("--list-geometries", list_geometries,
                "List the known geometries")
      ->trigger_on_parse();
    app.add_option("-i,--index", po.index_basename, "Input index prefix")
      ->required();
    app
      .add_option("-1,--read1", po.left_read_filenames,
                  "Path to list of (comma separated) read 1 files")
      ->required()
      ->delimiter(',');
    app
      .add_option("-2,--read2", po.right_read_filenames,
                  "Path to list of (comma separated) read 2 files")
      ->required()
      ->delimiter(',');
    app
      .add_option("-o,--output", po.output_dirname, "Path to output directory")
      ->required();
    app
      .add_option("-g,--geometry", po.library_geometry,
                  "Geometry of barcode, umi and read")
      ->required();
    app
      .add_option("-t,--threads", po.nthread,
                  "An integer that specifies the number of threads to use")
      ->default_val(16);
    app.add_flag(
      "--no-poison", po.no_poison,
      "Do not filter reads for poison k-mers, even if a poison table "
      "exists for the index");
    app.add_flag("-c,--struct-constraints", po.enable_structural_constraints,
                 "Apply structural constraints when performing mapping");
    app
      .add_option(
        "--skipping-strategy", skipping_rule,
        "Which skipping rule to use for pseudoalignment ({strict, permissive})")
      ->default_val("permissive");
    app.add_flag("--quiet", po.quiet,
                 "Try to be quiet in terms of console output");
    auto ignore_ambig =
      app.add_flag("--ignore-ambig-hits", po.ignore_ambig_hits,
                   "Ignore the existence of highly-frequent hits in mapped "
                   "targets, rather than "
                   "falling back to a simplified strategy to count them.");
    app
      .add_option("--max-ec-card", po.max_ec_card,
                  "Determines the maximum cardinality equivalence class "
                  "(number of (txp, orientation status) pairs) to examine "
                  "if not ignoring highly ambiguous hits.")
      ->excludes(ignore_ambig)
      ->default_val(piscem::defaults::max_ec_card);

    auto agroup = app.add_option_group(
      "advanced options", "provide advanced options to control mapping");
    agroup
      ->add_option(
        "--max-hit-occ", po.max_hit_occ,
        "In the first pass, consider only k-mers having <= --max-hit-occ hits.")
      ->default_val(piscem::defaults::max_hit_occ);
    agroup
      ->add_option(
        "--max-hit-occ-recover", po.max_hit_occ_recover,
        "If all k-mers have > --max-hit-occ hits, then make a second "
        "pass and consider k-mers "
        "having <= --max-hit-occ-recover hits.")
      ->default_val(piscem::defaults::max_hit_occ_recover);
    agroup
      ->add_option("--max-read-occ", po.max_read_occ,
                   "Reads with more than this number of mappings will not have "
                   "their mappings reported.")
      ->default_val(piscem::defaults::max_read_occ);

    /*
     "Note, this can be greater than "
     "--max-hit-occ-recover because of ambiguous hit recovery "
     "and the --max-ec-card parameter.")
    */

    CLI11_PARSE(app, argc, argv);

    spdlog_piscem::drop_all();
    auto logger =
      spdlog_piscem::create<spdlog_piscem::sinks::stderr_color_sink_mt>("");
    logger->set_pattern("%+");
    spdlog_piscem::set_default_logger(logger);

    if (po.quiet) {
      spdlog_piscem::set_level(spdlog_piscem::level::warn);
    }

    spdlog_piscem::info("enable structural constraints : {}",
                        po.enable_structural_constraints);

    po.attempt_occ_recover = (po.max_hit_occ_recover > po.max_hit_occ);

    // start the timer
    auto start_t = std::chrono::high_resolution_clock::now();

    std::optional<mindex::SkippingStrategy> skip_strat_opt =
      mindex::SkippingStrategy::from_string(skipping_rule);
    if (!skip_strat_opt) {
      spdlog_piscem::critical("The skipping strategy must be one of \"strict\" "
                              "or \"permissive\", but \"{}\" was passed in",
                              skipping_rule);
      return 1;
    }
    po.skip_strat = skip_strat_opt.value();

    bool geom_ok = set_geometry(po.library_geometry, po.pt, po.p);
    if (!geom_ok) {
      spdlog_piscem::critical("could not set the library geometry properly.");
      return 1;
    }

    // load a poison map if we had one
    // load a poison map if we had one
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

    // load the main index
    bool attempt_load_ec_map = !po.ignore_ambig_hits;
    mindex::reference_index ri(po.index_basename, attempt_load_ec_map);

    // RAD file path
    ghc::filesystem::path output_path(po.output_dirname);
    ghc::filesystem::create_directories(output_path);

    ghc::filesystem::path rad_file_path = output_path / "map.rad";
    ghc::filesystem::path unmapped_bc_file_path =
      output_path / "unmapped_bc_count.bin";

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
      spdlog_piscem::critical("Could not open {} for writing.",
                              rad_file_path.string());
      throw std::runtime_error("error creating output file.");
    }

    if (!unmapped_bc_file.good()) {
      spdlog_piscem::critical("Could not open {} for writing.",
                              unmapped_bc_file_path.string());
      throw std::runtime_error("error creating output file.");
    }

    pesc_output_info out_info;
    out_info.rad_file = std::move(rad_file);
    out_info.unmapped_bc_file = std::move(unmapped_bc_file);

    size_t bc_length = bc_kmer_t::k();
    size_t umi_length = umi_kmer_t::k();
    size_t chunk_offset =
      rad::util::write_rad_header(ri, bc_length, umi_length, out_info.rad_file);

    std::mutex iomut;

    uint32_t np = 1;
    
    auto num_input_files = po.left_read_filenames.size();
    size_t additional_files = (num_input_files > 1) ? (num_input_files - 1) : 0;

    // start with 1 parsing thread, and one more for every
    // 6 threads, as long as there are additional input files
    // to parse.
    size_t remaining_threads = po.nthread;
    for (size_t i = 0; i < additional_files; ++i) {
      if (remaining_threads >= 6) {
        np += 1;
        po.nthread -= 1;
        remaining_threads -= 6;
      } else {
        break;
      }
    }

    fastx_parser::FastxParser<fastx_parser::ReadPair> rparser(
      po.left_read_filenames, po.right_read_filenames, po.nthread, np);
    rparser.start();

    // set the k-mer size for the
    // CanonicalKmer type.
    CanonicalKmer::k(ri.k());

    std::atomic<uint64_t> global_nr{0};
    std::atomic<uint64_t> global_nh{0};
    std::atomic<uint64_t> global_np{0};
    std::vector<std::thread> workers;
    for (size_t i = 0; i < po.nthread; ++i) {
      switch (po.pt) {
      case protocol_t::CHROM_V2: {
        workers.push_back(
          std::thread([&ri, &rparser, &ptab, &po, &global_nr, &global_nh,
                       &global_np, &iomut, &out_info]() {
            chromium_v2 prot;
            do_map_dispatch<decltype(prot)>(ri, rparser, ptab, prot, po,
                                            global_nr, global_nh, global_np,
                                            out_info, iomut);
          }));
        break;
      }
      case protocol_t::CHROM_V2_5P: {
        workers.push_back(
          std::thread([&ri, &rparser, &ptab, &po, &global_nr, &global_nh,
                       &global_np, &iomut, &out_info]() {
            chromium_v2_5p prot;
            do_map_dispatch<decltype(prot)>(ri, rparser, ptab, prot, po,
                                            global_nr, global_nh, global_np,
                                            out_info, iomut);
          }));
        break;
      }
      case protocol_t::CHROM_V3: {
        workers.push_back(
          std::thread([&ri, &rparser, &ptab, &po, &global_nr, &global_nh,
                       &global_np, &iomut, &out_info]() {
            chromium_v3 prot;
            do_map_dispatch<decltype(prot)>(ri, rparser, ptab, prot, po,
                                            global_nr, global_nh, global_np,
                                            out_info, iomut);
          }));
        break;
      }
      case protocol_t::CHROM_V3_5P: {
        workers.push_back(
          std::thread([&ri, &rparser, &ptab, &po, &global_nr, &global_nh,
                       &global_np, &iomut, &out_info]() {
            chromium_v3_5p prot;
            do_map_dispatch<decltype(prot)>(ri, rparser, ptab, prot, po,
                                            global_nr, global_nh, global_np,
                                            out_info, iomut);
          }));
        break;
      }
      case protocol_t::CHROM_V4_3P: {
        workers.push_back(
          std::thread([&ri, &rparser, &ptab, &po, &global_nr, &global_nh,
                       &global_np, &iomut, &out_info]() {
            chromium_v4_3p prot;
            do_map_dispatch<decltype(prot)>(ri, rparser, ptab, prot, po,
                                            global_nr, global_nh, global_np,
                                            out_info, iomut);
          }));
        break;
      }
      case protocol_t::CUSTOM: {
        workers.push_back(
          std::thread([&ri, &rparser, &ptab, &po, &global_nr, &global_nh,
                       &global_np, &iomut, &out_info]() {
            // static_assert(std::is_same<custom_protocol,
            // decltype(po.p)::element_type>::value, "types are not the same!");
            // for the tempalate parameter below, we can use:
            // custom_protocol
            // decltype(po.p)::element_type
            // std::remove_reference<decltype(*(po.p))>::value
            //
            // but, we CAN NOT USE:
            // decltype(*(po.p)), as this resolves to custom_protocol& because
            // C++ hates reason
            do_map_dispatch<custom_protocol>(ri, rparser, ptab, *(po.p), po,
                                             global_nr, global_nh, global_np,
                                             out_info, iomut);
          }));
        break;
      }
      }
    }

    for (auto &w : workers) {
      w.join();
    }
    rparser.stop();

    spdlog_piscem::info("finished mapping.");

    // rewind to the start of the file and write the number of
    // chunks that we actually produced.
    out_info.rad_file.seekp(chunk_offset);
    uint64_t nc = out_info.num_chunks.load();
    out_info.rad_file.write(reinterpret_cast<char *>(&nc), sizeof(nc));

    out_info.rad_file.close();

    // We want to check if the RAD file stream was written to
    // properly. While we likely would have caught this earlier,
    // it is possible the badbit may not be set until the stream
    // actually flushes (perhaps even at close), so we check here
    // one final time that the status of the stream is as
    // expected. see :
    // https://stackoverflow.com/questions/28342660/error-handling-in-stdofstream-while-writing-data
    if (!out_info.rad_file) {
      spdlog_piscem::critical("The RAD file stream had an invalid status after "
                              "close; so some operation(s) may"
                              "have failed!\nA common cause for this is lack "
                              "of output disk space.\n"
                              "Consequently, the output may be corrupted.\n\n");
      return 1;
    }

    out_info.unmapped_bc_file.close();
    // Same as the RAD file stream above, check to make
    // sure the output is properly written.
    if (!out_info.unmapped_bc_file) {
      spdlog_piscem::critical(
        "The unmapped barcode file stream had an invalid status after "
        "close; so some operation(s) may"
        "have failed!\nA common cause for this is lack "
        "of output disk space.\n"
        "Consequently, the output may be corrupted.\n\n");
      return 1;
    }

    auto end_t = std::chrono::high_resolution_clock::now();
    auto num_sec =
      std::chrono::duration_cast<std::chrono::seconds>(end_t - start_t);
    piscem::meta_info::run_stats rs;
    rs.cmd_line(cmdline);
    rs.mode(piscem::RunMode::scrna);
    rs.num_reads(global_nr.load());
    rs.num_hits(global_nh.load());
    rs.num_poisoned(global_np.load());
    rs.num_seconds(num_sec.count());
    rs.ref_sig_info(ri.ref_sig_info());

    ghc::filesystem::path map_info_file_path = output_path / "map_info.json";
    bool info_ok = piscem::meta_info::write_map_info(rs, map_info_file_path);
    if (!info_ok) {
      spdlog_piscem::critical("failed to write map_info.json file");
    }
  }
  if (po.p) {
    spdlog_piscem::info(
      "used custom geometry with total bc length {} and total umi length {}",
      po.p->get_bc_len(), po.p->get_umi_len());
  }
  return 0;
}
