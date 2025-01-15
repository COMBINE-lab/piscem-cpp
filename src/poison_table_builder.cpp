#include "../external/sshash/include/util.hpp"
#include "../include/CanonicalKmerIterator.hpp"
#include "../include/FastxParser.hpp"
#include "../include/Kmer.hpp"
#include "../include/cli11/CLI11.hpp"
#include "../include/ghc/filesystem.hpp"
#include "../include/mapping/utils.hpp"
#include "../include/parallel_hashmap/phmap.h"
#include "../include/parallel_hashmap/phmap_dump.h"
#include "../include/poison_table.hpp"
#include "../include/reference_index.hpp"
#include "../include/spdlog_piscem/sinks/stdout_color_sinks.h"
#include "../include/spdlog_piscem/spdlog.h"
#include "../include/streaming_query.hpp"
#include "../include/util_piscem.hpp"

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <iterator>
#include <numeric>
#include <vector>

#ifdef __cplusplus
extern "C" {
#endif
int run_build_poison_table(int argc, char **argv);
#ifdef __cplusplus
}
#endif

using pufferfish::CanonicalKmerIterator;
using sshash::labeled_poison_occ_t;
using sshash::poison_occ_t;
using poison_map_t =
  phmap::flat_hash_map<uint64_t, uint64_t, sshash::RobinHoodHash>;

struct poison_kmer_state_strict {
public:
  void reset() {
    first = true;
    predecessor_present = false;
    predecessor_kmer = 0;
  }

  // given the current poison k-mer state and the current k-mer iterator,
  // determine if the current or former k-mer is poison.
  inline bool
  inspect_and_update(CanonicalKmerIterator &kit, mindex::reference_index &ri,
                     piscem::streaming_query<false> &cache,
                     std::vector<labeled_poison_occ_t> &poison_occs) {
    constexpr int64_t zero = 0;
    bool added_poison = false;
    // current canonical k-mer
    auto kmer = kit->first;
    auto phits = ri.query(kit, cache);
    uint32_t k = ri.k();

    // The condition for being poison is that either
    // (i) : The predecessor was present but the current k-mer
    // is absent. In this case the current k-mer is poison.
    // (ii) : The predecessor was absent but the current k-mer
    // is present. In this case, the predecessor k-mer is poison.
    bool present = !phits.empty();
    if (!first) {
      if (predecessor_present and !present) {
        // (case i) -> current k-mer is poison
        // If the current k-mer is missing, but the previous one
        // was present:
        // check if the previous k-mer hit the contig in the forward orientation
        bool isFw = predecessor_phits.hit_fw_on_contig();
        // if it was forward, then the current k-mer is at the position one
        // greater on the contig, otherwise it's at the position one lesser.
        int offset = isFw ? 1 : -1;
        int64_t cpos = predecessor_phits.contig_pos() + offset;
        // we can't go below the start the the contig or above the
        // last k-mer position
        cpos = std::max(zero, cpos);
        cpos = std::min(
          cpos, static_cast<int64_t>(predecessor_phits.contig_len() - k));
        // push the poison k-mer occurrence.
        poison_occs.push_back({kmer.getCanonicalWord(),
                               predecessor_phits.contig_id(),
                               static_cast<uint32_t>(cpos)});
      } else if (!predecessor_present and present) {
        // (case ii) -> predecessor k-mer is poison
        // If the current k-mer is present, but the previous one
        // was missing:
        // check if the current k-mer hit the contig in the forward orientation
        bool isFw = phits.hit_fw_on_contig();
        // if it was forward, then the previous k-mer is at the position one
        // lesser on the contig, otherwise it's at the position one greater.
        int offset = isFw ? -1 : 1;
        int64_t cpos = phits.contig_pos() + offset;
        // we can't go below the start the the contig or above the
        // last k-mer position
        cpos = std::max(zero, cpos);
        cpos = std::min(cpos, static_cast<int64_t>(phits.contig_len() - k));
        poison_occs.push_back(
          {predecessor_kmer, phits.contig_id(), static_cast<uint32_t>(cpos)});
      }
    } else {
      first = false;
    }
    // the current hit becomes the new predecessor
    predecessor_kmer = kmer.getCanonicalWord();
    predecessor_present = present;
    predecessor_phits = phits;

    return added_poison;
  }

private:
  bool first{true};
  bool predecessor_present{false};
  uint64_t predecessor_kmer{0};
  projected_hits predecessor_phits;
};

template <typename poison_state_t>
void find_poison_kmers(
  mindex::reference_index &ri,
  fastx_parser::FastxParser<fastx_parser::ReadSeq> &rparser,
  std::atomic<uint64_t> &global_nk,
  std::vector<labeled_poison_occ_t> &poison_kmer_occs) {
  pufferfish::CanonicalKmerIterator kit_end;

  piscem::streaming_query<false> cache(ri.get_dict());

  poison_state_t pstate;
  pstate.reset();

  // Get the read group by which this thread will
  // communicate with the parser (*once per-thread*)
  auto rg = rparser.getReadGroup();
  while (rparser.refill(rg)) {
    // Here, rg will contain a chunk of read pairs
    // we can process.
    for (auto &record : rg) {
      spdlog_piscem::info("processing {}", record.name);
      pstate.reset();
      cache.reset_state();

      pufferfish::CanonicalKmerIterator kit(record.seq);
      while (kit != kit_end) {
        bool inserted_locally =
          pstate.inspect_and_update(kit, ri, cache, poison_kmer_occs);
        (void)inserted_locally;
        ++kit;
        ++global_nk;
      }
      spdlog_piscem::info("finished processing {}", record.name);
    }
  }
}

struct decoy_opts {
  std::string index_basename;
  std::vector<std::string> decoy_seq_paths;
  std::string output_file;
  std::string poison_method;
  size_t nthreads;
  bool overwrite{false};
  bool quiet{false};
};

int run_build_poison_table(int argc, char *argv[]) {
  std::ios_base::sync_with_stdio(false);

  decoy_opts po;
  CLI::App app{"Extract poison k-mers."};
  app.add_option("-i,--index", po.index_basename, "input index prefix")
    ->required();
  app
    .add_option("-d,--decoys", po.decoy_seq_paths,
                "path to (\',\' separated) list of decoy files")
    ->required()
    ->delimiter(',');
  app
    .add_option("-m,--method", po.poison_method,
                "method to determine poison k-mers (one of {edge, node})")
    ->default_val("edge");
  app.add_flag("--overwrite", po.overwrite,
               "overwrite the poison file if it exists");
  app
    .add_option("-t,--threads", po.nthreads,
                "An integer that specifies the number of threads to use")
    ->default_val(16);
  app.add_flag("--quiet", po.quiet,
               "try to be quiet in terms of console output");
  CLI11_PARSE(app, argc, argv);

  spdlog_piscem::drop_all();
  auto logger =
    spdlog_piscem::create<spdlog_piscem::sinks::stderr_color_sink_mt>("");
  logger->set_pattern("%+");
  spdlog_piscem::set_default_logger(logger);

  if (po.quiet) {
    spdlog_piscem::set_level(spdlog_piscem::level::warn);
  }

  if (!(po.poison_method == "edge" or po.poison_method == "node")) {
    spdlog_piscem::critical("The -m/--method flag must be one of \"edge\" or "
                            "\"node\", but {} was provided",
                            po.poison_method);
    return 1;
  }

  po.output_file = po.index_basename + ".poison";
  if (ghc::filesystem::exists(po.output_file)) {
    if (!po.overwrite) {
      spdlog_piscem::critical(
        "The poison file {} already exists. If you wish to overwite it, please "
        "use the --overwrite flag.",
        po.output_file);
      return 1;
    }
  }

  size_t np = 1;
  std::vector<std::thread> workers;
  if ((po.decoy_seq_paths.size() > 1) and (po.nthreads >= 6)) {
    np += 1;
    po.nthreads -= 1;
  }

  std::vector<std::vector<labeled_poison_occ_t>> poison_kmer_occs;
  poison_kmer_occs.resize(po.nthreads);

  std::atomic<uint64_t> global_nk{0};

  {
    fastx_parser::FastxParser<fastx_parser::ReadSeq> rparser(
      po.decoy_seq_paths, po.nthreads, np, 1);
    rparser.start();
    mindex::reference_index ri(po.index_basename);
    CanonicalKmer::k(ri.k());

    if (po.poison_method == "edge") {
      // set the k-mer size for the
      // CanonicalKmer type.
      for (size_t i = 0; i < po.nthreads; ++i) {
        workers.push_back(
          std::thread([i, &poison_kmer_occs, &ri, &rparser, &global_nk]() {
            auto &pkoccs = poison_kmer_occs[i];
            find_poison_kmers<poison_kmer_state_strict>(ri, rparser, global_nk,
                                                        pkoccs);
          }));
      }
    } else {
      spdlog_piscem::warn("node mode is not yet implemented.");
      rparser.stop();
      return 1;
      /*
    for (size_t i = 0; i < po.nthreads; ++i) {
      workers.push_back(
        std::thread([i, &poison_kmer_occs, &ri, &rparser, &global_nk]() {
          auto& pkoccs = poison_kmer_occs[i];
          find_poison_kmers<poison_kmer_state_permissive>(ri, rparser,
    global_nk, pkoccs);
        }));
    }
    */
    }

    for (auto &w : workers) {
      w.join();
    }
    rparser.stop();
  } // scope to allow deallocating the reference index, since we no logner need
    // it.

  // combine the poison vectors
  size_t total_occs = 0;
  std::vector<std::pair<size_t, size_t>> occs_by_size;
  occs_by_size.reserve(poison_kmer_occs.size());
  for (size_t i = 0; i < poison_kmer_occs.size(); ++i) {
    size_t occ_size = poison_kmer_occs[i].size();
    occs_by_size.push_back(std::make_pair(occ_size, i));
    total_occs += occ_size;
  }
  std::sort(occs_by_size.begin(), occs_by_size.end());

  // start with the smallest one, and add to it.
  std::vector<labeled_poison_occ_t> poison_occs;
  poison_occs.reserve(total_occs);

  // when we are done with a vector we can clear it and
  // reclaim the memory.
  for (auto &si : occs_by_size) {
    size_t i = si.second;
    poison_occs.insert(poison_occs.end(), poison_kmer_occs[i].begin(),
                       poison_kmer_occs[i].end());
    poison_kmer_occs[i].clear();
    poison_kmer_occs[i].shrink_to_fit();
  }

  poison_table ptab;
  if (!ptab.build_from_occs(poison_occs)) {
    spdlog_piscem::critical(
      "could not build poision table from poison k-mer occurrences.");
    return 1;
  }
  if (!ptab.save_to_file(po.output_file, global_nk.load())) {
    spdlog_piscem::critical("could not succesfully save poison table to file.");
    return 1;
  }

  /*
  std::ofstream ofile("poison_kmers.fa");
  size_t pnum{0};
  CanonicalKmer kmer;
  for (auto& km : ptab.poison_map()) {
    kmer.fromNum(km.first);
    ofile << ">"
          << "poison." << pnum << "\n"
          << kmer.to_str() << "\n";
    ++pnum;
  }
  ofile.close();
  */

  return 0;
}
