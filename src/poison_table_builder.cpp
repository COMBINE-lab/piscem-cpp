#include "../include/CanonicalKmerIterator.hpp"
#include "../include/FastxParser.hpp"
#include "../include/Kmer.hpp"
#include "../include/cli11/CLI11.hpp"
#include "../include/ghc/filesystem.hpp"
#include "../include/reference_index.hpp"
#include "../include/util.hpp"
#include "../include/mapping/utils.hpp"
#include "../include/parallel_hashmap/phmap.h"
#include "../include/parallel_hashmap/phmap_dump.h"
#include "spdlog_piscem/spdlog.h"

#include <iostream>
#include <vector>
#include <numeric>
#include <cstdio>

#ifdef __cplusplus
extern "C" {
#endif
int run_build_poison_table(int argc, char** argv);
#ifdef __cplusplus
}
#endif

using poison_map_t = phmap::flat_hash_set<uint64_t, sshash::RobinHoodHash>;
using pufferfish::CanonicalKmerIterator;

struct poison_kmer_state_strict {
public:
  void reset() {
    predecessor_present = false;
    predecessor_kmer = 0;
  }

  // given the current poison k-mer state and the current k-mer iterator,
  // determine if the current or former k-mer is poison.
  inline bool inspect_and_update(CanonicalKmerIterator& kit, mindex::reference_index& ri, 
                          sshash::streaming_query_canonical_parsing& cache, 
                          poison_map_t& poison_kmers) {
    bool added_poison = false;
    // current canonical k-mer
    auto kmer = kit->first;
    auto phits = ri.query(kit, cache);

    // The condition for being poison is that either
    // (i) : The predecessor was present but the current k-mer 
    // is absent. In this case the current k-mer is poison.
    // (ii) : The predecessor was absent but the current k-mer 
    // is present. In this case, the predecessor k-mer is poison.

    bool present = !phits.empty();
    if (predecessor_present and !present) {
      // (case i) -> current k-mer is poison
      auto it = poison_kmers.insert(kmer.getCanonicalWord());
      added_poison = it.second;
    } else if (!predecessor_present and present) {
      // (case ii) -> predecessor k-mer is poison
      auto it = poison_kmers.insert(predecessor_kmer);
      added_poison = it.second;
    } 

    predecessor_kmer = kmer.getCanonicalWord();
    predecessor_present = present;
    return added_poison;
  }

private:
  bool predecessor_present{false};
  uint64_t predecessor_kmer{0};
};

struct poison_kmer_state_permissive {
public:
  void reset() { }

  // given the current poison k-mer state and the current k-mer iterator,
  // determine if the current or former k-mer is poison.
  inline bool inspect_and_update(CanonicalKmerIterator& kit, mindex::reference_index& ri, 
                          sshash::streaming_query_canonical_parsing& cache, 
                          poison_map_t& poison_kmers) {
    auto* dict = ri.get_dict();
    // current canonical k-mer
    auto kmer = kit->first;
    
    // given the current kmer, check if it is in the index or not
    // if it *is* in the index, then it cannot be a poison kmer
    auto phits = ri.query(kit, cache);
    bool present = !phits.empty();
    
    if (present) { return false; }

    // if we got here, the current k-mer is not in the index. Now
    // we can systematically check all possible neighbors.
    const std::array<char, 4> nucs = {'A', 'C', 'G', 'T'};
    for (const auto n : nucs) {
      CanonicalKmer neighbor = kmer;
      neighbor.shiftFw(n);
      bool neighbor_present = dict->is_member_uint64(neighbor.getCanonicalWord());
      if (neighbor_present) {
        poison_kmers.insert(kmer.getCanonicalWord());
        return true;
      }
    }

    // append to the other side
    for (const auto n : nucs) {
      CanonicalKmer neighbor = kmer;
      neighbor.shiftBw(n);
      bool neighbor_present = dict->is_member_uint64(neighbor.getCanonicalWord());
      if (neighbor_present) {
        auto it = poison_kmers.insert(kmer.getCanonicalWord());
        return it.second;
      }
    }
    return false;
  }
};

template <typename poison_state_t>
void find_poison_kmers(mindex::reference_index& ri, 
                      fastx_parser::FastxParser<fastx_parser::ReadSeq>& rparser,
                      std::atomic<uint64_t>& global_nk,
                      poison_map_t& poison_kmers) {
  pufferfish::CanonicalKmerIterator kit_end;
  sshash::streaming_query_canonical_parsing cache(ri.get_dict());
  
  poison_state_t pstate;
  pstate.reset();
  
  // Get the read group by which this thread will
  // communicate with the parser (*once per-thread*)
  auto rg = rparser.getReadGroup();
  while (rparser.refill(rg)) {
    // Here, rg will contain a chunk of read pairs
    // we can process.
    for (auto& record : rg) {
      spdlog_piscem::info("processing {}", record.name);
      pstate.reset();
      cache.start();
      
      pufferfish::CanonicalKmerIterator kit(record.seq);
      while (kit != kit_end) {

        bool inserted_locally = pstate.inspect_and_update(kit, ri, cache, poison_kmers);
        (void) inserted_locally;
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

int run_build_poison_table(int argc, char* argv[]) {

	std::ios_base::sync_with_stdio(false);
  spdlog_piscem::set_level(spdlog_piscem::level::info);

  decoy_opts po;
  CLI::App app{"Extract poison k-mers."};
  app.add_option("-i,--index", po.index_basename, "input index prefix")->required();
  app.add_option("-d,--decoys", po.decoy_seq_paths, "path to list of decoy files")
    ->required()
    ->delimiter(',');
  app.add_option("-m,--method", po.poison_method, "method to determine poison k-mers (one of {edge, node})")
    ->default_val("edge");
  app.add_flag("--overwrite", po.overwrite, "overwrite the poison file if it exists");
  app.add_option("-t,--threads", po.nthreads,
                 "An integer that specifies the number of threads to use")
    ->default_val(16);
  app.add_flag("--quiet", po.quiet, "try to be quiet in terms of console output");
  CLI11_PARSE(app, argc, argv);

  if ( !(po.poison_method == "edge" or po.poison_method == "node") ) {
    spdlog_piscem::critical("The -m/--method flag must be one of \"edge\" or \"node\", but {} was provided", 
                            po.poison_method);
    return 1;
  }

  po.output_file = po.index_basename + ".poison";
  if (ghc::filesystem::exists(po.output_file)) {
    if (!po.overwrite) {
      spdlog_piscem::critical("The poison file {} already exists. If you wish to overwite it, please use the --overwrite flag.",
                              po.output_file);
      return 1;
    }
  }

	mindex::reference_index ri(po.index_basename);
  // set the k-mer size for the
  // CanonicalKmer type.
  CanonicalKmer::k(ri.k());

  size_t np = 1;
  std::vector<std::thread> workers;
  if ((po.decoy_seq_paths.size() > 1) and (po.nthreads >= 6)) {
    np += 1;
    po.nthreads -= 1;
  }

  std::vector<poison_map_t> poison_kmer_maps;
  poison_kmer_maps.resize(po.nthreads);

  std::atomic<uint64_t> global_nk{0};
  std::atomic<uint64_t> global_np{0};

  fastx_parser::FastxParser<fastx_parser::ReadSeq> rparser(po.decoy_seq_paths, po.nthreads,
                                                           np, 1);
  rparser.start();

  if (po.poison_method == "edge") {
    for (size_t i = 0; i < po.nthreads; ++i) {
      workers.push_back(
        std::thread([i, &poison_kmer_maps, &ri, &rparser, &global_nk]() {
          auto& pkmap = poison_kmer_maps[i];
          find_poison_kmers<poison_kmer_state_strict>(ri, rparser, global_nk, pkmap);
        }));
    }
  } else {
    for (size_t i = 0; i < po.nthreads; ++i) {
      workers.push_back(
        std::thread([i, &poison_kmer_maps, &ri, &rparser, &global_nk]() {
          auto& pkmap = poison_kmer_maps[i];
          find_poison_kmers<poison_kmer_state_permissive>(ri, rparser, global_nk, pkmap);
        }));
    }
  }

  for (auto& w : workers) { w.join(); }
  rparser.stop();

  poison_map_t poison_map;
  while (!poison_kmer_maps.empty()) {
    auto& local_map = poison_kmer_maps.back();

    for (auto elem : local_map) {
      poison_map.insert(elem);
    }
    poison_kmer_maps.pop_back();
  }

  spdlog_piscem::info("FINAL: Examined {} total decoy k-mers, recorded {} poison k-mers.", global_nk, poison_map.size()); 

  phmap::BinaryOutputArchive ar_out(po.output_file.c_str());
  poison_map.phmap_dump(ar_out);

  return 0;
}
