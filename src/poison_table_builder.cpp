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

using pufferfish::CanonicalKmerIterator;
using sshash::poison_occ_t;
using poison_map_t = phmap::flat_hash_map<uint64_t, uint64_t, sshash::RobinHoodHash>;

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
                          std::vector<poison_occ_t>& poison_occs,
                          poison_map_t& poison_kmers) {
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
    if (predecessor_present and !present) {
      // (case i) -> current k-mer is poison
      //auto it = poison_kmers.insert({kmer.getCanonicalWord(), 1});
      //added_poison = it.second;
      int64_t cpos = phits.contig_pos();
      cpos = std::min(cpos + 1, static_cast<int64_t>(phits.contig_len() - k));
      poison_occs.push_back({kmer.getCanonicalWord(), phits.contig_id(), static_cast<uint32_t>(cpos)});
    } else if (!predecessor_present and present) {
      // (case ii) -> predecessor k-mer is poison
      //auto it = poison_kmers.insert({predecessor_kmer, 1});
      //added_poison = it.second;
      int64_t cpos = phits.contig_pos();
      cpos = std::max(static_cast<int64_t>(0), cpos - 1);
      poison_occs.push_back({predecessor_kmer, phits.contig_id(), static_cast<uint32_t>(cpos)});
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
                          std::vector<poison_occ_t>& poison_occs,
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
        poison_kmers.insert({kmer.getCanonicalWord(), 1});
        return true;
      }
    }

    // append to the other side
    for (const auto n : nucs) {
      CanonicalKmer neighbor = kmer;
      neighbor.shiftBw(n);
      bool neighbor_present = dict->is_member_uint64(neighbor.getCanonicalWord());
      if (neighbor_present) {
        auto it = poison_kmers.insert({kmer.getCanonicalWord(), 1});
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
                      std::vector<poison_occ_t>& poison_kmer_occs,
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

        bool inserted_locally = pstate.inspect_and_update(kit, ri, cache, poison_kmer_occs, poison_kmers);
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
  std::vector<std::vector<poison_occ_t>> poison_kmer_occs;

  poison_kmer_maps.resize(po.nthreads);
  poison_kmer_occs.resize(po.nthreads);

  std::atomic<uint64_t> global_nk{0};
  std::atomic<uint64_t> global_np{0};

  fastx_parser::FastxParser<fastx_parser::ReadSeq> rparser(po.decoy_seq_paths, po.nthreads,
                                                           np, 1);
  rparser.start();

  if (po.poison_method == "edge") {
    for (size_t i = 0; i < po.nthreads; ++i) {
      workers.push_back(
        std::thread([i, &poison_kmer_maps, &poison_kmer_occs, &ri, &rparser, &global_nk]() {
          auto& pkmap = poison_kmer_maps[i];
          auto& pkoccs = poison_kmer_occs[i];
          find_poison_kmers<poison_kmer_state_strict>(ri, rparser, global_nk, pkoccs, pkmap);
        }));
    }
  } else {
    for (size_t i = 0; i < po.nthreads; ++i) {
      workers.push_back(
        std::thread([i, &poison_kmer_maps, &poison_kmer_occs, &ri, &rparser, &global_nk]() {
          auto& pkmap = poison_kmer_maps[i];
          auto& pkoccs = poison_kmer_occs[i];
          find_poison_kmers<poison_kmer_state_permissive>(ri, rparser, global_nk, pkoccs, pkmap);
        }));
    }
  }

  for (auto& w : workers) { w.join(); }
  rparser.stop();

  // combine the vectors
  size_t smallest_idx = 0;
  size_t smallest = std::numeric_limits<size_t>::max();
  for (size_t i = 0; i < poison_kmer_occs.size(); ++i) {
    if (poison_kmer_occs[i].size() < smallest) {
      smallest_idx = i;
    }
  }

  std::vector<poison_occ_t> poison_occs(
    poison_kmer_occs[smallest_idx].begin(), 
    poison_kmer_occs[smallest_idx].end());

  poison_kmer_occs[smallest_idx].clear();
  poison_kmer_occs[smallest_idx].shrink_to_fit();
  for (size_t i = 0; i < poison_kmer_occs.size(); ++i) {
    poison_occs.insert(poison_occs.end(), poison_kmer_occs[i].begin(), poison_kmer_occs[i].end());
    poison_kmer_occs[i].clear();
    poison_kmer_occs[i].shrink_to_fit();
  }

  // we will build a map from each k-mer to the list of 
  // unitigs and positions where it occurs.  Therefore, here 
  // we want to sort by k-mer, then by unitig, then position
  std::sort(poison_occs.begin(), poison_occs.end(),
            [](const poison_occ_t& a, const poison_occ_t& b) -> bool {
              if (a.canonical_kmer == b.canonical_kmer) {
                if (a.unitig_id == b.unitig_id) {
                  return a.unitig_pos < b.unitig_pos;
                } else {
                  return a.unitig_id < b.unitig_id;
                }
              } else {
                return a.canonical_kmer < b.canonical_kmer;
              }
            });
  // remove duplicates
  poison_occs.erase( std::unique(poison_occs.begin(), poison_occs.end()), poison_occs.end());
  spdlog_piscem::info("Total number of distinct poison k-mer occs: {}", poison_occs.size());

  // build the overall map
  poison_map_t poison_map;
  std::vector<uint64_t> offsets;
  offsets.reserve(poison_occs.size()+1);

  size_t max_range = 0;
  auto occ_it = poison_occs.begin();
  auto range_start_it = occ_it;
  offsets.push_back(0);
  poison_map[occ_it->canonical_kmer] = 0;
  while (occ_it != poison_occs.end()) {
    // we started a new range, push back the starting point of 
    // the this range.
    if (occ_it->canonical_kmer != range_start_it->canonical_kmer) {
      size_t range_len = std::distance(range_start_it, occ_it);
      max_range = std::max(range_len, max_range);

      uint64_t dist_from_start = std::distance(poison_occs.begin(), occ_it);
      offsets.push_back(dist_from_start);
      poison_map[occ_it->canonical_kmer] = dist_from_start;
      range_start_it = occ_it;
    }
    occ_it++;
  }
  // don't forget the last one
  uint64_t dist_from_start = std::distance(poison_occs.begin(), occ_it);
  if (!offsets.empty() and (dist_from_start != offsets.back())) {
    offsets.push_back(dist_from_start); 
  }
  offsets.shrink_to_fit();

  spdlog_piscem::info("The most frequently occuring poison k-mer appeared in {} distinct unitig positions.", max_range);

  {
    std::string poc_filename = po.output_file + "_occs";
    std::ofstream poc_file(poc_filename, std::ios::binary);
    if (!poc_file.good()) {
      spdlog_piscem::critical("could not open occ output file {}", poc_filename);
    }

    size_t s = offsets.size();
    poc_file.write(reinterpret_cast<char*>(&s), sizeof(s));
    poc_file.write(reinterpret_cast<char*>(offsets.data()), s * sizeof(uint64_t));

    s = poison_occs.size();
    poc_file.write(reinterpret_cast<char*>(&s), sizeof(s));
    poc_file.write(reinterpret_cast<char*>(poison_occs.data()), s * sizeof(poison_occ_t));
  }

  spdlog_piscem::info("FINAL: Examined {} total decoy k-mers, recorded {} poison k-mers.", global_nk, poison_map.size()); 

  phmap::BinaryOutputArchive ar_out(po.output_file.c_str());
  poison_map.phmap_dump(ar_out);

  std::ofstream ofile("poison_kmers.fa");
  size_t pnum{0};
  CanonicalKmer kmer;
  for (auto& km : poison_map) {
    kmer.fromNum(km.first);
    ofile << ">"
          << "poison." << pnum << "\n"
          << kmer.to_str() << "\n";
    ++pnum;
  }
  ofile.close();

  return 0;
}