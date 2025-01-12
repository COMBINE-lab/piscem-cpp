#include "../include/CanonicalKmerIterator.hpp"
#include "../include/FastxParser.hpp"
#include "../include/Kmer.hpp"
#include "../include/cli11/CLI11.hpp"
#include "../include/ghc/filesystem.hpp"
#include "../include/reference_index.hpp"
#include "../include/util_piscem.hpp"
#include "../external/sshash/include/util.hpp"
#include "../include/mapping/utils.hpp"
#include "../include/parallel_hashmap/phmap.h"
#include "../include/parallel_hashmap/phmap_dump.h"
#include "spdlog_piscem/spdlog.h"

#include <iostream>
#include <vector>
#include <numeric>
#include <cstdio>

using poison_map_t = phmap::flat_hash_map<uint64_t, uint64_t, sshash::RobinHoodHash>;

void filter_poison_reads(poison_map_t& poison_map, 
                         std::vector<std::string>& read_filenames,
                         const std::string& output_file) {
  (void) output_file;
  fastx_parser::FastxParser<fastx_parser::ReadSeq> rparser(read_filenames, 1, 1);
  rparser.start();

  pufferfish::CanonicalKmerIterator kit_end;
  
  size_t reads_processed = {0};
  size_t poison_reads = {0};
  auto pmap_end = poison_map.end();

  // Get the read group by which this thread will
  // communicate with the parser (*once per-thread*)
  auto rg = rparser.getReadGroup();
  while (rparser.refill(rg)) {
    // Here, rg will contain a chunk of read pairs
    // we can process.
    for (auto& record : rg) {
      
      pufferfish::CanonicalKmerIterator kit(record.seq);
      ++reads_processed;
      while (kit != kit_end) {
        // current canonical k-mer
        if (kit->first.is_homopolymer()) { 
          kit++;
          continue; 
        }

        auto pmap_it = poison_map.find(kit->first.getCanonicalWord());
        if (pmap_it != pmap_end) {
          poison_reads++;
          break;;
        }
        ++kit;
      }
      if (reads_processed > 1 and reads_processed % 100000 == 0) {
        spdlog_piscem::info("processed {} reads: total poison reads = {}", reads_processed, poison_reads); 
      }
    }
  }

  spdlog_piscem::info("FINAL: processed {} reads: total poison reads = {}", reads_processed, poison_reads); 
  rparser.stop();
}

struct decoy_opts {
  std::string poison_map_name;
  std::vector<std::string> read_filenames;
  std::string output_file;
  size_t nthreads;
  bool quiet{false};
};

int main(int argc, char* argv[]) {

	std::ios_base::sync_with_stdio(false);
  spdlog_piscem::set_level(spdlog_piscem::level::info);

  decoy_opts po;
  CLI::App app{"Filter poison reads."};
  app.add_option("-i,--index", po.poison_map_name, "poison index file")->required();
  app.add_option("-r,--reads", po.read_filenames, "path to list of read_files")
    ->required()
    ->delimiter(',');
  app.add_option("-o,--output", po.output_file, "path to output file")->required();
  app.add_option("-t,--threads", po.nthreads,
                 "An integer that specifies the number of threads to use")
    ->default_val(16);
  app.add_flag("--quiet", po.quiet, "try to be quiet in terms of console output");
  CLI11_PARSE(app, argc, argv);


  spdlog_piscem::info("Loading poison k-mer map...");
  poison_map_t poison_map;
  phmap::BinaryInputArchive ar_in(po.poison_map_name.c_str());
  poison_map.phmap_load(ar_in);
  spdlog_piscem::info("done");

  // set the k-mer size for the
  // CanonicalKmer type.
  CanonicalKmer::k(31);


  filter_poison_reads(poison_map, po.read_filenames, po.output_file);
}
