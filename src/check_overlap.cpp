#include "../external/sshash/external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "../include/ghc/filesystem.hpp"
#include "../include/sc/util.hpp"
#include "../include/spdlog_piscem/spdlog.h"
#include "../include/spdlog_piscem/sinks/stdout_color_sinks.h"
#include "../include/meta_info.hpp"

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
#include <fstream>

using namespace std;

namespace check_overlap {
enum class protocol_t : uint8_t { CHROM_V2, CHROM_V3, CUSTOM };

enum MateTypeOverlap {
  doveTail,
  overlap,
  noOverlap,
  both
};

std::ostream &operator << ( std::ostream& strm, MateTypeOverlap tt )
{
   const std::string nameTT[] = { "doveTail", "overlap", "noOverlap", "both"};
   return strm << nameTT[tt];
}

struct MateOverlap {
  int frag_length;
  MateTypeOverlap ov_type;
  std::string frag;
  MateOverlap() {
    frag_length=-1;
    ov_type=noOverlap;
    frag="";
  }
};

std::string reverseComplement(std::string sequence) {
    std::string reverse = sequence;
    std::reverse(reverse.begin(), reverse.end()); // reverse the sequence
    for (uint32_t i = 0; i < reverse.size(); i++) {
        if (reverse[i] == 'A') reverse[i] = 'T';
        else if (reverse[i] == 'T') reverse[i] = 'A';
        else if (reverse[i] == 'C') reverse[i] = 'G';
        else if (reverse[i] == 'G') reverse[i] = 'C';
    }
    return reverse;
}

// modified from https://github.com/haowenz/chromap/blob/29bc02d02671ebfb6f71c4b24cafb7548c2ca901/src/chromap.cc#L109
inline std::string getOverlap(std::string& seq1, std::string& seq2, bool dovetail, const int32_t min_overlap_length, const int32_t error_threshold_for_merging) {

    const uint32_t raw_read1_length = seq1.length();
    const uint32_t raw_read2_length = seq2.length();
    
    std::string raw_negative_read1 = reverseComplement(seq1);
    std::string raw_negative_read2 = reverseComplement(seq2);

    std::string &read1 = raw_read1_length <= raw_read2_length ? seq1 : seq2;
    std::string negative_read2 = raw_read1_length <= raw_read2_length
                                          ? raw_negative_read2 : raw_negative_read1;
    std::string frag_seq="";
    bool is_merged = false;
    
    uint32_t overlap_length = 0;
    const int32_t seed_length = min_overlap_length / 2;
    
    std::string suffix_read = dovetail ? negative_read2 : read1;
    std::string prefix_read = dovetail ? read1 : negative_read2;
    const int32_t suff_length = suffix_read.length();
    const int32_t pref_length = prefix_read.length();
    // std::cout << "suffix " << suffix_read << std::endl;
    // std::cout << "prefix " << prefix_read << std::endl;  
    for (int si = 0; si < error_threshold_for_merging + 1; ++si) {
      int32_t seed_start_position =
          suffix_read.find(prefix_read.c_str() + size_t(si*seed_length), 0, size_t(seed_length));
      // std::cout << "seed_pos " << seed_start_position << std::endl;
      while (size_t(seed_start_position) != std::string::npos) {
        const bool before_seed_is_enough_long =
          seed_start_position >= si * seed_length;
        const bool overlap_is_enough_long =
          (suff_length - seed_start_position + seed_length * si) >=
          min_overlap_length;

        if (!before_seed_is_enough_long || !overlap_is_enough_long) {
          seed_start_position = suffix_read.find(
              prefix_read.c_str() + size_t(si * seed_length),size_t(seed_start_position + 1), size_t(seed_length));
          continue;
        }

        bool can_merge = true;
        int num_errors = 0;

        // The bases before the seed.
        for (int i = 0; i < seed_length * si; ++i) {
          if (suffix_read[seed_start_position - si * seed_length + i] !=
              prefix_read[i]) {
            ++num_errors;
            return frag_seq;
          }
          if (num_errors > error_threshold_for_merging) {
            can_merge = false;
            break;
          }
        }
        // std::cout << "num_errors " << num_errors << std::endl;
        // The bases after the seed.
        for (int32_t i = seed_length; i + seed_start_position < suff_length &&
                                      si * seed_length + i < pref_length;
            ++i) {
          // std::cout << "suff " << suffix_read[seed_start_position + i] << " pref " << prefix_read[si * seed_length + i] << std::endl;
          if (suffix_read[seed_start_position + i] !=
              prefix_read[si * seed_length + i]) {
            ++num_errors;
          }
          if (num_errors > error_threshold_for_merging) {
            can_merge = false;
            break;
          }
        }
        // std::cout << "num_errors " << num_errors << std::endl;
        if (can_merge) {
          is_merged = true;
          overlap_length =
              std::min(suff_length - seed_start_position + si * seed_length, pref_length);
          break;
        }

        seed_start_position = suffix_read.find(
            prefix_read.c_str() + size_t(si * seed_length), size_t(seed_start_position + 1), size_t(seed_length));
    }
      if (is_merged) {
          break;
      }
    }
    // std::cout << "overlap length " << overlap_length << std::endl;
    if (overlap_length != 0) {
      if (dovetail) {
        frag_seq = prefix_read.substr(0, overlap_length);
      } else {
        frag_seq = suffix_read + prefix_read.substr(overlap_length, overlap_length-pref_length);
      }
      
    }
  return frag_seq;
}

inline void findOverlapBetweenPairedEndReads(std::string &seq1, std::string &seq2, MateOverlap &mate_ov, int32_t min_overlap_length, const int32_t error_threshold_for_merging) {
  std::string st_ov = "";
  st_ov = getOverlap(seq1, seq2, true, min_overlap_length, error_threshold_for_merging); // try with dovetail
  // std::string type_ov = getOverlap(seq1, seq2, false, min_overlap_length);

  if(st_ov != "") {
    mate_ov.ov_type = MateTypeOverlap::doveTail;
  } else {
    st_ov = getOverlap(seq1, seq2, false, min_overlap_length, error_threshold_for_merging); // try without dovetail
    mate_ov.ov_type = st_ov != "" ? MateTypeOverlap::overlap : MateTypeOverlap::noOverlap;
  }
  mate_ov.frag_length = st_ov.length();
  mate_ov.frag = st_ov;
}

// void write_mate_match(fastx_parser::FastxParser<fastx_parser::ReadPair>& parser,
//             std::atomic<uint64_t>& global_nr, std::atomic<uint64_t>& global_nhits,
//             std::mutex& iomut, ofstream& myfile) {
//   auto rg = parser.getReadGroup();
//   // uint32_t count = 0;
  
//   uint32_t max_chunk_reads = 100000;
//   vector<int32_t> frag_length(max_chunk_reads, -1);
//   vector<int32_t> map_type(max_chunk_reads, -1);
//   vector<int32_t> read1_length(max_chunk_reads, 0);
//   vector<int32_t> read2_length(max_chunk_reads, 0);
//   std::map<MateTypeOverlap, int> map;
//   map[doveTail] = 0;
//   map[overlap] = 1;
//   map[noOverlap] = 2;
//   map[both] = 3;
//   uint64_t read_num = 0;
//     // reserve the space to later write
//     // down the number of reads in the
//     // first chunk.
//   uint32_t num_reads_in_chunk{0};
  
//   // std::cout << "yp11" << std::endl;
//   while (parser.refill(rg)) {
//     for (auto& record : rg) {
//       ++global_nr;
//       ++read_num;
//       auto rctr = global_nr.load();
//       auto hctr = global_nhits.load();
//       MateOverlap mov;
      
//       findOverlapBetweenPairedEndReads(record.first.seq, record.second.seq, mov, 30);
//       frag_length[num_reads_in_chunk] = mov.frag_length;
//       map_type[num_reads_in_chunk] = map[mov.ov_type];
//       read1_length[num_reads_in_chunk] = static_cast<int32_t>(record.first.seq.length());
//       read2_length[num_reads_in_chunk] = static_cast<int32_t>(record.second.seq.length());
//       ++num_reads_in_chunk;
      
//       if (num_reads_in_chunk > max_chunk_reads) {
//         iomut.lock();
//         for(int i=0; i < frag_length.size(); i++){
//           myfile << frag_length[i] << "\t" << map_type[i] << "\t" << read1_length[i] << "\t" << read2_length[i]<< std::endl;
//         }
//         iomut.unlock();
//         frag_length.clear();
//         map_type.clear();
//         num_reads_in_chunk = 0;
//                 // std::cout<<"so"<<std::endl;
//                 // reserve space for headers of next chunk
//       }
    
//     }
//   }
//   if (num_reads_in_chunk > 0) { 
//     iomut.lock();
//     for(uint32_t i=0; i < num_reads_in_chunk; i++){
//       myfile << frag_length[i] << "\t" << map_type[i] << "\t" << read1_length[i] << "\t" << read2_length[i] << std::endl;
//     }
//     iomut.unlock();
//     frag_length.clear();
//     map_type.clear();
//     num_reads_in_chunk = 0;
//   }
// }
}
// int main() {
//   // std::string seq1 = "CCTTCTCCTAATTCCGCAAATGTGAAGGGTAGGGGGACGTTAAGGACCTG";
//   // std::string seq2 = "GTCCTTAACGTCCCCCTACCCTTCACACAGAGTGCTAAAGGAGAAGGCT";
//   // MateOverlap mov {-1, both};
//   // findOverlapBetweenPairedEndReads(seq1, seq2, mov);
//   std::ios_base::sync_with_stdio(false);

//   // pesc_options po;
//   // CLI::App app{"PESC — single-cell ATAC-seq  mapper"};
//   // app.add_option("-i,--index", po.index_basename, "input index prefix")->required();
//   // app.add_option("-1,--read1", po.left_read_filenames, "path to list of read 1 files")
//   //     ->required()
//   //     ->delimiter(',');
//   // app.add_option("-2,--read2", po.right_read_filenames, "path to list of read 2 files")
//   //     ->required()
//   //     ->delimiter(',');
//   // app.add_option("-b,--barcode_file", po.right_read_filenames, "path to list of barcode files")
//   //     ->required()
//   //     ->delimiter(',');
//   // app.add_option("-o,--output", po.output_dirname, "path to output directory")->required();
//   // // app.add_option("-g,--geometry", po.library_geometry, "geometry of barcode, umi and read")
//   // //     ->required();
//   // app.add_option("-t,--threads", po.nthread,
//   //                 "An integer that specifies the number of threads to use")
//   //     ->default_val(16);
//   // app.add_flag("--quiet", po.quiet, "try to be quiet in terms of console output");
//   // auto check_ambig =
//   //     app.add_flag("--check-ambig-hits", po.check_ambig_hits,
//   //                   "check the existence of highly-frequent hits in mapped targets, rather than "
//   //                   "ignoring them.");
//   // app.add_option("--max-ec-card", po.max_ec_card,
//   //                 "determines the maximum cardinality equivalence class "
//   //                 "(number of (txp, orientation status) pairs) to examine "
//   //                 "if performing check-ambig-hits.")
//   //     ->needs(check_ambig)
//   //     ->default_val(256);
//   // CLI11_PARSE(app, argc, argv);

//   ofstream myfile("../mate_info.txt");
  
//   spdlog::drop_all();
//   auto logger = spdlog::create<spdlog::sinks::stderr_color_sink_mt>("");
//   logger->set_pattern("%+");
  
//   spdlog::set_default_logger(logger);
//   spdlog::info("starting.");
//   // std::vector<std::string> read_file1 = {"/fs/cbcb-lab/rob/students/noor/Atacseq/atac_v1_pbmc_10k_fastqs/atac_v1_pbmc_10k_S1_L001_R1_001.fastq.gz"};
//   // std::vector<std::string> read_file2 = {"/fs/cbcb-lab/rob/students/noor/Atacseq/atac_v1_pbmc_10k_fastqs/atac_v1_pbmc_10k_S1_L001_R3_001.fastq.gz"};
//   std::vector<std::string> read_file1 = {"/fs/cbcb-lab/rob/students/noor/Atacseq/l1.fastq"};
//   std::vector<std::string> read_file2 = {"/fs/cbcb-lab/rob/students/noor/Atacseq/l2.fastq"};
//   std::vector<std::string> read_file3 = {"/fs/cbcb-lab/rob/students/noor/Atacseq/l3.fastq"};
//   std::mutex iomut;
//   std::atomic<uint64_t> global_nr{0};
//   std::atomic<uint64_t> global_nh{0};
//   // std::cout << read_file1 << "\n";

//   auto start_t = std::chrono::high_resolution_clock::now();

//   uint32_t np = 1;
//   size_t nthread{16};
//   fastx_parser::FastxParser<fastx_parser::ReadPair> rparser(
//           read_file1, read_file2, nthread, np);
//   rparser.start();
//   if (nthread >= 6) {
//         np += 1;
//         nthread -= 1;
//   }
//   std::vector<std::thread> workers;
//   for (size_t i = 0; i < nthread; ++i) {
//     workers.push_back(std::thread(
//       [&rparser, &global_nr, &global_nh, &iomut, &myfile]() {
//           write_mate_match(rparser, global_nr, global_nh, iomut, myfile);
//       }));
//   }
//   for (auto& w : workers) { w.join(); }
//     rparser.stop();
//   myfile.close();
//   spdlog::info("finished matching.");
  
//   auto end_t = std::chrono::high_resolution_clock::now();
//   auto num_sec = std::chrono::duration_cast<std::chrono::seconds>(end_t - start_t);
//   // spdlog::info(std::to_string(num_sec));
//     // auto rg = rparser.getReadGroup();

//     // std::atomic<uint64_t> global_nr=0;
//     // uint16_t max_l = 0;
//     // while (rparser.refill(rg)) {
//     //     for (auto& record : rg) {
//     //         ++global_nr;
//     //         if(static_cast<uint16_t>(record.second.seq.length()) > max_l)
//     //             max_l = static_cast<uint16_t>(record.second.seq.length());

//     //         if(global_nr >= 10000)
//     //             break;
//     //     }
//     //     if(global_nr >= 10)
//     //         break;
//     // }
//     // std::cout << "max is" << max_l << "\n";
//     // rparser.stop();
//     // std::cout << "ass\n";
//     return 0;

// }