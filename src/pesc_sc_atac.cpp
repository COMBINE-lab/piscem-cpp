#include "../external/sshash/external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "../external/sshash/include/util.hpp"
#include "../include/CanonicalKmerIterator.hpp"
#include "../include/FastxParser.hpp"
#include "../include/Kmer.hpp"
#include "../include/cli11/CLI11.hpp"
#include "../include/ghc/filesystem.hpp"
#include "../include/mapping/utils.hpp"
#include "../include/mapping/utils_bin.hpp"
#include "../include/meta_info.hpp"
#include "../include/parallel_hashmap/phmap.h"
#include "../include/projected_hits.hpp"
#include "../include/rad/rad_header.hpp"
#include "../include/rad/rad_writer.hpp"
#include "../include/rad/util.hpp"
#include "../include/reference_index.hpp"
#include "../include/sc/util.hpp"
#include "../include/spdlog_piscem/sinks/stdout_color_sinks.h"
#include "../include/spdlog_piscem/spdlog.h"
#include "../include/streaming_query.hpp"
#include "../include/util_piscem.hpp"
#include "../include/boost/unordered/concurrent_flat_map.hpp"
#include "check_overlap.cpp"
// #include "FastxParser.cpp"
// #include "hit_searcher.cpp"
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
#include <vector>

using namespace klibpp;
using BarCodeRecovered = single_cell::util::BarCodeRecovered;
using bc_kmer_t = rad::util::bc_kmer_t;
using mapping::util::bin_pos;
using mapping::util::mapping_cache_info;
using mapping::util::poison_state_t;

struct pesc_atac_options {
  std::string index_basename;
  std::vector<std::string> single_read_filenames;
  std::vector<std::string> left_read_filenames;
  std::vector<std::string> right_read_filenames;
  std::vector<std::string> barcode_filenames;
  std::string output_dirname;
  bool no_poison{true};
  bool use_sam_format{false};
  bool use_bed_format{false};
  bool check_kmers_orphans{false};
  bool use_chr{false};
  bool tn5_shift{true};
  bool enable_structural_constraints{false};
  float thr{0.7};
  bool quiet{false};
  bool check_ambig_hits{false};
  uint16_t blen{16};
  uint32_t max_ec_card{256};
  uint64_t bin_size{1000};
  uint64_t bin_overlap{300};
  size_t nthread{16};
};

template <typename mapping_cache_info_t>
bool map_fragment(
  fastx_parser::ReadTrip &record, poison_state_t &poison_state,
  mapping_cache_info_t &map_cache_left, mapping_cache_info_t &map_cache_right,
  mapping_cache_info_t &map_cache_out, std::atomic<uint64_t> &k_match,
  std::atomic<uint64_t> &l_match, std::atomic<uint64_t> &r_match,
  std::atomic<uint64_t> &dove_num, std::atomic<uint64_t> &dove_match,
  std::atomic<uint64_t> &ov_num, std::atomic<uint64_t> &ov_match,
  std::atomic<uint64_t> &r_orphan, std::atomic<uint64_t> &l_orphan,
  bool check_kmers_orphans, mapping::util::bin_pos &binning, bool use_chr) {

  bool km = false; // kmatch checker
  map_cache_out.clear();
  poison_state.clear();
  poison_state.set_fragment_end(mapping::util::fragment_end::LEFT);

  check_overlap::MateOverlap mate_ov;
  check_overlap::findOverlapBetweenPairedEndReads(
    record.first.seq, record.second.seq, mate_ov, 30, 0);
  if (mate_ov.frag != "") {
    bool exit = mapping::util::map_read(&mate_ov.frag, map_cache_out,
                                        poison_state, binning, km, use_chr);
    if (km) {
      ++k_match;
    }

    // std::cout << map_cache_out.accepted_hits.size() << std::endl;
    if (!map_cache_out.accepted_hits.empty()) {
      uint32_t max_num_hits = map_cache_out.accepted_hits.front().num_hits;
      std::sort(map_cache_out.accepted_hits.begin(),
                map_cache_out.accepted_hits.end(),
                mapping::util_bin::simple_hit_less_bins);
      mapping::util_bin::remove_duplicate_hits(map_cache_out, max_num_hits);
      map_cache_out.map_type = (map_cache_out.accepted_hits.size() > 0)
                                 ? mapping::util::MappingType::MAPPED_PAIR
                                 : mapping::util::MappingType::UNMAPPED;
      for (auto &hit : map_cache_out.accepted_hits) {
        hit.fragment_length = mate_ov.frag_length;
        int32_t r2_len = record.first.seq.length() <= record.second.seq.length()
                           ? record.second.seq.length()
                           : record.first.seq.length();
        int32_t r1_len = record.first.seq.length() <= record.second.seq.length()
                           ? record.first.seq.length()
                           : record.second.seq.length();
        const int32_t ref_len =
          static_cast<int32_t>(map_cache_out.hs.get_index()->ref_len(hit.tid));
        hit.mate_pos = hit.is_fw ? hit.pos + hit.fragment_length - r2_len - 1
                                 : hit.pos + r1_len - hit.fragment_length - 1;
        if (hit.mate_pos > ref_len) {
          hit.mate_pos = hit.pos;
        }
        if (mate_ov.ov_type == check_overlap::MateTypeOverlap::doveTail) {
          hit.mate_pos = hit.pos;
        }
      }
      map_cache_out.frag_seq = mate_ov.frag;
      map_cache_out.read1 = mate_ov.frag_fw;
    }
    // add remove max_hits
    if (mate_ov.ov_type == check_overlap::MateTypeOverlap::doveTail) {
      dove_match += map_cache_out.accepted_hits.empty() ? 0 : 1;
      dove_num += 1;
    } else {
      ov_match += map_cache_out.accepted_hits.empty() ? 0 : 1;
      ov_num += 1;
    }
    return exit;
  }

  bool early_exit_left = mapping::util::map_read(
    &record.first.seq, map_cache_left, poison_state, binning, km, use_chr);
  // bool l_map=false,r_map=false;
  // if (map_cache_left.accepted_hits.size() > 0 &&
  // map_cache_left.accepted_hits.size() < 5) {
  //     l_map=true;
  // }
  // std::cout << "record is " << record.first.name << std::endl;
  // std::cout << " left\n";
  // mapping::util::print_hits(map_cache_left.accepted_hits);

  if (poison_state.is_poisoned()) {
    return false;
  }
  bool right_km = false;
  poison_state.set_fragment_end(mapping::util::fragment_end::RIGHT);
  bool early_exit_right =
    mapping::util::map_read(&record.second.seq, map_cache_right, poison_state,
                            binning, right_km, use_chr);
  // if (map_cache_right.accepted_hits.size() > 0 &&
  // map_cache_right.accepted_hits.size() < 5) {
  //     r_map=true;
  // }
  // std::cout << " right\n";
  // mapping::util::print_hits(map_cache_right.accepted_hits);

  // if(l_map && r_map) {
  //     std::cout << " left\n";
  //     mapping::util::print_hits(map_cache_left.accepted_hits);
  //     std::cout << " right\n";
  //     mapping::util::print_hits(map_cache_right.accepted_hits);
  // }
  // mapping::util::print_hits(map_cache_right.accepted_hits);
  if (poison_state.is_poisoned()) {
    return false;
  }
  if (km or right_km) {
    ++k_match;
  }

  int32_t left_len = static_cast<int32_t>(record.first.seq.length());
  int32_t right_len = static_cast<int32_t>(record.second.seq.length());

  l_match += map_cache_left.accepted_hits.empty() ? 0 : 1;
  r_match += map_cache_right.accepted_hits.empty() ? 0 : 1;
  mapping::util_bin::merge_se_mappings(map_cache_left, map_cache_right,
                                       left_len, right_len, check_kmers_orphans,
                                       map_cache_out);
  // if (l_map && r_map && map_cache_out.accepted_hits.empty()) {
  //     std::cout << record.first.name << std::endl;
  //     std::cout << "merge not mapping\n";
  // }

  // std::cout << "merged\n";
  // mapping::util::print_hits(map_cache_out.accepted_hits);
  l_orphan +=
    map_cache_out.map_type == mapping::util::MappingType::MAPPED_FIRST_ORPHAN
      ? 1
      : 0;
  r_orphan +=
    map_cache_out.map_type == mapping::util::MappingType::MAPPED_SECOND_ORPHAN
      ? 1
      : 0;
  if (map_cache_out.map_type ==
        mapping::util::MappingType::MAPPED_FIRST_ORPHAN ||
      map_cache_out.map_type ==
        mapping::util::MappingType::MAPPED_SECOND_ORPHAN) {
    for (auto &hit : map_cache_out.accepted_hits) {
      hit.fragment_length = map_cache_out.map_type ==
                                mapping::util::MappingType::MAPPED_FIRST_ORPHAN
                              ? record.first.seq.length()
                              : record.second.seq.length();
    }
  }
  return (early_exit_left or early_exit_right);
}

template <typename mapping_cache_info_t>
bool map_fragment(
  fastx_parser::ReadPair &record, poison_state_t &poison_state,
  mapping_cache_info_t &map_cache_left, mapping_cache_info_t &map_cache_right,
  mapping_cache_info_t &map_cache_out, std::atomic<uint64_t> &k_match,
  std::atomic<uint64_t> &l_match, std::atomic<uint64_t> &r_match,
  std::atomic<uint64_t> &dove_num, std::atomic<uint64_t> &dove_match,
  std::atomic<uint64_t> &ov_num, std::atomic<uint64_t> &ov_match,
  std::atomic<uint64_t> &r_orphan, std::atomic<uint64_t> &l_orphan,
  bool check_kmers_orphans, mapping::util::bin_pos &binning, bool use_chr) {

  (void)map_cache_left;
  (void)map_cache_right;
  (void)l_match;
  (void)r_match;
  (void)dove_match;
  (void)dove_num;
  (void)ov_num;
  (void)ov_match;
  (void)l_orphan;
  (void)r_orphan;
  (void)check_kmers_orphans;
  bool km = false; // kmatch checker
  map_cache_out.clear();
  poison_state.clear();
  poison_state.set_fragment_end(mapping::util::fragment_end::LEFT);

  bool early_exit_left = mapping::util::map_read(
    &record.first.seq, map_cache_out, poison_state, binning, km, use_chr);
  if (poison_state.is_poisoned()) {
    return false;
  }
  if (km) {
    ++k_match;
  }

  int32_t left_len = static_cast<int32_t>(record.first.seq.length());

  if (!map_cache_out.accepted_hits.empty()) {
    uint32_t max_num_hits = map_cache_out.accepted_hits.front().num_hits;
    std::sort(map_cache_out.accepted_hits.begin(),
              map_cache_out.accepted_hits.end(),
              mapping::util_bin::simple_hit_less_bins);
    mapping::util_bin::remove_duplicate_hits(map_cache_out, max_num_hits);
    map_cache_out.map_type = (map_cache_out.accepted_hits.size() > 0)
                               ? mapping::util::MappingType::SINGLE_MAPPED
                               : mapping::util::MappingType::UNMAPPED;
    for (auto &hit : map_cache_out.accepted_hits) {
      hit.fragment_length = left_len;
    }
    map_cache_out.frag_seq = record.first.seq;
  }
  return early_exit_left;
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

  std::ofstream sam_file;
  std::mutex sam_mutex;
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

  // the file where the sequence from first fastq and string mapping will be
  // stored
};

template <typename mapping_cache_info_t>
inline void
write_sam_mappings(mapping_cache_info_t &map_cache_out, bc_kmer_t &bck,
                   phmap::flat_hash_map<uint64_t, uint32_t> &unmapped_bc_map,
                   fastx_parser::ReadPair &record, std::string &workstr_left,
                   std::atomic<uint64_t> &global_nhits,
                   std::ostringstream &osstream) {
  (void)workstr_left;
  constexpr uint16_t is_secondary = 256;
  constexpr uint16_t is_rc = 16;

  if (!map_cache_out.accepted_hits.empty()) {
    ++global_nhits;
    if (map_cache_out.map_type == mapping::util::MappingType::UNMAPPED) {
      unmapped_bc_map[bck.word(0)] += 1;
    }
    bool secondary = false;
    for (auto &ah : map_cache_out.accepted_hits) {
      uint16_t flag = secondary ? is_secondary : 0;
      // flag += 2;
      flag += ah.is_fw ? 0 : is_rc;
      // flag += first_seg;

      std::string *sptr = nullptr;
      if (is_rc) {
        combinelib::kmers::reverseComplement(record.first.seq, workstr_left);
        sptr = &workstr_left;
      } else {
        sptr = &record.first.seq;
      }
      osstream << record.first.name << "\t" << flag << "\t"
               << map_cache_out.hs.get_index()->ref_name(ah.tid) << "\t"
               << ah.pos + 1 << "\t255\t*\t*\t0\t" << record.first.seq.length()
               << "\t" << *sptr << "\t*\n";
      secondary = true;
    }
  } else {
    osstream << record.first.name << "\t" << 4 << "\t"
             << "*\t0\t0\t*\t*\t0\t0\t" << record.first.seq << "\t*\n";
  }
}

template <typename mapping_cache_info_t>
inline void
write_sam_mappings(mapping_cache_info_t &map_cache_out, bc_kmer_t &bck,
                   phmap::flat_hash_map<uint64_t, uint32_t> &unmapped_bc_map,
                   fastx_parser::ReadTrip &record, std::string &workstr_left,
                   std::string &workstr_right,
                   std::atomic<uint64_t> &global_nhits,
                   std::ostringstream &osstream) {

  (void)workstr_right;
  constexpr uint16_t is_secondary = 256;
  constexpr uint16_t is_rc = 16;
  constexpr uint16_t mate_rc = 32;
  constexpr uint16_t unmapped = 4;
  constexpr uint16_t mate_unmapped = 8;

  auto map_type = map_cache_out.map_type;

  bool mated_before_mapping = !map_cache_out.frag_seq.empty();

  if (map_type == mapping::util::MappingType::UNMAPPED) {
    unmapped_bc_map[bck.word(0)] += 1;
  }

  if (map_type == mapping::util::MappingType::SINGLE_MAPPED) {
    if (!map_cache_out.accepted_hits.empty()) {
      ++global_nhits;
      bool secondary = false;
      for (auto &ah : map_cache_out.accepted_hits) {
        uint16_t flag = secondary ? is_secondary : 0;
        // flag += 2;
        flag += ah.is_fw ? 0 : is_rc;
        // flag += first_seg;

        std::string *sptr = nullptr;
        if (is_rc) {
          combinelib::kmers::reverseComplement(record.first.seq, workstr_left);
          sptr = &workstr_left;
        } else {
          sptr = &record.first.seq;
        }
        osstream << record.first.name << "\t" << flag << "\t"
                 << map_cache_out.hs.get_index()->ref_name(ah.tid) << "\t"
                 << ah.pos + 1 << "\t255\t*\t*\t0\t"
                 << record.first.seq.length() << "\t" << *sptr << "\t*\n";
        secondary = true;
      }
    } else {
      osstream << record.first.name << "\t" << 4 << "\t"
               << "*\t0\t0\t*\t*\t0\t0\t" << record.first.seq << "\t*\n";
    }
    return;
  }
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
    for (auto &ah : map_cache_out.accepted_hits) {
      uint16_t flag_first =
        secondary ? base_flag_first + is_secondary : base_flag_first;
      uint16_t flag_second =
        secondary ? base_flag_second + is_secondary : base_flag_second;

      std::string *sptr_first = nullptr;
      std::string *sptr_second = nullptr;
      int32_t pos_first = 0;
      int32_t pos_second = 0;

      // if both reads are mapped
      if (map_type == mapping::util::MappingType::MAPPED_PAIR) {
        pos_first = ah.pos + 1;
        pos_second = ah.mate_pos + 1;

        if (ah.is_fw) {
          flag_first += mate_rc;
          sptr_first =
            mated_before_mapping ? &map_cache_out.frag_seq : &record.first.seq;

          flag_second += is_rc;
          if (!have_rc_second) {
            have_rc_second = true;
            if (mated_before_mapping) {
              combinelib::kmers::reverseComplement(map_cache_out.frag_seq,
                                                   workstr_right);
            } else {
              combinelib::kmers::reverseComplement(record.second.seq,
                                                   workstr_right);
            }
          }
          sptr_second = &workstr_right;

        } else {
          flag_first += is_rc;
          if (!have_rc_first) {
            have_rc_first = true;
            if (mated_before_mapping) {
              combinelib::kmers::reverseComplement(map_cache_out.frag_seq,
                                                   workstr_left);
            } else {
              combinelib::kmers::reverseComplement(record.first.seq,
                                                   workstr_left);
            }
          }
          sptr_first = &workstr_left;

          flag_second += mate_rc;
          sptr_second =
            mated_before_mapping ? &map_cache_out.frag_seq : &record.second.seq;
        }
      } else if (map_type == mapping::util::MappingType::MAPPED_FIRST_ORPHAN) {
        pos_first = ah.pos + 1;
        pos_second = 0;

        sptr_first = &record.first.seq;
        sptr_second = &record.second.seq;

        if (!ah.is_fw) { // if the mapped read is rc
          flag_first += is_rc;
          if (!have_rc_first) {
            have_rc_first = true;
            combinelib::kmers::reverseComplement(record.first.seq,
                                                 workstr_left);
          }
          sptr_first = &workstr_left;

          flag_second += mate_rc;
        }

      } else if (map_type == mapping::util::MappingType::MAPPED_SECOND_ORPHAN) {
        pos_first = 0;
        pos_second = ah.pos + 1;

        sptr_first = &record.first.seq;
        sptr_second = &record.second.seq;
        if (!ah.is_fw) {
          flag_first += mate_rc;
          flag_second += is_rc;
          if (!have_rc_second) {
            have_rc_second = true;
            combinelib::kmers::reverseComplement(record.second.seq,
                                                 workstr_right);
          }
          sptr_second = &workstr_right;
        }
      }

      auto print_pos_mapq_cigar = [](bool mapped, int32_t pos, int32_t read_len,
                                     int32_t ref_len, std::ostream &os) {
        if (!mapped) {
          os << "0\t255\t*\t";
          return;
        }
        int32_t pad_start = 0;
        int32_t pad_end = 0;
        int32_t m_len = read_len;

        if (pos + read_len >= ref_len) {
          pad_end = (pos + read_len) - ref_len + 1;
          m_len -= pad_end;
        }

        if (pos <= 0) {
          pad_start = (-pos) + 1;
          m_len -= pad_start;
          pos = 1;
        }

        os << pos << "\t255\t";

        if (pad_start > 0) {
          os << pad_start << "S";
        }
        if (pad_end > 0) {
          os << m_len << 'M' << pad_end << "S\t";
        } else {
          os << m_len << "M\t";
        }
      };

      const auto ref_name = map_cache_out.hs.get_index()->ref_name(ah.tid);
      const int32_t ref_len =
        static_cast<int32_t>(map_cache_out.hs.get_index()->ref_len(ah.tid));
      std::string r1name = record.first.name;
      std::string r2name = record.second.name;

      if (mated_before_mapping && !map_cache_out.read1) {
        r1name = record.second.name;
        r2name = record.first.name;
      }
      //   std::string r1name = (mated_before_mapping && !map_cache_out.read1) ?
      //   record.second.name :
      int32_t r1len = mated_before_mapping ? map_cache_out.frag_seq.length()
                                           : record.first.seq.length();
      //   std::string r2name = (mated_before_mapping && map_cache_out.read1) ?
      //   record.second.name : record.first.name;
      int32_t r2len = mated_before_mapping ? map_cache_out.frag_seq.length()
                                           : record.second.seq.length();
      // if (tn5_shift) {
      //     if (pos_first <= pos_second) {
      //         pos_first += 4;
      //     } else {
      //         pos_second += 4;
      //     }
      // }
      // auto frag_len = tn5_shift ? ah.frag_len()-9 : ah.frag_len();
      auto frag_len = ah.frag_len();
      osstream << r1name << "\t" << flag_first << "\t"
               << ((flag_first & unmapped) ? "*" : ref_name)
               << '\t'; // if mapped RNAME, else *

      print_pos_mapq_cigar(!(flag_first & unmapped), pos_first,
                           static_cast<int32_t>(r1len), ref_len, osstream);

      osstream << ((flag_first & mate_unmapped) ? '*' : '=') << '\t' // RNEXT
               << ((flag_first & mate_unmapped) ? 0 : std::max(1, pos_second))
               << '\t' // PNEXT
               << frag_len * (ah.pos < ah.mate_pos ? 1 : -1) << '\t'
               << *sptr_first << "\t*\n";
      osstream << r2name << "\t" << flag_second << "\t"
               << ((flag_second & unmapped) ? "*" : ref_name)
               << '\t'; // if mapped RNAME, else *
      print_pos_mapq_cigar(!(flag_second & unmapped), pos_second,
                           static_cast<int32_t>(r2len), ref_len, osstream);
      osstream << ((flag_second & mate_unmapped) ? '*' : '=') << '\t' // RNEXT
               << ((flag_second & mate_unmapped) ? 0 : std::max(1, pos_first))
               << '\t' // PNEXT
               << -frag_len * (ah.pos < ah.mate_pos ? 1 : -1) << '\t'
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
struct RadT {};
struct SamT {};
template <typename FragT, typename SketchHitT, typename OutputT = RadT>
void do_map(mindex::reference_index &ri,
            fastx_parser::FastxParser<FragT> &parser,
            mapping::util::bin_pos &binning, poison_table &poison_map,
            std::atomic<uint64_t> &global_nr,
            std::atomic<uint64_t> &global_nhits,
            std::atomic<uint64_t> &global_nmult, std::atomic<uint64_t> &k_match,
            std::atomic<uint64_t> &l_match, std::atomic<uint64_t> &r_match,
            std::atomic<uint64_t> &dove_num, std::atomic<uint64_t> &dove_match,
            std::atomic<uint64_t> &ov_num, std::atomic<uint64_t> &ov_match,
            std::atomic<uint64_t> &r_orphan, std::atomic<uint64_t> &l_orphan,
            std::atomic<uint64_t> &global_npoisoned, pesc_output_info &out_info,
            std::mutex &iomut, bool write_bed, bool check_kmers_orphans,
            bool tn5_shift, bool use_chr, 
            boost::concurrent_flat_map<uint64_t, sshash::lookup_result>& unitig_end_cache,
            RAD::RAD_Writer &rw,
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

  constexpr bool paired_end_frags =
    std::is_same_v<fastx_parser::ReadTrip, FragT>;
  // the reads are paired
  if constexpr (paired_end_frags) {
    poison_state.paired_for_mapping = true;
  }
  // SAM output
  uint64_t processed = 0;
  uint64_t buff_size = 10000;
  // these don't really belong here
  std::string workstr_left;
  std::string workstr_right;
  std::ostringstream osstream;

  mapping_cache_info<SketchHitT, piscem::streaming_query<true>> map_cache_left(ri, &unitig_end_cache);
  mapping_cache_info<SketchHitT, piscem::streaming_query<true>> map_cache_right(ri, &unitig_end_cache);
  mapping_cache_info<SketchHitT, piscem::streaming_query<true>> map_cache_out(ri, &unitig_end_cache);

  size_t max_chunk_reads = 5000;

  auto rg = parser.getReadGroup();
  uint32_t num_reads_in_chunk{0};

  // for(int32_t i = 0; i < ri.num_refs(); i++) {
  //     std::cout << "i is " <<  i << " ref len is " << ri.ref_len(i) << " name
  //     is " << ri.ref_name(i) << std::endl;
  // }

  mindex::hit_searcher hs(&ri);
  uint64_t read_num = 0;
  (void)read_num;

  std::string temp_buff = "";
  while (parser.refill(rg)) {
    for (auto &record : rg) {
      ++global_nr;
      ++read_num;
      auto rctr = global_nr.load();
      auto hctr = global_nhits.load();

      if (write_mapping_rate and (rctr % 500000 == 0)) {
        iomut.lock();
        std::cerr << "\rprocessed (" << rctr << ") reads; (" << hctr
                  << ") had mappings.";
        iomut.unlock();
      }

      std::string *bc{nullptr};
      if constexpr (paired_end_frags) {
        bc = &record.third.seq;
      } else {
        bc = &record.second.seq;
      }
      bc_kmer_t bc_kmer;

      auto recovered = single_cell::util::recover_barcode(*bc);
      if (recovered == BarCodeRecovered::NOT_RECOVERED) {
        continue;
      }

      bool bc_ok = bc_kmer.fromChars(*bc);
      if (!bc_ok) {
        continue;
      }

      // if we couldn't correct it with 1 `N`, then skip.
      bool had_early_stop = map_fragment(
        record, poison_state, map_cache_left, map_cache_right, map_cache_out,
        k_match, l_match, r_match, dove_num, dove_match, ov_num, ov_match,
        r_orphan, l_orphan, check_kmers_orphans, binning, use_chr);
      (void)had_early_stop;
      // if (had_early_stop) {

      // }
      if (poison_state.is_poisoned()) {
        global_npoisoned++;
      }

      global_nmult += map_cache_out.accepted_hits.size() > 1 ? 1 : 0;

      if constexpr (std::is_same_v<OutputT, RadT>) {
        global_nhits += map_cache_out.accepted_hits.empty() ? 0 : 1;
        rad::util::write_to_rad_stream_atac(
          bc_kmer, map_cache_out.map_type, map_cache_out.accepted_hits,
          map_cache_out.unmapped_bc_map, num_reads_in_chunk, temp_buff, *bc, ri,
          rw, token, tn5_shift);
      } else {
        // silence warning
        (void)token;
        (void)tn5_shift;
      }

      if constexpr (std::is_same_v<OutputT, SamT>) {
        ++processed;
        // mapping::util::print_hits(map_cache_out.accepted_hits);
        if constexpr (std::is_same_v<fastx_parser::ReadTrip, FragT>) {
          write_sam_mappings(
            map_cache_out, bc_kmer, map_cache_out.unmapped_bc_map, record,
            workstr_left, workstr_right, global_nhits, osstream);
        } else {
          write_sam_mappings(map_cache_out, bc_kmer,
                             map_cache_out.unmapped_bc_map, record,
                             workstr_left, global_nhits, osstream);
        }

        if (processed >= buff_size) {

          std::string o = osstream.str();
          out_info.sam_mutex.lock();
          out_info.sam_file << o;
          out_info.sam_mutex.unlock();
          osstream.clear();
          osstream.str("");
          processed = 0;
        }
      }
      // dump buffer
      if (num_reads_in_chunk > max_chunk_reads) {
        if (write_bed) {
          out_info.bed_mutex.lock();
          out_info.bed_file << temp_buff;
          out_info.bed_mutex.unlock();
        }
        out_info.num_chunks++;
        temp_buff = "";
        num_reads_in_chunk = 0;
      }
    }
  }

  if (num_reads_in_chunk > 0) {
    if (write_bed) {
      out_info.bed_mutex.lock();
      out_info.bed_file << temp_buff;
      out_info.bed_mutex.unlock();
    }
    out_info.num_chunks++;
    temp_buff = "";
    num_reads_in_chunk = 0;
  }

  if (processed > 0) {
    std::string o = osstream.str();
    out_info.sam_mutex.lock();
    out_info.sam_file << o;
    out_info.sam_mutex.unlock();
    osstream.clear();
    osstream.str("");
    processed = 0;
  }

  // // unmapped barcode writer
  { // make a scope and dump the unmapped barcode counts
    std::string ubcw = "";
    for (auto &kv : map_cache_out.unmapped_bc_map) {
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
  // std::cerr << "left non matches : " << map_cache_left.hs.non_matches_cnt <<
  // "\n";
}

#ifdef __cplusplus
extern "C" {
#endif
int run_pesc_sc_atac(int argc, char **argv);
#ifdef __cplusplus
}
#endif

void print_header(mindex::reference_index &ri, std::ostringstream &osstream) {
  //   cmdline << "@HD\tVN:1.0\tSO:unsorted\n";
  for (uint64_t i = 0; i < ri.num_refs(); ++i) {
    osstream << "@SQ\tSN:" << ri.ref_name(i) << "\tLN:" << ri.ref_len(i)
             << "\n";
  }
  //   cmdline << "@PG\tID:mindex_map\tPN:mapper\tVN:0.0.1\t"
  //             << "CL:" << "\n";
}

int run_pesc_sc_atac(int argc, char **argv) {
  std::ios_base::sync_with_stdio(false);
  pesc_atac_options po;
  std::string skipping_rule;

  pesc_output_info out_info;
  size_t nthread{16};
  CLI::App app{"Single cell Atac Seq mapper"};
  app.add_option("-i,--index", po.index_basename, "input index prefix")
    ->required();
  auto ogroup = app.add_option_group("input reads", "provide input reads");

  CLI::Option *read_opt =
    ogroup
      ->add_option("-r,--reads", po.single_read_filenames,
                   "Path to list (comma separated) of single-end files")
      ->delimiter(',');

  CLI::Option *paired_left_opt =
    ogroup
      ->add_option("-1,--read1", po.left_read_filenames,
                   "Path to list (comma separated) of read 1 files")
      ->delimiter(',');
  CLI::Option *paired_right_opt =
    ogroup
      ->add_option("-2,--read2", po.right_read_filenames,
                   "Path to list (comma separated) of read 2 files")
      ->delimiter(',');

  paired_left_opt->excludes(read_opt);
  paired_right_opt->excludes(read_opt);
  read_opt->excludes(paired_left_opt, paired_right_opt);
  paired_left_opt->needs(paired_right_opt);
  paired_right_opt->needs(paired_left_opt);

  ogroup->require_option(1, 2);
  app
    .add_option("-b,--barcode", po.barcode_filenames,
                "path to list of barcodes")
    ->required()
    ->delimiter(',');
  app.add_option("-o,--output", po.output_dirname, "path to output directory")
    ->required();
  app
    .add_option("-t,--threads", po.nthread,
                "An integer that specifies the number of threads to use")
    ->default_val(16);
  app.add_flag("--sam-format", po.use_sam_format,
               "Write SAM format output rather than RAD.");
  app.add_flag("--kmers-orphans", po.check_kmers_orphans,
               "Check if any mapping kmer exist for a mate, if there exists "
               "mapping for the other read (default false)");
  app.add_flag("--bed-format", po.use_bed_format, "Dump output to bed.");
  app.add_flag("--use-chr", po.use_chr, "use chromosomes as virtual color.");
  app.add_option("--tn5-shift", po.tn5_shift, "Tn5 shift");
  app
    .add_option("--no-poison", po.no_poison,
                "Do not filter reads for poison k-mers, even if a poison table "
                "exists for the index")
    ->default_val(true);
  app.add_flag("-c,--struct-constraints", po.enable_structural_constraints,
               "Apply structural constraints when performing mapping");
  app
    .add_option("--skipping-strategy", skipping_rule,
                "Which skipping rule to use for pseudoalignment ({strict, "
                "permissive, strict})")
    ->default_val("permissive");
  app.add_flag("--quiet", po.quiet,
               "Try to be quiet in terms of console output");
  app.add_option("--thr", po.thr, "threshold for psa")->default_val(0.7);
  app.add_option("--bclen", po.blen, "length for barcode")->default_val(16);
  app.add_option("--bin-size", po.bin_size, "size for binning")
    ->default_val(1000);
  app.add_option("--bin-overlap", po.bin_overlap, "size for bin overlap")
    ->default_val(300);
  auto check_ambig = app.add_flag("--check-ambig-hits", po.check_ambig_hits,
                                  "check the existence of highly-frequent hits "
                                  "in mapped targets, rather than "
                                  "ignoring them.");
  (void)check_ambig; // currently unused in atacseq mode

  CLI11_PARSE(app, argc, argv);

  bool paired_end = !po.right_read_filenames.empty();

  spdlog_piscem::drop_all();
  auto logger =
    spdlog_piscem::create<spdlog_piscem::sinks::stdout_color_sink_mt>("");
  logger->set_pattern("%+");

  if (po.quiet) {
    spdlog_piscem::set_level(spdlog_piscem::level::warn);
  }
  spdlog_piscem::set_default_logger(logger);

  spdlog_piscem::info("enable structural constraints : {}",
                      po.enable_structural_constraints);

  std::optional<mindex::SkippingStrategy> skip_strat_opt =
    mindex::SkippingStrategy::from_string(skipping_rule);
  if (!skip_strat_opt) {
    spdlog_piscem::critical(
      "The skipping strategy must be one of \"strict\", "
      " \"every\" or \"permissive\", but \"{}\" was passed in",
      skipping_rule);
    return 1;
  }

  nthread = po.nthread;

  auto start_t = std::chrono::high_resolution_clock::now();

  ghc::filesystem::path output_path(po.output_dirname);
  ghc::filesystem::create_directories(output_path);

  ghc::filesystem::path bed_file_path = output_path / "map.bed";
  ghc::filesystem::path sam_file_path = output_path / "map.sam";
  ghc::filesystem::path unmapped_bc_file_path =
    output_path / "unmapped_bc_count.bin";
  ghc::filesystem::path mapped_bc_file_path = output_path / "mapped_bc.txt";

  if (po.use_bed_format) {
    std::ofstream bed_file(bed_file_path.string());
    if (!bed_file.good()) {
      spdlog_piscem::critical("Could not open {} for writing.",
                              bed_file_path.string());
      throw std::runtime_error("error creating bed file.");
    }
    out_info.bed_file = std::move(bed_file);
  }

  if (po.use_sam_format) {
    std::ofstream sam_file(sam_file_path.string());
    if (!sam_file.good()) {
      spdlog_piscem::critical("Could not open {} for writing.",
                              bed_file_path.string());
      throw std::runtime_error("error creating sam file.");
    }
    out_info.sam_file = std::move(sam_file);
  }

  std::ofstream unmapped_bc_file(unmapped_bc_file_path.string());
  if (!unmapped_bc_file.good()) {
    spdlog_piscem::critical("Could not open {} for writing.",
                            unmapped_bc_file_path.string());
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

  out_info.unmapped_bc_file = std::move(unmapped_bc_file);

  mindex::reference_index ri(po.index_basename, po.check_ambig_hits);
  std::string rad_file_path = output_path;
  rad_file_path.append("/map.rad");
  RAD::Tag_Defn tag_defn;
  RAD::Tag_List file_tag_vals;
  file_tag_vals.add(RAD::Type::u16(po.blen));
  std::vector<uint64_t> len;
  len.reserve(ri.num_refs());
  for (decltype(ri.num_refs()) i = 0; i < ri.num_refs(); i++) {
    len.push_back(ri.ref_len(i));
  }
  file_tag_vals.add(RAD::Type::v_u64(len));

  std::vector<std::string> refs;

  bc_kmer_t::k(16);
  rad::util::write_rad_header_atac(ri, refs, tag_defn);
  mapping::util::bin_pos binning(&ri, po.thr, po.bin_size, po.bin_overlap);

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
  std::ostringstream val;
  if (po.use_sam_format) {
    print_header(ri, val);
    out_info.sam_file << val.str();
  }

  uint32_t np = 1;

  std::atomic<uint64_t> global_nr{0};
  std::atomic<uint64_t> global_nh{0};
  std::atomic<uint64_t> global_nmult{0}; // number of multimapping
  std::atomic<uint64_t> k_match{
    0}; // whether the kmer exists in the unitig table
  std::atomic<uint64_t> l_match{0};
  std::atomic<uint64_t> r_match{0};
  std::atomic<uint64_t> dove_match{0};
  std::atomic<uint64_t> dove_num{0};
  std::atomic<uint64_t> ov_num{0};
  std::atomic<uint64_t> ov_match{0};
  std::atomic<uint64_t> l_orphan{0};
  std::atomic<uint64_t> r_orphan{0};
  std::atomic<uint64_t> global_np{
    0}; // whether the kmer exists in the unitig table
  std::mutex iomut;
  constexpr size_t unitig_end_cache_size{5000000};

  if (paired_end) {
    using FragmentT = fastx_parser::ReadTrip;
    
    auto num_input_files = po.left_read_filenames.size();
    size_t additional_files = (num_input_files > 1) ? (num_input_files - 1) : 0;

    // start with 1 parsing thread, and one more for every
    // 6 threads, as long as there are additional input files
    // to parse.
    size_t remaining_threads = nthread;
    for (size_t i = 0; i < additional_files; ++i) {
      if (remaining_threads >= 6) {
        np += 1;
        nthread -= 1;
        remaining_threads -= 6;
      } else {
        break;
      }
    }

    fastx_parser::FastxParser<fastx_parser::ReadTrip> rparser(
      po.left_read_filenames, po.right_read_filenames, po.barcode_filenames,
      nthread, np);

    rparser.start();
    boost::concurrent_flat_map<uint64_t, sshash::lookup_result> unitig_end_cache(unitig_end_cache_size);
    std::vector<std::thread> workers;
    for (size_t i = 0; i < nthread; ++i) {
      workers.push_back(std::thread([&ri, &po, &rparser, &binning, &ptab,
                                     &global_nr, &global_nh, &global_nmult,
                                     &k_match, &global_np, &out_info, &iomut,
                                     &rw, &l_match, &r_match, &dove_match,
                                     &dove_num, &ov_num, &ov_match, &r_orphan,
                                     &l_orphan, &unitig_end_cache]() {
        const auto token = rw.get_token();
        if (!po.enable_structural_constraints) {
          using SketchHitT =
            mapping::util::sketch_hit_info_no_struct_constraint;
          if (po.use_sam_format) {
            do_map<FragmentT, SketchHitT, SamT>(
              ri, rparser, binning, ptab, global_nr, global_nh, global_nmult,
              k_match, l_match, r_match, dove_num, dove_match, ov_num, ov_match,
              r_orphan, l_orphan, global_np, out_info, iomut, po.use_bed_format,
              po.check_kmers_orphans, po.tn5_shift, po.use_chr, unitig_end_cache, rw, token);
          } else {
            do_map<FragmentT, SketchHitT, RadT>(
              ri, rparser, binning, ptab, global_nr, global_nh, global_nmult,
              k_match, l_match, r_match, dove_num, dove_match, ov_num, ov_match,
              r_orphan, l_orphan, global_np, out_info, iomut, po.use_bed_format,
              po.check_kmers_orphans, po.tn5_shift, po.use_chr, unitig_end_cache, rw, token);
          }
        } else {
          using SketchHitT = mapping::util::sketch_hit_info;
          if (po.use_sam_format) {
            do_map<FragmentT, SketchHitT, SamT>(
              ri, rparser, binning, ptab, global_nr, global_nh, global_nmult,
              k_match, l_match, r_match, dove_num, dove_match, ov_num, ov_match,
              r_orphan, l_orphan, global_np, out_info, iomut, po.use_bed_format,
              po.check_kmers_orphans, po.tn5_shift, po.use_chr, unitig_end_cache, rw, token);
          } else {
            do_map<FragmentT, SketchHitT, RadT>(
              ri, rparser, binning, ptab, global_nr, global_nh, global_nmult,
              k_match, l_match, r_match, dove_num, dove_match, ov_num, ov_match,
              r_orphan, l_orphan, global_np, out_info, iomut, po.use_bed_format,
              po.check_kmers_orphans, po.tn5_shift, po.use_chr, unitig_end_cache, rw, token);
          }
        }
      }));
    }

    for (auto &w : workers) {
      w.join();
    }
    rparser.stop();
  } else {
    using FragmentT = fastx_parser::ReadPair;
 
    auto num_input_files = po.single_read_filenames.size();
    size_t additional_files = (num_input_files > 1) ? (num_input_files - 1) : 0;

    // start with 1 parsing thread, and one more for every
    // 6 threads, as long as there are additional input files
    // to parse.
    size_t remaining_threads = nthread;
    for (size_t i = 0; i < additional_files; ++i) {
      if (remaining_threads >= 6) {
        np += 1;
        nthread -= 1;
        remaining_threads -= 6;
      } else {
        break;
      }
    }

    fastx_parser::FastxParser<fastx_parser::ReadPair> rparser(
      po.single_read_filenames, po.barcode_filenames, nthread, np);

    rparser.start();
    boost::concurrent_flat_map<uint64_t, sshash::lookup_result> unitig_end_cache(unitig_end_cache_size);
    std::vector<std::thread> workers;
    for (size_t i = 0; i < nthread; ++i) {
      workers.push_back(std::thread([&ri, &po, &rparser, &binning, &ptab,
                                     &global_nr, &global_nh, &global_nmult,
                                     &k_match, &global_np, &out_info, &iomut,
                                     &rw, &l_match, &r_match, &dove_match,
                                     &dove_num, &ov_num, &ov_match, &r_orphan,
                                     &l_orphan, &unitig_end_cache]() {
        const auto token = rw.get_token();
        if (!po.enable_structural_constraints) {
          using SketchHitT =
            mapping::util::sketch_hit_info_no_struct_constraint;
          if (po.use_sam_format) {
            do_map<FragmentT, SketchHitT, SamT>(
              ri, rparser, binning, ptab, global_nr, global_nh, global_nmult,
              k_match, l_match, r_match, dove_num, dove_match, ov_num, ov_match,
              r_orphan, l_orphan, global_np, out_info, iomut, po.use_bed_format,
              po.check_kmers_orphans, po.tn5_shift, po.use_chr, unitig_end_cache, rw, token);
          } else {
            do_map<FragmentT, SketchHitT, RadT>(
              ri, rparser, binning, ptab, global_nr, global_nh, global_nmult,
              k_match, l_match, r_match, dove_num, dove_match, ov_num, ov_match,
              r_orphan, l_orphan, global_np, out_info, iomut, po.use_bed_format,
              po.check_kmers_orphans, po.tn5_shift, po.use_chr, unitig_end_cache, rw, token);
          }
        } else {
          using SketchHitT = mapping::util::sketch_hit_info;
          if (po.use_sam_format) {
            do_map<FragmentT, SketchHitT, SamT>(
              ri, rparser, binning, ptab, global_nr, global_nh, global_nmult,
              k_match, l_match, r_match, dove_num, dove_match, ov_num, ov_match,
              r_orphan, l_orphan, global_np, out_info, iomut, po.use_bed_format,
              po.check_kmers_orphans, po.tn5_shift, po.use_chr, unitig_end_cache, rw, token);
          } else {
            do_map<FragmentT, SketchHitT, RadT>(
              ri, rparser, binning, ptab, global_nr, global_nh, global_nmult,
              k_match, l_match, r_match, dove_num, dove_match, ov_num, ov_match,
              r_orphan, l_orphan, global_np, out_info, iomut, po.use_bed_format,
              po.check_kmers_orphans, po.tn5_shift, po.use_chr, unitig_end_cache, rw, token);
          }
        }
      }));
    }

    for (auto &w : workers) {
      w.join();
    }
    rparser.stop();
  }
  spdlog_piscem::info("finished mapping.");
  rw.close();

  if (po.use_sam_format) {
    out_info.sam_file.close();

    if (!out_info.sam_file) {
      spdlog_piscem::critical("The SAM file stream had an invalid status after "
                              "close; so some operation(s) may"
                              "have failed!\nA common cause for this is lack "
                              "of output disk space.\n"
                              "Consequently, the output may be corrupted.\n\n");
      return 1;
    }
  }

  if (po.use_bed_format) {
    out_info.bed_file.close();

    if (!out_info.bed_file) {
      spdlog_piscem::critical("The BED file stream had an invalid status after "
                              "close; so some operation(s) may"
                              "have failed!\nA common cause for this is lack "
                              "of output disk space.\n"
                              "Consequently, the output may be corrupted.\n\n");
      return 1;
    }
  }

  out_info.unmapped_bc_file.close();

  auto end_t = std::chrono::high_resolution_clock::now();
  auto num_sec =
    std::chrono::duration_cast<std::chrono::seconds>(end_t - start_t);
  piscem::meta_info::run_stats rs;
  rs.cmd_line(cmdline);
  rs.mode(piscem::RunMode::scatac);
  rs.num_reads(global_nr.load());
  rs.num_hits(global_nh.load());
  rs.num_multihits(global_nmult.load());
  rs.num_seconds(num_sec.count());
  rs.num_kmatch(k_match.load());
  rs.num_lmatch(l_match.load());
  rs.num_rmatch(r_match.load());
  rs.num_dovenum(dove_num.load());
  rs.num_dovematch(dove_match.load());
  rs.num_ovnum(ov_num.load());
  rs.num_ovmatch(ov_match.load());
  rs.num_rorphan(r_orphan.load());
  rs.num_lorphan(l_orphan.load());

  ghc::filesystem::path map_info_file_path = output_path / "map_info.json";
  bool info_ok = piscem::meta_info::write_map_info(rs, map_info_file_path);
  if (!info_ok) {
    spdlog_piscem::critical("failed to write map_info.json file");
  }

  return 0;
}
