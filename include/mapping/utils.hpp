#pragma once

#include "../include/CanonicalKmerIterator.hpp"
#include "../include/FastxParser.hpp"
#include "../include/hit_searcher.hpp"
#include "../include/itlib/small_vector.hpp"
#include "../include/parallel_hashmap/phmap.h"
#include "../include/poison_table.hpp"
#include "../include/projected_hits.hpp"
#include "../include/streaming_query.hpp"
#include "../include/reference_index.hpp"


#include <algorithm>
#include <cassert>
#include <cmath> // for std::ceil on linux
#include <fstream>
#include <limits>
#include <numeric>
#include <type_traits>
#include <vector>

namespace mapping {

namespace util {

class bin_pos {
    public:
        static constexpr uint64_t invalid_bin_id{std::numeric_limits<uint64_t>::max()};

        explicit bin_pos(mindex::reference_index* pfi,
                        float thr = 0.7,
                        uint64_t bin_size = 20000,
                        uint64_t overlap = 300
                        ): thr(thr), bin_size(bin_size), overlap(overlap) { 
            pfi_ = pfi,
            thr = thr,
            bin_size = bin_size,
            overlap = overlap,
            compute_cum_rank();
        };

        uint64_t get_cum_bin(size_t i) const {
            return cum_bin_ids[i];
        }

        float get_thr() const {
            return thr;
        }

        uint64_t get_bin_size() const {
            return bin_size;
        }

        uint64_t get_overlap() const {
            return overlap;
        }
        
        inline bool is_valid(uint64_t bin_id) const {
          return bin_id != invalid_bin_id;
        }

        inline std::pair<uint64_t, uint64_t> get_bin_id(uint64_t tid, uint64_t pos) const {
            uint64_t first_bin_id = cum_bin_ids[tid];
            uint64_t first_bin_id_in_next = cum_bin_ids[tid+1];
            assert(bin_size > overlap);

            uint64_t rel_pos = pos;
            uint64_t rel_bin = rel_pos / bin_size;
            uint64_t bin1 = first_bin_id + rel_bin;

            bool bin2_on_same_ref = (bin1+1) < first_bin_id_in_next;
            uint64_t bin2_rel_start_pos = ((rel_bin+1) * bin_size);
            // std::cout << "bin2 start pos " << bin2_rel_start_pos << " bin1 " << bin1 << std::endl;
            // `invalid_bin_id indicates` that the kmer does not belong to the overlapping region
            uint64_t bin2 = (bin2_on_same_ref and (rel_pos > (bin2_rel_start_pos - overlap))) ? (bin1+1) : invalid_bin_id; 
            return {bin1, bin2};
        }

        mindex::reference_index* get_ref() { return pfi_;}

    private:
        mindex::reference_index* pfi_;
        std::vector<uint64_t> cum_bin_ids;

        float thr;
        uint64_t bin_size;
        uint64_t overlap;

        void compute_cum_rank() {
            size_t n_refs = static_cast<size_t>(pfi_->num_refs());
            cum_bin_ids.resize(n_refs+1);
            cum_bin_ids[0] = 0;
            for(size_t i = 1; i <= n_refs; i++) {
                uint64_t rlen = static_cast<uint64_t>(pfi_->ref_len(i-1));
                uint64_t bins_in_ref = rlen / bin_size;
                uint64_t last_bin_end_pos = bins_in_ref * bin_size;
                // if there is space left after the bin before the 
                // end of the reference, then add another bin.
                if (last_bin_end_pos < rlen) {
                  bins_in_ref += 1;
                }
                cum_bin_ids[i] = cum_bin_ids[(i-1)] + bins_in_ref;
            }
        }
};

constexpr int32_t invalid_frag_len = std::numeric_limits<int32_t>::min();
constexpr int32_t invalid_mate_pos = std::numeric_limits<int32_t>::min();

enum class orientation_filter : uint8_t { NONE, FORWARD_ONLY, RC_ONLY };

struct simple_hit {
  bool is_fw{false};
  bool mate_is_fw{false};
  int32_t pos{-1};
  float score{0.0};
  uint32_t num_hits{0};
  uint32_t tid{std::numeric_limits<uint32_t>::max()};

  int32_t mate_pos{std::numeric_limits<int32_t>::max()};
  int32_t fragment_length{std::numeric_limits<int32_t>::max()};
  uint64_t bin_id{std::numeric_limits<uint64_t>::max()};

    // bool operator<(const mapping::util::simple_hit& a
    //                 ) {
    //   if (a.is_fw != is_fw) { return is_fw > a.is_fw; }
    //   // orientations are the same
    //   if (a.tid != tid) { return tid < a.tid; }
    //   // orientations & txps are the same
    //   if (a.pos != pos) { return pos < a.pos; }
    //   // orientations, txps, and positions are the same
    //   if (a.num_hits != num_hits ) { return num_hits < a.num_hits ; }
    //   // orientations, txps, positions and num_hits are the same 
    //   return bin_id < a.bin_id;
    // }

  inline bool valid_pos(int32_t read_len, uint32_t txp_len, int32_t max_over) {
    int32_t signed_txp_len = static_cast<int32_t>(txp_len);
    return (pos > -max_over) and
           ((pos + read_len) < (signed_txp_len + max_over));
  }
  inline bool has_mate() const { return mate_pos != invalid_mate_pos; }
  inline bool mate_is_mapped() const { return mate_pos != invalid_mate_pos; }
  inline int32_t frag_len() const {
    return (fragment_length != invalid_frag_len) ? fragment_length : 0;
  }
};

inline void print_hits(const std::vector<mapping::util::simple_hit> &hits ) {
    for (const auto& hit : hits) {
        std::cout << "isFw: " << hit.is_fw << std::endl;
        std::cout << "pos: " << hit.pos << std::endl;
        std::cout << "mate pos: " << hit.mate_pos << std::endl;
        std::cout << "num hits: " << hit.num_hits << std::endl;
        std::cout << "tid: " << hit.tid << std::endl;
        std::cout << "bin_id: " << hit.bin_id << std::endl;
        std::cout << "frag_len: " << hit.fragment_length << std::endl;
        std::cout << "------------------------" << std::endl;
    }
}
enum class MappingType : uint8_t {
  UNMAPPED = 0,
  SINGLE_MAPPED = 1,
  MAPPED_FIRST_ORPHAN = 2,
  MAPPED_SECOND_ORPHAN = 3,
  MAPPED_PAIR = 4
};

enum class HitDirection : uint8_t { FW, RC, BOTH };

constexpr uint8_t max_distortion = std::numeric_limits<uint8_t>::max();

struct chain_state {
  int32_t read_start_pos{-1};
  int32_t prev_pos{-1};
  int32_t curr_pos{-1};
  uint8_t num_hits{0};
  uint8_t min_distortion{max_distortion};
};

inline bool compare_chains(const chain_state &a, const chain_state &b) {
  return a.prev_pos < b.prev_pos;
}

struct sketch_hit_info {
  static constexpr size_t max_num_chains = 8;
  // add a hit to the current target that occurs in the forward
  // orientation with respect to the target.
  bool add_fw(int32_t ref_pos, int32_t read_pos, int32_t rl, int32_t k,
              int32_t max_stretch, float score_inc) {
    (void)rl;
    (void)k;
    bool added{false};

    int32_t approx_map_pos = ref_pos - read_pos;
    // If structural constraints have been disabled
    // then simply count the number of hits we see in
    // the given orientation (being careful to count
    // a k-mer of a given rank only one time).
    if (ignore_struct_constraints_fw) {
      if (read_pos > last_read_pos_fw) {
        if (last_read_pos_fw == -1) {
          approx_pos_fw = approx_map_pos;
        }
        last_ref_pos_fw = ref_pos;
        last_read_pos_fw = read_pos;
        fw_score += score_inc;
        ++fw_hits;
        added = true;
      }
      return added;
    }
    // NO STRUCTURAL CONSTRAINTS

    // If this is a k-mer of a new rank
    if (read_pos > last_read_pos_fw) {
      // if this is the first k-mer we are seeing
      if (last_read_pos_fw == -1) {
        approx_pos_fw = approx_map_pos;
      } else {
        // we are seeing a k-mer of a new rank.
        // at this point, we've seen all k-mers of the
        // previous rank, so copy over any valid chains
        // for the next search.
        compact_chains(fw_chains, fw_hits);
      }

      // update the current query position
      // (i.e. rank) of the seed we are processing
      last_read_pos_fw = read_pos;
      ++fw_rank;
    }

    // if this is still a hit for the *first*
    // k-mer of the chains — *note* the first
    // k-mer will be of rank 0 because fw_rank
    // is initialized to -1.
    if (fw_rank == 0) {
      process_rank0_hit(approx_map_pos, ref_pos, fw_chains, approx_pos_fw,
                        ignore_struct_constraints_fw, fw_hits, added);
    } else {
      // this is a hit for a k-mer of rank > 0, so we already
      // have a set of active chains.
      process_hit(true, approx_map_pos, ref_pos, max_stretch, fw_chains,
                  fw_hits, added);
    }
    return added;

    ///// orig
    /*
    // since hits are collected by moving _forward_ in the
    // read, if this is a fw hit, it should be moving
    // forward in the reference. Only add it if this is
    // the case.  This ensure that we don't
    // double-count a k-mer that might occur twice on
    // this target.
    if (ref_pos > last_ref_pos_fw and read_pos > last_read_pos_fw) {
        if (last_read_pos_fw == -1) {
            approx_pos_fw = ref_pos - read_pos;
        } else {
            if ((ref_pos - approx_pos_fw) > max_stretch) { return false; }
        }
        last_ref_pos_fw = ref_pos;
        last_read_pos_fw = read_pos;
        fw_score += score_inc;
        ++fw_hits;
        added = true;
    }
    return added;
    */
  }

  // add a hit to the current target that occurs in the forward
  // orientation with respect to the target.
  bool add_rc(int32_t ref_pos, int32_t read_pos, int32_t rl, int32_t k,
              int32_t max_stretch, float score_inc) {
    bool added{false};
    int32_t approx_map_pos = (ref_pos - (rl - (read_pos + k)));

    // NO STRUCTURAL CONSTRAINTS
    if (ignore_struct_constraints_rc) {
      if (read_pos > last_read_pos_rc) {
        approx_pos_rc = approx_map_pos;
        if (last_read_pos_rc == -1) {
          approx_end_pos_rc = ref_pos + read_pos;
          first_read_pos_rc = read_pos;
        }
        rc_score += score_inc;
        ++rc_hits;

        // new
        rightmost_bound_rc = last_ref_pos_rc;

        last_ref_pos_rc = ref_pos;
        last_read_pos_rc = read_pos;
        added = true;
      }
      return added;
    }
    // NO STRUCTURAL CONSTRAINTS

    // If this is a k-mer of a new rank
    if (read_pos > last_read_pos_rc) {
      approx_pos_rc = approx_map_pos;
      // if this is the first k-mer we are seeing
      if (last_read_pos_rc == -1) {
        approx_end_pos_rc = ref_pos + read_pos;
        first_read_pos_rc = read_pos;
      } else {
        // point, we've seen all k-mers of the previous rank, so copy over any
        // valid chains for the next search.
        compact_chains(rc_chains, rc_hits);
      }

      // update the current query position
      // (i.e. rank) of the seed we are processing
      ++rc_rank;
      // new
      rightmost_bound_rc = last_ref_pos_rc;

      last_ref_pos_rc = ref_pos;
      last_read_pos_rc = read_pos;
    }

    // if this is still a hit for the *first*
    // k-mer of the chains
    if (rc_rank == 0) {
      process_rank0_hit(approx_map_pos, ref_pos, rc_chains, approx_pos_rc,
                        ignore_struct_constraints_rc, rc_hits, added);
    } else {
      // this is a hit for a k-mer of rank > 0, so we already
      // have a set of active chains.
      process_hit(false, approx_map_pos, ref_pos, max_stretch, rc_chains,
                  rc_hits, added);
      if (added) {
        approx_pos_rc = approx_map_pos;
      }
    }
    return added;
    /*
    // since hits are collected by moving _forward_ in the
    // read, if this is an rc hit, it should be moving
    // backwards in the reference.
    // In general, we only add the hit if this is the case.
    // This ensures that we don't double-count a k-mer that
    // might occur twice on this target.

    // we have a special case here; what if the same exact
    // k-mer (i.e. not just the same sequence but same position
    // on the query) occurs more than one time on this refernece?
    //
    // In that case, the GENERAL case code will already have
    // processed and seen a k-mer with the read position
    // equal to `read_pos`.  In the case below, we see
    // a hit with the *same* read pos again (one or more times).
    //
    // Here, we swap out the previous hit having this read_pos
    // if the position of the current hit on the read is
    // the same and the position on the reference is greater
    // (this is a heuristic to help in the case of tandem repeats or
    // highly-repetitive subsequence).
    // NOTE: consider if a similar heuristic should be
    // adopted for the forward case.
    if ((read_pos == last_read_pos_rc) and (ref_pos > last_ref_pos_rc) and
        (ref_pos < rightmost_bound_rc)) {

        last_ref_pos_rc = ref_pos;
        // if the read_pos was the same as the first read pos
        // then also update the approx_end_pos_rc accordingly
        // NOTE: for the time being don't mess with this position
        // empirically this does better, but if we really want
        // to optimize this for accuracy we need a better general
        // heuristic.
        // if (read_pos == first_read_pos_rc) {
        //   approx_end_pos_rc = ref_pos + read_pos;
        //}

        return added;
    }

    // GENERAL case
    if (ref_pos < last_ref_pos_rc and read_pos > last_read_pos_rc) {
        approx_pos_rc = (ref_pos - (rl - (read_pos + k)));
        if (last_read_pos_rc == -1) {
            approx_end_pos_rc = ref_pos + read_pos;
            first_read_pos_rc = read_pos;
        } else {
            if (approx_end_pos_rc - approx_pos_rc > max_stretch) { return false;
    }
        }
        rc_score += score_inc;
        ++rc_hits;

        // new
        rightmost_bound_rc = last_ref_pos_rc;

        last_ref_pos_rc = ref_pos;
        last_read_pos_rc = read_pos;
        added = true;
    }
    return added;
    */
  }

  // for directly incrementing the number of hits
  // even when we are not building chains (e.g. in the case
  // of filtering based on occurrences of ambiguous seeds).
  inline void inc_fw_hits() { ++fw_hits; }
  inline void inc_rc_hits() { ++rc_hits; }

  inline uint32_t max_hits_for_target() { return std::max(fw_hits, rc_hits); }

  // true if forward, false if rc
  // second element is score
  inline HitDirection best_hit_direction() {
    int32_t fw_minus_rc =
      static_cast<int32_t>(fw_hits) - static_cast<int32_t>(rc_hits);
    return (fw_minus_rc > 0)
             ? HitDirection::FW
             : ((fw_minus_rc < 0) ? HitDirection::RC : HitDirection::BOTH);
  }

  inline simple_hit get_fw_hit() {
    return simple_hit{true,     false,   approx_pos_fw,
                      fw_score, fw_hits, std::numeric_limits<uint32_t>::max()};
  }

  inline simple_hit get_rc_hit() {
    return simple_hit{false,    false,   approx_pos_rc,
                      rc_score, rc_hits, std::numeric_limits<uint32_t>::max()};
  }

  inline simple_hit get_best_hit() {
    auto best_direction = best_hit_direction();
    return (best_direction != HitDirection::RC)
             ? simple_hit{true,          false,
                          approx_pos_fw, fw_score,
                          fw_hits,       std::numeric_limits<uint32_t>::max()}
             : simple_hit{false,         false,
                          approx_pos_rc, rc_score,
                          rc_hits,       std::numeric_limits<uint32_t>::max()};
  }

  inline std::string to_string() {
    std::stringstream ss;
    ss << "fw_hits: " << fw_hits << ", fw_score : " << fw_score
       << ", fw_pos : " << approx_pos_fw << " || rc_hits: " << rc_hits
       << ", rc_score: " << rc_score << ", rc_pos: " << approx_pos_rc;
    return ss.str();
  }

  int32_t last_read_pos_fw{-1};
  int32_t last_read_pos_rc{-1};
  int32_t rightmost_bound_rc{std::numeric_limits<int32_t>::max()};

  // marks the read position (key) of the
  // first hit we see in the rc direction
  int32_t first_read_pos_rc{-1};

  int32_t last_ref_pos_fw{-1};
  int32_t last_ref_pos_rc{std::numeric_limits<int32_t>::max()};

  int32_t approx_pos_fw{-1};
  int32_t approx_pos_rc{-1};
  int32_t approx_end_pos_rc{-1};

  uint32_t fw_hits{0};
  uint32_t rc_hits{0};
  float fw_score{0.0};
  float rc_score{0.0};

  bool ignore_struct_constraints_fw{false};
  bool ignore_struct_constraints_rc{false};

  int32_t tid{std::numeric_limits<int32_t>::max()};
  int32_t fw_rank{-1};
  int32_t rc_rank{-1};
  itlib::small_vector<chain_state, max_num_chains> fw_chains;
  itlib::small_vector<chain_state, max_num_chains> rc_chains;

private:
  inline void
  compact_chains(itlib::small_vector<chain_state, max_num_chains> &chains,
                 const uint32_t required_hits) {
    chains.erase(std::remove_if(chains.begin(), chains.end(),
                                [required_hits](chain_state &s) -> bool {
                                  // remove this chain if it doesn't satisfy
                                  // the hit constraint.
                                  if (s.num_hits < required_hits) {
                                    return true;
                                  }
                                  // the current position becomes the
                                  // previous position for kmers of the
                                  // next rank, and the curr pos gets reset
                                  // to -1.
                                  s.prev_pos = s.curr_pos;
                                  s.curr_pos = -1;
                                  return false;
                                }),
                 chains.end());
  }

  inline void
  process_rank0_hit(int32_t approx_map_pos, int32_t hit_pos,
                    itlib::small_vector<chain_state, max_num_chains> &chains,
                    int32_t &approx_pos_out, bool &ignore_struct_constraints,
                    uint32_t &num_hits, bool &added) {
    // if there are too many possible chains, just punt
    // and turn off structural constraints for this
    // query, reference, orientation tuple.
    if (chains.size() == chains.capacity()) {
      ignore_struct_constraints = true;
      approx_pos_out = chains.front().read_start_pos;
      num_hits = chains.front().num_hits;
      added = true;
    } else {
      // otherise add this hit
      chains.push_back({approx_map_pos, -1, hit_pos, 1, max_distortion});
      added = true;
    }
    num_hits = 1;
  }

  inline void
  process_hit(bool is_fw_hit, // we are processing a hit in the forward
                              // orientation (otherwise, RC)
              int32_t read_start_pos, int32_t next_hit_pos, int32_t max_stretch,
              itlib::small_vector<chain_state, max_num_chains> &chains,
              uint32_t &num_hits, bool &added) {
    (void)max_stretch;
    // find the chain that best matches this k-mer.
    chain_state predecessor_probe{read_start_pos, next_hit_pos, -1, 0,
                                  max_distortion};
    auto chain_pos = std::lower_bound(chains.begin(), chains.end(),
                                      predecessor_probe, compare_chains);

    // If this is a forward hit, lower bound will return
    // the first element >= the key (the current hit), so back up
    // by one element.
    // Otherwise, if this is a reverse complement hit, we actually
    // want the first element >= the current hit, so keep it.
    if (is_fw_hit) {
      if (chain_pos > chains.begin()) {
        chain_pos--;
      } else {
        added = false;
        return;
      }
    }

    // FOR DEBUG
    // std::cerr << " === examining hit at ref pos : " << next_hit_pos << ",
    // lower bound search found chain ending at : " << chain_pos->prev_pos << "
    // =========\n";

    constexpr int32_t max_chain_stretch = 31;
    // if we found a valid chain to extend
    int32_t tries = 0;
    // if we don't find a valid chain in the expected place, look at *one* more
    // predecessor (or successor in the reverse complement case) to see if we
    // can extend that chain. NOTE: Think if # of tries should be an advanced
    // parameter available to the user.
    while ((chain_pos < chains.end()) and (chain_pos >= chains.begin()) and
           !added and (tries < 2)) {
      auto stretch = std::abs(chain_pos->read_start_pos - read_start_pos);
      stretch =
        std::min(stretch, static_cast<decltype(stretch)>(max_distortion));
      // and if it hasn't yet been extended
      if (chain_pos->curr_pos == -1) {
        if (stretch < max_chain_stretch) {
          // then extend this chain
          chain_pos->curr_pos = next_hit_pos;
          chain_pos->min_distortion = static_cast<uint8_t>(stretch);
          // and increment the number of hits
          ++(chain_pos->num_hits);
          uint32_t curr_max_hits = num_hits;
          num_hits =
            std::max(curr_max_hits,
                     static_cast<decltype(curr_max_hits)>(chain_pos->num_hits));
          added = true;
        }
      } else if (stretch < chain_pos->min_distortion) {
        chain_pos->min_distortion = static_cast<uint8_t>(stretch);
        added = true;
      }
      // the chain to check for forward strand hits is the predecessor
      // and for reverse complement strand hits it's the successor.
      if (is_fw_hit) {
        chain_pos--;
      } else {
        chain_pos++;
      }
      tries++;
    }
  }
};

//
//
// The sketch hit info we want to maintain when we are
// not interested in
//
//

struct sketch_hit_info_no_struct_constraint {
  // add a hit to the current target that occurs in the forward
  // orientation with respect to the target.
  bool add_fw(int32_t ref_pos, int32_t read_pos, int32_t rl, int32_t k,
              int32_t max_stretch, float score_inc) {
    (void)rl;
    (void)k;
    (void)max_stretch;
    bool added{false};

    int32_t approx_map_pos = ref_pos - read_pos;
    // Simply count the number of hits we see in
    // the given orientation (being careful to count
    // a k-mer of a given rank only one time).
    if (read_pos > last_read_pos_fw) {
      if (last_read_pos_fw == -1) {
        approx_pos_fw = approx_map_pos;
      }
      last_ref_pos_fw = ref_pos;
      last_read_pos_fw = read_pos;
      fw_score += score_inc;
      ++fw_hits;
      added = true;
    }
    return added;
    // NO STRUCTURAL CONSTRAINTS
  }

  // add a hit to the current target that occurs in the forward
  // orientation with respect to the target.
  bool add_rc(int32_t ref_pos, int32_t read_pos, int32_t rl, int32_t k,
              int32_t max_stretch, float score_inc) {
    (void)rl;
    (void)k;
    (void)max_stretch;
    bool added{false};
    int32_t approx_map_pos = (ref_pos - (rl - (read_pos + k)));

    if (read_pos > last_read_pos_rc) {
      approx_pos_rc = approx_map_pos;
      if (last_read_pos_rc == -1) {
        approx_end_pos_rc = ref_pos + read_pos;
        first_read_pos_rc = read_pos;
      }
      rc_score += score_inc;
      ++rc_hits;

      // new
      rightmost_bound_rc = last_ref_pos_rc;

      last_ref_pos_rc = ref_pos;
      last_read_pos_rc = read_pos;
      added = true;
    }
    return added;
  }

  // for directly incrementing the number of hits
  // even when we are not building chains (e.g. in the case
  // of filtering based on occurrences of ambiguous seeds).
  inline void inc_fw_hits() { ++fw_hits; }
  inline void inc_rc_hits() { ++rc_hits; }

  inline uint32_t max_hits_for_target() { return std::max(fw_hits, rc_hits); }

  // true if forward, false if rc
  // second element is score
  inline HitDirection best_hit_direction() {
    int32_t fw_minus_rc =
      static_cast<int32_t>(fw_hits) - static_cast<int32_t>(rc_hits);
    return (fw_minus_rc > 0)
             ? HitDirection::FW
             : ((fw_minus_rc < 0) ? HitDirection::RC : HitDirection::BOTH);
  }

  inline simple_hit get_fw_hit() {
    return simple_hit{true,     false,   approx_pos_fw,
                      fw_score, fw_hits, std::numeric_limits<uint32_t>::max()};
  }

  inline simple_hit get_rc_hit() {
    return simple_hit{false,    false,   approx_pos_rc,
                      rc_score, rc_hits, std::numeric_limits<uint32_t>::max()};
  }

  inline simple_hit get_best_hit() {
    auto best_direction = best_hit_direction();
    return (best_direction != HitDirection::RC)
             ? simple_hit{true,          false,
                          approx_pos_fw, fw_score,
                          fw_hits,       std::numeric_limits<uint32_t>::max()}
             : simple_hit{false,         false,
                          approx_pos_rc, rc_score,
                          rc_hits,       std::numeric_limits<uint32_t>::max()};
  }

  inline std::string to_string() {
    std::stringstream ss;
    ss << "fw_hits: " << fw_hits << ", fw_score : " << fw_score
       << ", fw_pos : " << approx_pos_fw << " || rc_hits: " << rc_hits
       << ", rc_score: " << rc_score << ", rc_pos: " << approx_pos_rc;
    return ss.str();
  }

  int32_t last_read_pos_fw{-1};
  int32_t last_read_pos_rc{-1};
  int32_t rightmost_bound_rc{std::numeric_limits<int32_t>::max()};

  // marks the read position (key) of the
  // first hit we see in the rc direction
  int32_t first_read_pos_rc{-1};

  int32_t last_ref_pos_fw{-1};
  int32_t last_ref_pos_rc{std::numeric_limits<int32_t>::max()};

  int32_t approx_pos_fw{-1};
  int32_t approx_pos_rc{-1};
  int32_t approx_end_pos_rc{-1};

  uint32_t fw_hits{0};
  uint32_t rc_hits{0};
  float fw_score{0.0};
  float rc_score{0.0};
  int32_t tid{std::numeric_limits<int32_t>::max()};

};

enum class fragment_end : uint8_t { LEFT, RIGHT };

struct poison_state_t {
  inline void clear() {
    poisoned_left = false;
    poisoned_right = false;
    fend = fragment_end::LEFT;
  }

  inline bool is_poisoned() const {
    if (paired_for_mapping) {
      return poisoned_left or poisoned_right;
    } else {
      return poisoned_left;
    }
  }

  // Returns true if the mapping was poisoned, false otherwise. This uses the
  // notion of "poison" k-mers (aka "distinguishing flanking k-mers" or DFKs[1])
  // to determine if a mapping should be discarded because it is contaminated by
  // decoy k-mers that might suggest a different mapping if other decoy sequence
  // was included in the index.
  //
  // [1] Hjörleifsson, K. E., Sullivan, D. K., Holley, G., Melsted, P. &
  // Pachter, L. Accurate quantification of single-nucleus and single-cell
  // RNA-seq transcripts. bioRxiv.
  // https://www.biorxiv.org/content/early/2022/12/02/2022. 12.02.518832 (2022)
  bool scan_raw_hits(std::string &s, uint32_t k, mindex::reference_index *index,
                     std::vector<std::pair<int, projected_hits>> &h,
                     mindex::SkippingStrategy strat) {
    // a read that didn't map can't be poisoned
    if (h.empty()) {
      return false;
    }

    constexpr bool verbose = false;
    bool strict_mode = (strat == mindex::SkippingStrategy::STRICT);
    bool was_poisoned = false;
    auto first_pos = h.front().first;
    pufferfish::CanonicalKmerIterator kit_end;

    pufferfish::CanonicalKmerIterator kit(s);

    auto &last_phit = h.back().second;
    uint32_t last_uni = last_phit.contig_id();
    int64_t dist_to_contig_end = -1;
    // fw ori
    if (last_phit.hit_fw_on_contig()) {
      dist_to_contig_end = static_cast<int64_t>(last_phit.contig_len()) -
                           (static_cast<int64_t>(last_phit.contig_pos() + k));
    } else { // rc ori
      dist_to_contig_end = static_cast<int64_t>(last_phit.contig_pos());
    }
    auto terminal_pos = h.back().first + dist_to_contig_end;

    if constexpr (verbose) {
      for (auto &hit : h) {
        std::cerr << "read_pos: " << hit.first << ", phit: " << hit.second
                  << "\n";
        auto &refs = hit.second.refRange;

        for (auto v : refs) {
          const auto &ref_pos_ori = hit.second.decode_hit(v);
          uint32_t tid = sshash::util::transcript_id(v);
          int32_t pos = static_cast<int32_t>(ref_pos_ori.pos);
          bool ori = ref_pos_ori.isFW;

          if constexpr (verbose) {
            auto &tname = index->ref_name(tid);
            std::cerr << "kmer: " << s.substr(hit.first, k) << ", tid: " << tid
                      << " (name: " << tname << "), pos: " << pos
                      << ", dir: " << (ori ? "fw" : "rc") << "\n";
          }
        }
      }
    }

    // scan up to the first hit looking for any poison k-mer
    while ((kit != kit_end) and (kit->second < first_pos)) {
      if (ptab->key_exists(kit->first.getCanonicalWord())) {
        was_poisoned = true;
        if constexpr (verbose) {
          std::cerr << "[[[was poisoned (" << kit->second << ","
                    << s.substr(kit->second, k) << ")]]]\n";
        }
        return was_poisoned;
      }
      ++kit;
    }

    // at this point, we got to the first hit on the read
    // now we scan intervals.
    auto start_it = h.begin();
    auto end_it = start_it + 1;
    int last_pos = 0;
    (void)last_pos;
    (void)last_uni;
    while ((end_it != h.end()) and (kit != kit_end)) {
      // the first unitig to which the poison kmer can belong
      auto u1 = start_it->second.contig_id();
      auto u2 = end_it->second.contig_id();

      auto cp1 = start_it->second.contig_pos();
      auto cp2 = end_it->second.contig_pos();

      // auto p1 = start_it->first;
      auto p2 = end_it->first;
      bool right_bound_resulted_from_open_search =
        end_it->second.resulted_from_open_search;
      last_pos = p2;
      last_uni = u2;

      auto lb = std::min(cp1, cp2);
      lb = (lb > 0) ? lb + 1 : 0;
      auto ub = std::max(cp1, cp2);
      ub = (ub < start_it->second.contigLen_ - k)
             ? ub - 1
             : start_it->second.contigLen_ - k;

      while ((kit != kit_end) and (kit->second < p2)) {
        if (!kit->first.is_low_complexity()) {
          if (!right_bound_resulted_from_open_search) {
            // we shouldn't even have to check in strict mode.
            was_poisoned = (strict_mode)
                             ? false
                             : ptab->key_occurs_in_unitig_between(
                                 kit->first.getCanonicalWord(), u1, lb, ub);
          } else {
            was_poisoned = ptab->key_exists(kit->first.getCanonicalWord());
          }
          if (was_poisoned) {
            if constexpr (verbose) {
              std::cerr << "[[[was poisoned (" << kit->second << ","
                        << s.substr(kit->second, k) << ")]]]\n";
            }
            return was_poisoned;
          }
        }
        ++kit;
      }
      ++start_it;
      ++end_it;
    }

    // for any remaining k-mers in the read after the end of the last
    // matching interval.
    while (kit != kit_end) {
      // was_poisoned = kit->first.is_low_complexity() ? false :
      // ptab->key_exists(kit->first.getCanonicalWord());//
      was_poisoned =
        (kit->second < terminal_pos)
          ? ptab->key_occurs_in_unitig(kit->first.getCanonicalWord(), last_uni)
          : ptab->key_exists(kit->first.getCanonicalWord());
      if (was_poisoned) {
        if constexpr (verbose) {
          ptab->print_occs(kit->first.getCanonicalWord());
        }
        break;
      }
      ++kit;
    }
    if constexpr (verbose) {
      if (was_poisoned) {
        std::cerr << "[[[was poisoned (" << kit->second << ","
                  << s.substr(kit->second, k) << ")]]]\n";
      } else {
        std::cerr << "\n\n[[[was not poisoned]]]\n\n";
      }
    }
    return was_poisoned;
  }

  // returns true if this poison state object is valid
  // (i.e. if there is a valid posion table associated with
  // it) and false otherwise.
  inline bool is_valid() const { return ptab != nullptr; }

  // poison whichever fragment end is currently active.
  inline void poison_read() {
    if (fend == fragment_end::LEFT) {
      poisoned_left = true;
    } else {
      poisoned_right = true;
    }
  }

  // set the active fragment end (either LEFT or RIGHT)
  inline void set_fragment_end(fragment_end e) { fend = e; }
  // get the current fragment end (either LEFT or RIGHT)
  inline fragment_end get_fragment_end() const { return fend; }

  fragment_end fend{fragment_end::LEFT};
  bool poisoned_left{false};
  bool poisoned_right{false};
  bool paired_for_mapping{false};
  poison_table *ptab{nullptr};
};

template <typename sketch_hit_info_t> struct mapping_cache_info {
public:
  mapping_cache_info(mindex::reference_index &ri)
    : k(ri.k()), q(ri.get_dict()), hs(&ri) {}

  inline void clear() {
    map_type = mapping::util::MappingType::UNMAPPED;
    q.start();
    hs.clear();
    hit_map.clear();
    accepted_hits.clear();
    has_matching_kmers = false;
    ambiguous_hit_indices.clear();
    frag_seq = "";
    read1 = true;
  }

  // will store how the read mapped
  mapping::util::MappingType map_type{mapping::util::MappingType::UNMAPPED};

  // map from reference id to hit info
  phmap::flat_hash_map<uint32_t, sketch_hit_info_t> hit_map;
  std::vector<mapping::util::simple_hit> accepted_hits;

  // map to recall the number of unmapped reads we see
  // for each barcode
  phmap::flat_hash_map<uint64_t, uint32_t> unmapped_bc_map;

  size_t max_hit_occ = 256;
  size_t max_hit_occ_recover = 1024;
  bool attempt_occ_recover = (max_hit_occ_recover > max_hit_occ);
  size_t max_read_occ = 2500;
  size_t k{0};

  // to perform queries
  piscem::streaming_query q;
  // implements the PASC algorithm
  mindex::hit_searcher hs;
  size_t max_chunk_reads = 5000;
  // regardless of having full mappings, did any k-mers match
  bool has_matching_kmers{false};
  // holds the indices of k-mers too ambiguous to chain, but which
  // we might later want to check the existence of
  itlib::small_vector<uint32_t, 255> ambiguous_hit_indices;
  // max ec card
  uint32_t max_ec_card{4096};
  // Paired end reads when trying to concatenate before mapping
  std::string frag_seq = "";
  bool read1 = true;
};

template <typename mapping_cache_info_t>
inline bool
map_read(std::string *read_seq, mapping_cache_info_t &map_cache,
         poison_state_t &poison_state,
         mindex::SkippingStrategy strat = mindex::SkippingStrategy::STRICT,
         bool verbose = false) {
  map_cache.clear();
  // rebind map_cache variables to
  // local names
  auto &q = map_cache.q;
  auto &hs = map_cache.hs;
  auto &hit_map = map_cache.hit_map;
  auto &accepted_hits = map_cache.accepted_hits;
  auto &map_type = map_cache.map_type;
  const bool attempt_occ_recover = map_cache.attempt_occ_recover;
  const bool perform_ambig_filtering = map_cache.hs.get_index()->has_ec_table();
  auto k = map_cache.k;
  bool apply_poison_filter = poison_state.is_valid();

  map_cache.has_matching_kmers =
    hs.get_raw_hits_sketch(*read_seq, q, strat, true, false);
  bool early_stop = false;

  // if we are checking ambiguous hits, the maximum EC
  // size we will consider.
  const size_t max_ec_ambig = map_cache.max_ec_card;

  // if there were hits
  if (map_cache.has_matching_kmers) {
    uint32_t num_valid_hits{0};
    uint64_t total_occs{0};
    uint64_t largest_occ{0};
    auto &raw_hits = hs.get_left_hits();

    // if we are applying a poison filter, do it here.
    if (apply_poison_filter) {
      bool was_poisoned = poison_state.scan_raw_hits(
        *read_seq, k, map_cache.hs.get_index(), raw_hits, strat);
      if (was_poisoned) {
        poison_state.poison_read();
        map_type = mapping::util::MappingType::UNMAPPED;
        return true;
      }
    }

    // SANITY
    decltype(raw_hits[0].first) prev_read_pos = -1;
    // the maximum span the supporting k-mers of a
    // mapping position are allowed to have.
    // NOTE this is still > read_length b/c the stretch is measured wrt the
    // START of the terminal k-mer.
    int32_t max_stretch = static_cast<int32_t>(read_seq->length() * 1.0);

    // a raw hit is a pair of read_pos and a projected hit

    // the least frequent hit for this fragment.
    uint64_t min_occ = std::numeric_limits<uint64_t>::max();

    // this is false by default and will be set to true
    // if *every* collected hit for this fragment occurs
    // max_hit_occ times or more.
    int32_t signed_rl = static_cast<int32_t>(read_seq->length());
    auto collect_mappings_from_hits =
      [&max_stretch, &min_occ, &hit_map, &num_valid_hits, &total_occs,
       &largest_occ, signed_rl, k, perform_ambig_filtering,
       /* &map_cache, */ // currently unused
       verbose](auto &raw_hits, auto &prev_read_pos, auto &max_allowed_occ,
                auto &ambiguous_hit_indices) -> bool {
      (void)verbose;
      int32_t hit_idx{0};
      hit_map.clear();

      for (auto &raw_hit : raw_hits) {
        auto &read_pos = raw_hit.first;
        auto &proj_hits = raw_hit.second;
        auto &refs = proj_hits.refRange;

        uint64_t num_occ = static_cast<uint64_t>(refs.size());
        min_occ = std::min(min_occ, num_occ);

        bool still_have_valid_target = false;
        prev_read_pos = read_pos;

        if (num_occ <= max_allowed_occ) {
          total_occs += num_occ;
          largest_occ = (num_occ > largest_occ) ? num_occ : largest_occ;
          float score_inc = 1.0;

          for (auto v : refs) {
            const auto &ref_pos_ori = proj_hits.decode_hit(v);
            uint32_t tid = sshash::util::transcript_id(v);
            int32_t pos = static_cast<int32_t>(ref_pos_ori.pos);
            bool ori = ref_pos_ori.isFW;
            auto &target = hit_map[tid];

            /* FOR DEBUG
            *
            if (true or verbose) {
                auto& tname = map_cache.hs.get_index()->ref_name(tid);
                std::cerr << "\traw_hit [read_pos: " << read_pos << " ]:" <<
            tname
                          << ", " << pos << ", " << (ori ? "fw" : "rc") << "\n";
            }
            */

            // Why >= here instead of == ?
            // Because hits can happen on the same target in both the forward
            // and rc orientations, it is possible that we start the loop with
            // the target having num_valid_hits hits in a given orientation (o)
            // we see a new hit for this target in oriention o (now it has
            // num_valid_hits + 1) then we see a hit for this target in
            // orientation rc(o).  We still want to add / consider this hit, but
            // max_hits_for_target() > num_valid_hits. So, we must allow for
            // that here.
            if (target.max_hits_for_target() >= num_valid_hits) {
              if (ori) {
                target.add_fw(pos, static_cast<int32_t>(read_pos), signed_rl, k,
                              max_stretch, score_inc);
              } else {
                target.add_rc(pos, static_cast<int32_t>(read_pos), signed_rl, k,
                              max_stretch, score_inc);
              }

              still_have_valid_target |=
                (target.max_hits_for_target() >= num_valid_hits + 1);
            }

          } // DONE: for (auto &pos_it : refs)
          ++num_valid_hits;

          // if there are no targets reaching the valid hit threshold, then
          // break early
          if (!still_have_valid_target) {
            return true;
          }

        } else if (perform_ambig_filtering) { // HERE we have that num_occ >
                                              // max_allowed_occ
          ambiguous_hit_indices.push_back(hit_idx);
        }

        ++hit_idx;
      } // DONE : for (auto& raw_hit : raw_hits)

      return false;
    };

    auto mao_first_pass = map_cache.max_hit_occ - 1;
    early_stop = collect_mappings_from_hits(
      raw_hits, prev_read_pos, mao_first_pass, map_cache.ambiguous_hit_indices);

    // If our default threshold was too stringent, then fallback to a more
    // liberal threshold and look up the k-mers that occur the least frequently.
    // Specifically, if the min occuring hits have frequency <
    // max_hit_occ_recover (2500 by default) times, then collect the min
    // occuring hits to get the mapping.
    if (attempt_occ_recover and (min_occ >= map_cache.max_hit_occ) and
        (min_occ < map_cache.max_hit_occ_recover)) {
      num_valid_hits = 0;
      map_cache.ambiguous_hit_indices.clear();
      prev_read_pos = -1;
      uint64_t max_allowed_occ = min_occ;
      early_stop =
        collect_mappings_from_hits(raw_hits, prev_read_pos, max_allowed_occ,
                                   map_cache.ambiguous_hit_indices);
    }

    // Further filtering of mappings by ambiguous k-mers
    if (perform_ambig_filtering and !hit_map.empty() and
        !map_cache.ambiguous_hit_indices.empty()) {
      phmap::flat_hash_set<uint64_t> observed_ecs;
      size_t min_cardinality_ec_size = std::numeric_limits<size_t>::max();
      uint64_t min_cardinality_ec = std::numeric_limits<size_t>::max();
      size_t min_cardinality_index = 0;
      size_t visited = 0;
      auto &ec_table = map_cache.hs.get_index()->get_ec_table();

      auto visit_ec = [&hit_map](uint64_t ent, bool fw_on_contig) -> bool {
        uint32_t tid = (ent >> 2);
        auto hm_it = hit_map.find(tid);
        bool found = false;
        if (hm_it != hit_map.end()) {
          // we found this target, now:
          // (1) check the orientation
          uint32_t ori = (ent & 0x3);
          // (2) add hits in the appropriate way
          switch (ori) {
          case 0: // fw
            (fw_on_contig) ? hm_it->second.inc_fw_hits()
                           : hm_it->second.inc_rc_hits();
            break;
          case 1: // rc
            (fw_on_contig) ? hm_it->second.inc_rc_hits()
                           : hm_it->second.inc_fw_hits();
            break;
          default: // both
            hm_it->second.inc_fw_hits();
            hm_it->second.inc_rc_hits();
          }
          found = true;
        }
        return found;
      };

      // for each ambiguous hit
      for (auto hit_idx : map_cache.ambiguous_hit_indices) {
        auto &proj_hit = raw_hits[hit_idx].second;
        uint32_t contig_id = proj_hit.contig_id();
        bool fw_on_contig = proj_hit.hit_fw_on_contig();

        // put the combination of the eq and the k-mer orientation
        // into the map.
        uint64_t ec = ec_table.ec_for_tile(contig_id);
        uint64_t ec_key = ec | (fw_on_contig ? 0 : 0x8000000000000000);
        // if we've already seen this ec, no point in processing
        // it again.
        if (observed_ecs.contains(ec_key)) {
          continue;
        }
        // otherwise, insert it.
        observed_ecs.insert(ec_key);

        auto ec_entries = ec_table.entries_for_ec(ec);
        if (ec_entries.size() < min_cardinality_ec_size) {
          min_cardinality_ec_size = ec_entries.size();
          min_cardinality_ec = ec;
          min_cardinality_index = hit_idx;
        }
        if (ec_entries.size() > max_ec_ambig) {
          continue;
        }
        ++visited;
        for (const auto ent : ec_entries) {
          visit_ec(ent, fw_on_contig);
        } // all target oritentation pairs in this eq class
        ++num_valid_hits;
      } // all ambiguous hits

      // if we haven't visited *any* equivalence classes (they were all)
      // too ambiguous, then make a last-ditch effort to just visit the
      // the one with smallest cardinality.
      if (visited == 0) {
        auto hit_idx = min_cardinality_index;
        auto &proj_hit = raw_hits[hit_idx].second;
        bool fw_on_contig = proj_hit.hit_fw_on_contig();

        uint64_t ec = min_cardinality_ec;
        auto ec_entries = ec_table.entries_for_ec(ec);
        for (const auto ent : ec_entries) {
          visit_ec(ent, fw_on_contig);
        } // all target oritentation pairs in this eq class
        ++num_valid_hits;
      } // done visiting the last-ditch ec

    } // if we are processing ambiguous hits

    uint32_t best_alt_hits = 0;
    // int32_t signed_read_len = static_cast<int32_t>(record.seq.length());
    for (auto &kv : hit_map) {
      auto best_hit_dir = kv.second.best_hit_direction();

      // if the best direction is FW or BOTH, add the fw hit
      // otherwise add the RC.
      auto simple_hit = (best_hit_dir != mapping::util::HitDirection::RC)
                          ? kv.second.get_fw_hit()
                          : kv.second.get_rc_hit();

      if (simple_hit.num_hits >= num_valid_hits) {
        simple_hit.tid = kv.first;
        accepted_hits.emplace_back(simple_hit);
        // if we had equally good hits in both directions
        // add the rc hit here (since we added the fw)
        // above if the best hit was either FW or BOTH
        if (best_hit_dir == mapping::util::HitDirection::BOTH) {
          auto second_hit = kv.second.get_rc_hit();
          second_hit.tid = kv.first;
          accepted_hits.emplace_back(second_hit);
        }
      } else {
        // best_alt_score = simple_hit.score > best_alt_score ? simple_hit.score
        // : best_alt_score;
        best_alt_hits = simple_hit.num_hits > best_alt_hits
                          ? simple_hit.num_hits
                          : best_alt_hits;
      }
    }

    // max_read_occ = had_max_read_occ ? accepted_hits.size() : max_hit_occ;

    /*
     * This rule; if enabled, allows through mappings missing a single hit, if
    there
     * was no mapping with all hits. NOTE: this won't work with the current
    early-exit
     * optimization however.
    if (accepted_hits.empty() and (num_valid_hits > 1) and (best_alt_hits >=
    num_valid_hits
    - 1)) { for (auto& kv : hit_map) { auto simple_hit =
    kv.second.get_best_hit(); if (simple_hit.num_hits >= best_alt_hits) {
          //if (simple_hit.valid_pos(signed_read_len,
    transcripts[kv.first].RefLength, 10)) { simple_hit.tid = kv.first;
    accepted_hits.emplace_back(simple_hit);
          //}
        }
      }
    }
    */
  } // DONE : if (rh)

  // If the read mapped to > maxReadOccs places, discard it
  if (accepted_hits.size() > map_cache.max_read_occ) {
    accepted_hits.clear();
    map_type = mapping::util::MappingType::UNMAPPED;
  } else if (!accepted_hits.empty()) {
    map_type = mapping::util::MappingType::SINGLE_MAPPED;
  }

  return early_stop;
}

template <typename mapping_cache_info_t>
inline bool
map_read(std::string *read_seq, mapping_cache_info_t &map_cache,
         poison_state_t &poison_state,
         const mapping::util::bin_pos& binning,
         bool &k_match,
         mindex::SkippingStrategy strat = mindex::SkippingStrategy::STRICT,
         bool verbose = false) {
  map_cache.clear();
  // rebind map_cache variables to
  // local names
  auto &q = map_cache.q;
  auto &hs = map_cache.hs;
  auto &hit_map = map_cache.hit_map;
  auto &accepted_hits = map_cache.accepted_hits;
  auto &map_type = map_cache.map_type;
  const bool attempt_occ_recover = map_cache.attempt_occ_recover;
  const bool perform_ambig_filtering = map_cache.hs.get_index()->has_ec_table();
  auto k = map_cache.k;
  bool apply_poison_filter = poison_state.is_valid();
  auto thr = binning.get_thr();

  map_cache.has_matching_kmers =
    hs.get_raw_hits_sketch_everykmer(*read_seq, q, true, false);
  
  bool early_stop = false;

  // if we are checking ambiguous hits, the maximum EC
  // size we will consider.
  const size_t max_ec_ambig = map_cache.max_ec_card;

  // if there were hits
  if (map_cache.has_matching_kmers) {
    k_match = true;
    uint32_t num_valid_hits{0};
    uint64_t total_occs{0};
    uint64_t largest_occ{0};
    auto &raw_hits = hs.get_left_hits();
    // if we are applying a poison filter, do it here.
    if (apply_poison_filter) {
      bool was_poisoned = poison_state.scan_raw_hits(
        *read_seq, k, map_cache.hs.get_index(), raw_hits, strat);
      if (was_poisoned) {
        poison_state.poison_read();
        map_type = mapping::util::MappingType::UNMAPPED;
        return true;
      }
    }

    // SANITY
    decltype(raw_hits[0].first) prev_read_pos = -1;
    // the maximum span the supporting k-mers of a
    // mapping position are allowed to have.
    // NOTE this is still > read_length b/c the stretch is measured wrt the
    // START of the terminal k-mer.
    int32_t max_stretch = static_cast<int32_t>(read_seq->length() * 1.0);

    // a raw hit is a pair of read_pos and a projected hit

    // the least frequent hit for this fragment.
    uint64_t min_occ = std::numeric_limits<uint64_t>::max();

    // this is false by default and will be set to true
    // if *every* collected hit for this fragment occurs
    // max_hit_occ times or more.
    int32_t signed_rl = static_cast<int32_t>(read_seq->length());
    auto collect_mappings_from_hits_thr =
      [&max_stretch, &min_occ, &hit_map, &num_valid_hits, &total_occs,
       &largest_occ, &binning, signed_rl, k, perform_ambig_filtering, thr,
       verbose](auto &raw_hits, auto &prev_read_pos, auto &max_allowed_occ,
                auto &ambiguous_hit_indices) -> bool {
      (void)verbose;
      int32_t hit_idx{0};
      hit_map.clear();

      for (auto &raw_hit : raw_hits) {
        // std::cout << "hit\n";
        auto &read_pos = raw_hit.first;
        auto &proj_hits = raw_hit.second;
        auto &refs = proj_hits.refRange;
        // std::cout << "refs len, read_pos " << refs.size() << " " << read_seq->substr(read_pos, 23) << std::endl;
        uint64_t num_occ = static_cast<uint64_t>(refs.size());
        min_occ = std::min(min_occ, num_occ);

        prev_read_pos = read_pos;
        if (num_occ <= max_allowed_occ) {
          total_occs += num_occ;
          largest_occ = (num_occ > largest_occ) ? num_occ : largest_occ;
          float score_inc = 1.0;

          for (auto v : refs) {
            const auto &ref_pos_ori = proj_hits.decode_hit(v);
            // std::cout << "ref pos " << ref_pos_ori.pos << std::endl;
            uint32_t tid = sshash::util::transcript_id(v);
            int32_t pos = static_cast<int32_t>(ref_pos_ori.pos);
            int32_t signed_read_pos = static_cast<int32_t>(read_pos);
            std::pair<uint64_t, uint64_t> bins = binning.get_bin_id(tid, pos);
            // std::cout << "tid bins " << tid << " " << bins.first << " " << bins.second << std::endl;
            bool ori = ref_pos_ori.isFW;

            auto& target1 = hit_map[bins.first];
            target1.tid = tid;
            if (ori) {
              target1.add_fw(pos, signed_read_pos, signed_rl, k, max_stretch, score_inc);
            } else {
              target1.add_rc(pos, signed_read_pos, signed_rl, k, max_stretch, score_inc);
            }
 
            typename std::remove_reference<decltype(target1)>::type* target2 = nullptr;
            if (binning.is_valid(bins.second)) {
              target2 = &hit_map[bins.second];
              target2->tid = tid;
              if (ori) {
                target2->add_fw(pos, signed_read_pos, signed_rl, k, max_stretch, score_inc);
              } else {
                target2->add_rc(pos, signed_read_pos, signed_rl, k, max_stretch, score_inc);
              }
            }
          } // DONE: for (auto &pos_it : refs)
          ++num_valid_hits;
        } else if (perform_ambig_filtering) { // HERE we have that num_occ >
                                              // max_allowed_occ
            ambiguous_hit_indices.push_back(hit_idx);
        }
        ++hit_idx;
      } // DONE : for (auto& raw_hit : raw_hits)

      return false;
    };

    auto map_first_pass = map_cache.max_hit_occ - 1;
    early_stop = collect_mappings_from_hits_thr(
      raw_hits, prev_read_pos, map_first_pass, map_cache.ambiguous_hit_indices);
    

    // If our default threshold was too stringent, then fallback to a more
    // liberal threshold and look up the k-mers that occur the least frequently.
    // Specifically, if the min occuring hits have frequency <
    // max_hit_occ_recover (2500 by default) times, then collect the min
    // occuring hits to get the mapping.
    if (attempt_occ_recover and (min_occ >= map_cache.max_hit_occ) and
        (min_occ < map_cache.max_hit_occ_recover)) {
      num_valid_hits = 0;
      map_cache.ambiguous_hit_indices.clear();
      prev_read_pos = -1;
      uint64_t max_allowed_occ = min_occ;
      early_stop =
        collect_mappings_from_hits_thr(raw_hits, prev_read_pos, max_allowed_occ,
                                   map_cache.ambiguous_hit_indices);
    }
    // Further filtering of mappings by ambiguous k-mers
    if (perform_ambig_filtering and !hit_map.empty() and
        !map_cache.ambiguous_hit_indices.empty()) {
      phmap::flat_hash_set<uint64_t> observed_ecs;
      size_t min_cardinality_ec_size = std::numeric_limits<size_t>::max();
      uint64_t min_cardinality_ec = std::numeric_limits<size_t>::max();
      size_t min_cardinality_index = 0;
      size_t visited = 0;
      auto &ec_table = map_cache.hs.get_index()->get_ec_table();

      auto visit_ec = [&hit_map](uint64_t ent, bool fw_on_contig) -> bool {
        uint32_t bin_id = (ent >> 2);
        auto hm_it = hit_map.find(bin_id);
        bool found = false;
        if (hm_it != hit_map.end()) {
          // we found this target, now:
          // (1) check the orientation
          uint32_t ori = (ent & 0x3);
          // (2) add hits in the appropriate way
          switch (ori) {
          case 0: // fw
            (fw_on_contig) ? hm_it->second.inc_fw_hits()
                           : hm_it->second.inc_rc_hits();
            break;
          case 1: // rc
            (fw_on_contig) ? hm_it->second.inc_rc_hits()
                           : hm_it->second.inc_fw_hits();
            break;
          default: // both
            hm_it->second.inc_fw_hits();
            hm_it->second.inc_rc_hits();
          }
          found = true;
        }
        return found;
      };

      // for each ambiguous hit
      for (auto hit_idx : map_cache.ambiguous_hit_indices) {
        auto &proj_hit = raw_hits[hit_idx].second;
        uint32_t contig_id = proj_hit.contig_id();
        bool fw_on_contig = proj_hit.hit_fw_on_contig();

        // put the combination of the eq and the k-mer orientation
        // into the map.
        uint64_t ec = ec_table.ec_for_tile(contig_id);
        uint64_t ec_key = ec | (fw_on_contig ? 0 : 0x8000000000000000);
        // if we've already seen this ec, no point in processing
        // it again.
        if (observed_ecs.contains(ec_key)) {
          continue;
        }
        // otherwise, insert it.
        observed_ecs.insert(ec_key);

        auto ec_entries = ec_table.entries_for_ec(ec);
        if (ec_entries.size() < min_cardinality_ec_size) {
          min_cardinality_ec_size = ec_entries.size();
          min_cardinality_ec = ec;
          min_cardinality_index = hit_idx;
        }
        if (ec_entries.size() > max_ec_ambig) {
          continue;
        }
        ++visited;
        for (const auto ent : ec_entries) {
          visit_ec(ent, fw_on_contig);
        } // all target oritentation pairs in this eq class
        ++num_valid_hits;
      } // all ambiguous hits

      // if we haven't visited *any* equivalence classes (they were all)
      // too ambiguous, then make a last-ditch effort to just visit the
      // the one with smallest cardinality.
      if (visited == 0) {
        auto hit_idx = min_cardinality_index;
        auto &proj_hit = raw_hits[hit_idx].second;
        bool fw_on_contig = proj_hit.hit_fw_on_contig();

        uint64_t ec = min_cardinality_ec;
        auto ec_entries = ec_table.entries_for_ec(ec);
        for (const auto ent : ec_entries) {
          visit_ec(ent, fw_on_contig);
        } // all target oritentation pairs in this eq class
        ++num_valid_hits;
      } // done visiting the last-ditch ec

    } // if we are processing ambiguous hits

    // we need at least threshold * the maximum number of hits
    uint32_t best_alt_hits = 0;
    num_valid_hits = std::max(static_cast<uint32_t>(1), static_cast<uint32_t>(std::ceil(num_valid_hits*thr))); 
    // std::cout << "hm size " << hit_map.size() << std::endl;  
    // std::cout << "n valid " << num_valid_hits << std::endl; 
    for (auto &kv : hit_map) {
      // if (hit_map.size() <= 10) {
      //   std::cout << "tid " << kv.first
      // }
      auto best_hit_dir = kv.second.best_hit_direction();

      // if the best direction is FW or BOTH, add the fw hit
      // otherwise add the RC.
      auto simple_hit = (best_hit_dir != mapping::util::HitDirection::RC)
                          ? kv.second.get_fw_hit()
                          : kv.second.get_rc_hit();
      if (simple_hit.num_hits >= num_valid_hits) {
        // std::cout << "yes\n";
        simple_hit.bin_id = kv.first;
        simple_hit.tid = kv.second.tid;
        accepted_hits.emplace_back(simple_hit);
        // if we had equally good hits in both directions
        // add the rc hit here (since we added the fw)
        // above if the best hit was either FW or BOTH
        if (best_hit_dir == mapping::util::HitDirection::BOTH) {
          auto second_hit = kv.second.get_rc_hit();
          second_hit.bin_id = kv.first;
          second_hit.tid = kv.second.tid;
          accepted_hits.emplace_back(second_hit);
        }
      } else {
          best_alt_hits = simple_hit.num_hits > best_alt_hits
                          ? simple_hit.num_hits
                          : best_alt_hits;
      }
    }

    // max_read_occ = had_max_read_occ ? accepted_hits.size() : max_hit_occ;
  } // DONE : if (rh)
  // If the read mapped to > maxReadOccs places, discard it
  if (accepted_hits.size() > map_cache.max_read_occ) {
    // std::cout << "multi\n";
    accepted_hits.clear();
    map_type = mapping::util::MappingType::UNMAPPED;
  } else if (!accepted_hits.empty()) {
    map_type = mapping::util::MappingType::SINGLE_MAPPED;
    // std::cout << "yp\n";
  }
  // std::cout << "es " << early_stop << std::endl;
  // std::cout << "accepted_hits " << accepted_hits.size() << std::endl;
  // if (accepted_hits.size() < 10) {
  //   print_hits(accepted_hits);
  // }
  // std::cout << "checking\n";
  return early_stop;
}

template <typename mapping_cache_info_t>
inline void merge_se_mappings(mapping_cache_info_t &map_cache_left,
                              mapping_cache_info_t &map_cache_right,
                              int32_t left_len, int32_t right_len,
                              mapping_cache_info_t &map_cache_out) {
  map_cache_out.clear();
  auto &accepted_left = map_cache_left.accepted_hits;
  auto &accepted_right = map_cache_right.accepted_hits;

  size_t had_matching_kmers_left = map_cache_left.has_matching_kmers;
  size_t had_matching_kmers_right = map_cache_right.has_matching_kmers;

  size_t num_accepted_left = accepted_left.size();
  size_t num_accepted_right = accepted_right.size();

  if ((num_accepted_left > 0) and (num_accepted_right > 0)) {
    // look for paired end mappings
    // so we have to sort our accepted hits
    struct {
      // sort first by orientation, then by transcript id, and finally by
      // position
      bool operator()(const mapping::util::simple_hit &a,
                      const mapping::util::simple_hit &b) {
        if (a.is_fw != b.is_fw) {
          return a.is_fw > b.is_fw;
        }
        // orientations are the same
        if (a.tid != b.tid) {
          return a.tid < b.tid;
        }
        return a.pos < b.pos;
      }
    } simple_hit_less;
    std::sort(accepted_left.begin(), accepted_left.end(), simple_hit_less);
    std::sort(accepted_right.begin(), accepted_right.end(), simple_hit_less);

    const mapping::util::simple_hit smallest_rc_hit = {false, false, -1,
                                                       0.0,   0,     0};
    // start of forward sub-list
    auto first_fw1 = accepted_left.begin();
    // end of forward sub-list is first non-forward hit
    auto last_fw1 = std::lower_bound(accepted_left.begin(), accepted_left.end(),
                                     smallest_rc_hit, simple_hit_less);
    // start of rc list
    auto first_rc1 = last_fw1;
    // end of rc list
    auto last_rc1 = accepted_left.end();

    // start of forward sub-list
    auto first_fw2 = accepted_right.begin();
    // end of forward sub-list is first non-forward hit
    auto last_fw2 =
      std::lower_bound(accepted_right.begin(), accepted_right.end(),
                       smallest_rc_hit, simple_hit_less);
    // start of rc list
    auto first_rc2 = last_fw2;
    // end of rc list
    auto last_rc2 = accepted_right.end();

    auto back_inserter = std::back_inserter(map_cache_out.accepted_hits);
    using iter_t = decltype(first_fw1);
    using out_iter_t = decltype(back_inserter);

    auto merge_lists = [left_len, right_len](iter_t first1, iter_t last1,
                                             iter_t first2, iter_t last2,
                                             out_iter_t out) -> out_iter_t {
      // https://en.cppreference.com/w/cpp/algorithm/set_intersection
      while (first1 != last1 && first2 != last2) {
        if (first1->tid < first2->tid) {
          ++first1;
        } else {
          if (!(first2->tid < first1->tid)) {
            // first1->tid == first2->tid have the same transcript.
            int32_t pos_fw = first1->is_fw ? first1->pos : first2->pos;
            int32_t pos_rc = first1->is_fw ? first2->pos : first1->pos;
            int32_t frag_len = (pos_rc - pos_fw);
            if ((-32 < frag_len) and (frag_len < 2000)) {
              // if left is fw and right is rc then
              // fragment length is (right_pos + right_len - left_pos) + 1
              // otherwise it is (left_pos + left_len - right_pos) + 1
              bool right_is_rc = !first2->is_fw;
              int32_t tlen = right_is_rc
                               ? ((first2->pos + right_len - first1->pos) + 1)
                               : ((first1->pos + left_len - first2->pos) + 1);
              *out++ = {first1->is_fw, first2->is_fw, first1->pos, 0.0, 0,
                        first1->tid,   first2->pos,   tlen};
              ++first1;
            }
          }
          ++first2;
        }
      }
      return out;
    };

    // find hits of form 1:fw, 2:rc
    merge_lists(first_fw1, last_fw1, first_rc2, last_rc2, back_inserter);
    // find hits of form 1:rc, 2:fw
    merge_lists(first_rc1, last_rc1, first_fw2, last_fw2, back_inserter);

    map_cache_out.map_type = (map_cache_out.accepted_hits.size() > 0)
                               ? MappingType::MAPPED_PAIR
                               : MappingType::UNMAPPED;
  } else if ((num_accepted_left > 0) and !had_matching_kmers_right) {
    // just return the left mappings
    std::swap(map_cache_left.accepted_hits, map_cache_out.accepted_hits);
    map_cache_out.map_type = (map_cache_out.accepted_hits.size() > 0)
                               ? MappingType::MAPPED_FIRST_ORPHAN
                               : MappingType::UNMAPPED;
  } else if ((num_accepted_right > 0) and !had_matching_kmers_left) {
    // just return the right mappings
    std::swap(map_cache_right.accepted_hits, map_cache_out.accepted_hits);
    map_cache_out.map_type = (map_cache_out.accepted_hits.size() > 0)
                               ? MappingType::MAPPED_SECOND_ORPHAN
                               : MappingType::UNMAPPED;
  } else {
    // return nothing
  }
}

} // namespace util

namespace util_bin {

template <typename mapping_cache_info_t>
inline void remove_duplicate_hits(mapping_cache_info_t &map_cache, uint32_t max_num_hits) {
  auto& accepted_hits = map_cache.accepted_hits;
  using hit_list_iter_t = decltype(accepted_hits.begin());

  struct {
    // hits are equivalent if they have the same orientation, reference, and position
    bool operator()(const mapping::util::simple_hit& a,
                    const mapping::util::simple_hit& b) {
      return ((a.is_fw == b.is_fw) and (a.tid == b.tid) and (a.pos == b.pos));
    }
  } equiv_hit;

  
  auto canonicalize_hits = [&equiv_hit](hit_list_iter_t it1, hit_list_iter_t it2, hit_list_iter_t end, uint32_t &max_nh) {
    while (it2 != end) {
      if (equiv_hit(*it1, *it2)) { 
        auto min_bin = std::min(it1->bin_id, it2->bin_id);
        auto max_num_hits = std::max(it1->num_hits, it2->num_hits);
        it1->bin_id = it2->bin_id = min_bin;
        it1->num_hits = it2->num_hits = max_num_hits;
        max_nh = std::max(max_nh, max_num_hits);
      }
      ++it1; ++it2;
    }
  };
  
  
  uint32_t m_hits = max_num_hits == 0 ? 0 : accepted_hits.front().num_hits;
  
  canonicalize_hits(accepted_hits.begin(), accepted_hits.begin() + 1, 
  accepted_hits.end(), m_hits);  
  auto last = std::unique(accepted_hits.begin(), accepted_hits.end(), equiv_hit);
  accepted_hits.erase(last, accepted_hits.end());
  

    
}

template <typename mapping_cache_info_t>
inline void merge_se_mappings(mapping_cache_info_t& map_cache_left,
                              mapping_cache_info_t& map_cache_right, int32_t left_len,
                              int32_t right_len, mapping_cache_info_t& map_cache_out 
                              ) {
  map_cache_out.clear();
  
  auto& accepted_left = map_cache_left.accepted_hits;
  auto& accepted_right = map_cache_right.accepted_hits;
  
  size_t had_matching_kmers_left = map_cache_left.has_matching_kmers;
  size_t had_matching_kmers_right = map_cache_right.has_matching_kmers;

  size_t num_accepted_left = accepted_left.size();
  size_t num_accepted_right = accepted_right.size();

  uint32_t max_num_hits = 0;
  struct {
    // sort first by orientation, then by transcript id, position, num_hits and, 
    // finaly bin_id.
    bool operator()(const mapping::util::simple_hit& a,
                    const mapping::util::simple_hit& b) {
      if (a.is_fw != b.is_fw) { return a.is_fw > b.is_fw; }
      // orientations are the same
      if (a.tid != b.tid) { return a.tid < b.tid; }
      // orientations & txps are the same
      if (a.pos != b.pos) { return a.pos < b.pos; }
      // orientations, txps, and positions are the same
      if (a.num_hits != b.num_hits ) { return a.num_hits < b.num_hits ; }
      // orientations, txps, positions and num_hits are the same 
      return a.bin_id < b.bin_id;
    }
  } simple_hit_less_bins;

/*  struct {
    // hits are equivalent if they have the same orientation, reference, and position
    bool operator()(const mapping::util::simple_hit& a,
                    const mapping::util::simple_hit& b) {
      return ((a.is_fw == b.is_fw) and (a.tid == b.tid) and (a.pos == b.pos));
    }
  } equiv_hit;*/

  // scan through and "canonicalize" duplicate hits.  In this 
  // process, non-duplcate hits are left untouched, and for duplicates
  // (defined by the captured `equiv_hit`) functor, each duplicate gets 
  // the smaller bin id and the larger score.
  // using hit_list_iter_t = decltype(accepted_left.begin());
  // arguments are a pair of adjacent iterators to the first and second elements (`it1` and `it2`)
  // and an `end` iterator (following the standard meaning, it is one past the last valid element).
  // The `max_nh` argument will record the max num_hits observed over all mappings and 
  // be updated accordingly.
/*  auto canonicalize_hits = [&equiv_hit](hit_list_iter_t it1, hit_list_iter_t it2, hit_list_iter_t end, uint32_t& max_nh) {
    while (it2 != end) {
      if (equiv_hit(*it1, *it2)) { 
        auto min_bin = std::min(it1->bin_id, it2->bin_id);
        auto max_num_hits = std::max(it1->num_hits, it2->num_hits);
        it1->bin_id = it2->bin_id = min_bin;
        it1->num_hits = it2->num_hits = max_num_hits;
        max_nh = std::max(max_nh, max_num_hits);
      }
      ++it1; ++it2;
    }
  };*/

  
  if ((num_accepted_left > 0) and (num_accepted_right > 0)) {
    // sort by orientation, tid, position, num_hits and finally bin_id
    std::sort(accepted_left.begin(), accepted_left.end(), simple_hit_less_bins);
    std::sort(accepted_right.begin(), accepted_right.end(), simple_hit_less_bins);

    // since these are paired-end mappings
    // we don't want to record individual end num_hits 
    // here so we pass a dummy variable.
    uint32_t discard = 0;
    // canonicalize the hits, giving each equivalent hit the smaller 
    // of the possible bin ids and the larger of the possible scores.
    // canonicalize_hits(accepted_left.begin(), accepted_left.begin() + 1, accepted_left.end(), discard);
    // canonicalize_hits(accepted_right.begin(), accepted_right.begin() + 1, accepted_right.end(), discard);
    remove_duplicate_hits(map_cache_left, discard);
    remove_duplicate_hits(map_cache_right, discard);
    // remove duplicates
    // auto last_left = std::unique(accepted_left.begin(), accepted_left.end(), equiv_hit);
    // auto last_right = std::unique(accepted_right.begin(), accepted_right.end(), equiv_hit);
    // accepted_left.erase(last_left, accepted_left.end());
    // accepted_right.erase(last_right, accepted_right.end());
  
    // std::cout << "left\n";
    // mapping::util::print_hits(accepted_left);
    // std::cout << "right\n";
    // mapping::util::print_hits(accepted_right);

    const mapping::util::simple_hit smallest_rc_hit = {false, false, -1, 0.0, 0, 0, 0, 0, 0};
    // start of forward sub-list

  
    auto first_fw1 = accepted_left.begin();
    // std::cout << "first " << first_fw1->bin_id << std::endl;
    // end of forward sub-list is first non-forward hit
    auto last_fw1 = std::lower_bound(accepted_left.begin(), accepted_left.end(),
                                     smallest_rc_hit, simple_hit_less_bins);
    // bool a = first_fw1 == last_fw1;
    // std::cout << a << std::endl;
    // std::cout << "last fw " << last_fw1->bin_id << std::endl;
    // start of rc list
    auto first_rc1 = last_fw1;
    // std::cout << "first rc " << first_rc1->bin_id << std::endl;
    // end of rc list
    auto last_rc1 = accepted_left.end();
    // std::cout << "last rc " << last_rc1->bin_id << std::endl;
    // start of forward sub-list
    auto first_fw2 = accepted_right.begin();
    // end of forward sub-list is first non-forward hit
    auto last_fw2 = std::lower_bound(accepted_right.begin(), accepted_right.end(),
                                     smallest_rc_hit, simple_hit_less_bins);
    // start of rc list
    auto first_rc2 = last_fw2;
    // end of rc list
    auto last_rc2 = accepted_right.end();

    auto back_inserter = std::back_inserter(map_cache_out.accepted_hits);
    using iter_t = decltype(first_fw1);
    using out_iter_t = decltype(back_inserter);
    auto merge_lists = [left_len, right_len, &max_num_hits](iter_t first1, iter_t last1, iter_t first2,
                                                            iter_t last2, out_iter_t out) -> out_iter_t {
        // https://en.cppreference.com/w/cpp/algorithm/set_intersection
        while (first1 != last1) {
          // we have to consider mappings that span a bin-bin boundary 
          // because the left read may reside in the first bin and the 
          // right read in the second bin
          // std::cout << "first1 " << first1->bin_id << std::endl;
          auto f2 = first2;
          while (f2 != last2) {
            // std::cout << "f2 " << f2->bin_id << " " << f2->pos << std::endl;
            // if (f2->bin_id - first1->bin_id > 1) {
            //   ++first1;
            //   break;
            // }
            // std::cout << f2->is_fw << " " << f2->bin_id << std::endl;
            if ((f2->tid == first1->tid) && (f2->is_fw != first1->is_fw)) {
              // std::cout << "dd\n";
              if ((f2->bin_id + 1 == first1->bin_id) || (f2->bin_id == first1->bin_id + 1) ||
                (f2->bin_id == first1->bin_id)) {
            
                int32_t pos_fw = first1->is_fw ? first1->pos : f2->pos;
                int32_t pos_rc = first1->is_fw ? f2->pos : first1->pos;
                int32_t frag_len = (pos_rc - pos_fw);
                // std::cout << "frag len " << frag_len << std::endl;
                if ((-20 < frag_len) and (frag_len < 1000)) {
                  // std::cout << "yes\n";
                  // if left is fw and right is rc then
                  // fragment length is (right_pos + right_len - left_pos) + 1
                  // otherwise it is (left_pos + left_len - right_pos) + 1
                  bool right_is_rc = !first2->is_fw;
                  int32_t tlen = right_is_rc
                    ? ((f2->pos + right_len - first1->pos) + 1)
                    : ((first1->pos + left_len - f2->pos) + 1);

                  uint32_t nhits = first1->num_hits + f2->num_hits;
                  max_num_hits = std::max(max_num_hits, nhits);
                  *out++ = {first1->is_fw, f2->is_fw, first1->pos, 0.0, nhits,
                    first1->tid, f2->pos, tlen, first1->bin_id};
                }
              }
            }
            ++f2;
          }
          ++first1;
        }
        return out;
      };

    // find hits of form 1:fw, 2:rc
    merge_lists(first_fw1, last_fw1, first_rc2, last_rc2, back_inserter);
    // find hits of form 1:rc, 2:fw
    merge_lists(first_rc1, last_rc1, first_fw2, last_fw2, back_inserter);

    map_cache_out.map_type = (map_cache_out.accepted_hits.size() > 0) ? util::MappingType::MAPPED_PAIR
      : util::MappingType::UNMAPPED;
  // } else if ((num_accepted_left > 0) and !had_matching_kmers_right) {
  } else if (num_accepted_left > 0) {
    // sort by orientation, tid, position, num_hits and finally bin_id
    std::sort(accepted_left.begin(), accepted_left.end(), simple_hit_less_bins);
    max_num_hits = accepted_left.front().num_hits;
    remove_duplicate_hits(map_cache_left, max_num_hits);
    // auto last_left = std::unique(accepted_left.begin(), accepted_left.end(), equiv_hit);
    // remove duplicates
    // accepted_left.erase(last_left, accepted_left.end());

    std::swap(map_cache_left.accepted_hits, map_cache_out.accepted_hits);
    map_cache_out.map_type = (map_cache_out.accepted_hits.size() > 0)
      ? util::MappingType::MAPPED_FIRST_ORPHAN
      : util::MappingType::UNMAPPED;
  // } else if ((num_accepted_right > 0) and !had_matching_kmers_left) {
   } else if ((num_accepted_right > 0)) {
    // sort by orientation, tid, position, num_hits and finally bin_id
    std::sort(accepted_right.begin(), accepted_right.end(), simple_hit_less_bins);
    // max_num_hits = accepted_right.front().num_hits;
    // canonicalize_hits(accepted_right.begin(), accepted_right.begin() + 1, accepted_right.end(), max_num_hits);
    // auto last_right = std::unique(accepted_right.begin(), accepted_right.end(), equiv_hit);
    // remove duplicates
    // accepted_right.erase(last_right, accepted_right.end());
    max_num_hits = accepted_right.front().num_hits;
    remove_duplicate_hits(map_cache_right, max_num_hits);
    std::swap(map_cache_right.accepted_hits, map_cache_out.accepted_hits);
    map_cache_out.map_type = (map_cache_out.accepted_hits.size() > 0)
      ? util::MappingType::MAPPED_SECOND_ORPHAN
      : util::MappingType::UNMAPPED;
  } else {
    //  if (((num_accepted_right > 0) and had_matching_kmers_left) || 
    //     (((num_accepted_left > 0) and had_matching_kmers_right))) {
    //     std::cout << "weird case\n";
    //     }
    // return nothing
  }

  // While even if we are thresholding, only top hits should be kept
  if (map_cache_out.accepted_hits.size() > 0) {
    auto accepted_hits_last = 
      std::remove_if(map_cache_out.accepted_hits.begin(), 
                     map_cache_out.accepted_hits.end(),
                     [max_num_hits](const mapping::util::simple_hit& h) -> bool {
                     return h.num_hits < max_num_hits;
                     });
    map_cache_out.accepted_hits.erase(accepted_hits_last, map_cache_out.accepted_hits.end());
  }
}

}
} // namespace mapping
