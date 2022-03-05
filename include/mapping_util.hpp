#pragma once

#include <vector>
#include <cassert>
#include <fstream>
#include <cmath>  // for std::ceil on linux
#include <numeric>


namespace mapping {

namespace util {

struct simple_hit {
    bool is_fw{false};
    int32_t pos{-1};
    float score{0.0};
    uint32_t num_hits{0};
    uint32_t tid{std::numeric_limits<uint32_t>::max()};
    bool valid_pos(int32_t read_len, uint32_t txp_len, int32_t max_over) {
      int32_t signed_txp_len = static_cast<int32_t>(txp_len);
      return (pos > -max_over) and ((pos + read_len) < (signed_txp_len + max_over));
    }
  };

  enum class HitDirection : uint8_t {FW, RC, BOTH};

  struct sketch_hit_info {

    // add a hit to the current target that occurs in the forward 
    // orientation with respect to the target.
    bool add_fw(int32_t ref_pos, int32_t read_pos, int32_t rl, int32_t k, int32_t max_stretch, float score_inc) {
      (void)rl;
      (void)k;
      bool added{false};
      
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
        //if (last_ref_pos_fw > -1 and (ref_pos > last_ref_pos_fw + 15)) { return false; }
        last_ref_pos_fw = ref_pos;
        last_read_pos_fw = read_pos;
        fw_score += score_inc;
        ++fw_hits;
        added = true;
      }   
      return added;
    }

    // add a hit to the current target that occurs in the forward 
    // orientation with respect to the target.
    bool add_rc(int32_t ref_pos, int32_t read_pos, int32_t rl, int32_t k, int32_t max_stretch, float score_inc) {

      bool added{false};
      // since hits are collected by moving _forward_ in the
      // read, if this is an rc hit, it should be moving 
      // backwards in the reference. Only add it if this is
      // the case.
      // This ensures that we don't double-count a k-mer that 
      // might occur twice on this target.
      if (ref_pos < last_ref_pos_rc and read_pos > last_read_pos_rc) {
        approx_pos_rc = (ref_pos - (rl - (read_pos + k)));
        if (last_read_pos_rc == -1) { 
          approx_end_pos_rc = ref_pos + read_pos; 
        } else {
          if (approx_end_pos_rc - approx_pos_rc > max_stretch) { return false;}
        }
        //if (last_ref_pos_rc > -1 and ref_pos < last_ref_pos_rc - 15) { return false; }
        last_ref_pos_rc = ref_pos;
        last_read_pos_rc = read_pos;
        rc_score += score_inc;
        ++rc_hits;
        added = true;
        
      }
      return added;
    }

    inline uint32_t max_hits_for_target() {
      return std::max(fw_hits, rc_hits);
    }

    // true if forward, false if rc
    // second element is score
    inline HitDirection best_hit_direction() {
      int32_t fw_minus_rc = static_cast<int32_t>(fw_hits) - static_cast<int32_t>(rc_hits);
      return (fw_minus_rc > 0) ? HitDirection::FW : 
             ((fw_minus_rc < 0) ? HitDirection::RC : HitDirection::BOTH);
    }

    inline simple_hit get_fw_hit() {
      return simple_hit{true, approx_pos_fw, fw_score, fw_hits, std::numeric_limits<uint32_t>::max()}; 
    }

    inline simple_hit get_rc_hit() {
      return simple_hit{false, approx_pos_rc, rc_score, rc_hits, std::numeric_limits<uint32_t>::max()};
    }

    inline simple_hit get_best_hit() {
      auto best_direction = best_hit_direction();
      return (best_direction != HitDirection::RC) ? 
        simple_hit{true, approx_pos_fw, fw_score, fw_hits, std::numeric_limits<uint32_t>::max()} : 
        simple_hit{false, approx_pos_rc, rc_score, rc_hits, std::numeric_limits<uint32_t>::max()};
    }

    inline std::string to_string() {
      std::stringstream ss;
      ss << "fw_hits: " << fw_hits << ", fw_score : " << fw_score << ", fw_pos : " << approx_pos_fw
         << " || rc_hits: " << rc_hits << ", rc_score: " << rc_score << ", rc_pos: " << approx_pos_rc;
      return ss.str();
    }

    int32_t last_read_pos_fw{-1};
    int32_t last_read_pos_rc{-1};

    int32_t last_ref_pos_fw{-1};
    int32_t last_ref_pos_rc{std::numeric_limits<int32_t>::max()};
    
    int32_t approx_pos_fw{-1};
    int32_t approx_pos_rc{-1};
    int32_t approx_end_pos_rc{-1};

    uint32_t fw_hits{0};
    uint32_t rc_hits{0};
    float fw_score{0.0};
    float rc_score{0.0}; 
  };

}

}