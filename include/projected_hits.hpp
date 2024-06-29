#pragma once

#include "util.hpp"
#include <iostream>
struct ref_pos {
    uint32_t pos;
    bool isFW;
};

struct projected_hits {
    uint32_t contigIdx_;
    // The relative position of the k-mer inducing this hit on the
    // contig
    uint32_t contigPos_;
    // How the k-mer inducing this hit maps to the contig
    // true for fw, false for rc
    bool contigOrientation_;
    uint32_t contigLen_;
    uint64_t globalPos_;
    uint32_t k_;

    sshash::util::contig_span refRange;
    
    bool resulted_from_open_search{false};

    inline bool empty() { return refRange.empty(); }

    inline uint32_t contig_id() const { return contigIdx_; }
    inline uint32_t contig_pos() const { return contigPos_; }
    inline uint32_t contig_len() const { return contigLen_; }
    inline bool hit_fw_on_contig() const { return contigOrientation_; }

    inline ref_pos decode_hit(uint64_t v) {
        // true if the contig is fowrard on the reference
        bool contigFW = sshash::util::orientation(v);
        // we are forward with respect to the reference if :
        // (1) contigFW and contigOrientation_
        // (2) !contigFW and !contigOrientation_
        // we are reverse complement with respect to the reference if :
        // (3) configFW and !contigOrientation_
        // (4) !configFW and contigOrientation_

        // if we're in the forward orientation, then our position is
        // just the contig offset plus or relative position
        uint32_t rpos;  //{0};
        bool rfw;       //{false};
        if (contigFW and contigOrientation_) {
            // kmer   :          AGC
            // contig :      ACTTAGC
            // ref    :  GCA[ACTTAGC]CA
            rpos = sshash::util::pos(v) + contigPos_;
            rfw = true;
        } else if (contigFW and !contigOrientation_) {
            // kmer   :          GCT
            // contig :      ACTTAGC
            // ref    :  GCA[ACTTAGC]CA
            rpos = sshash::util::pos(v) + contigPos_;
            rfw = false;
        } else if (!contigFW and contigOrientation_) {
            // kmer   :          AGT
            // contig :      GCTAAGT
            // ref    :  GCA[ACTTAGC]CA
            rpos = sshash::util::pos(v) + contigLen_ - (contigPos_ + k_);
            rfw = false;
        } else {  // if (!contigFW and !contigOrientation_) {
            // kmer   :          ACT
            // contig :      GCTAAGT
            // ref    :  GCA[ACTTAGC]CA
            rpos = sshash::util::pos(v) + contigLen_ - (contigPos_ + k_);
            rfw = true;
        }

        return {rpos, rfw};
    }

    // inline friend function :
    // https://stackoverflow.com/questions/381164/friend-and-inline-method-whats-the-point 
    // this helps to avoid duplicate symbol error.
    inline friend std::ostream& operator<<(std::ostream& os, projected_hits& h) {
        os << "{ proj_hit : "
           << "\t{ contig_idx : " << h.contigIdx_ << ", "
           << "contig_pos : " << h.contigPos_ << ", "
           << "contig_ori : " << (h.contigOrientation_ ? "fw" : "rc") << ", "
           << "contig_len : " << h.contigLen_ << ", "
           << "global_pos : " << h.globalPos_ << ", "
           << "ref_range_len : " << h.refRange.size() << "}\t}";
        return os;
    }

  inline bool operator!=(const projected_hits& p2) {
    return (contigIdx_ != p2.contigIdx_) or
           (contigPos_ != p2.contigPos_) or
           (contigOrientation_ != p2.contigOrientation_) or
           (contigLen_ != p2.contigLen_) or
           (globalPos_ != p2.globalPos_) or
           (k_ != p2.k_);
  }

  inline bool operator==(const projected_hits& p2) {
    return !operator!=(p2);
  }
};


