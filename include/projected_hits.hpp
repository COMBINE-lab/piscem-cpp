#pragma once

#include "nonstd/span.hpp"

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
    uint32_t k_;
    nonstd::span<sshash::util::Position> refRange;

    inline bool empty() { return refRange.empty(); }

    inline uint32_t contig_id() const { return contigIdx_; }

    inline ref_pos decode_hit(sshash::util::Position& p) {
        // true if the contig is fowrard on the reference
        bool contigFW = p.orientation();
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
            rpos = p.pos() + contigPos_;
            rfw = true;
        } else if (contigFW and !contigOrientation_) {
            // kmer   :          GCT
            // contig :      ACTTAGC
            // ref    :  GCA[ACTTAGC]CA
            rpos = p.pos() + contigPos_;
            rfw = false;
        } else if (!contigFW and contigOrientation_) {
            // kmer   :          AGT
            // contig :      GCTAAGT
            // ref    :  GCA[ACTTAGC]CA
            rpos = p.pos() + contigLen_ - (contigPos_ + k_);
            rfw = false;
        } else {  // if (!contigFW and !contigOrientation_) {
            // kmer   :          ACT
            // contig :      GCTAAGT
            // ref    :  GCA[ACTTAGC]CA
            rpos = p.pos() + contigLen_ - (contigPos_ + k_);
            rfw = true;
        }

        return {rpos, rfw};
    }
};