#pragma once

#include "util.hpp"

namespace sshash {

struct skew_index {
    skew_index()
        : min_log2(constants::min_l)
        , max_log2(constants::max_l)
        , max_num_strings_in_bucket(0)
        , log2_max_num_strings_in_bucket(0) {
        mphfs.resize(0);
        positions.resize(0);
    }

    /* Returns the number of k-mers in the index. */
    uint64_t print_info() const;

    bool empty() const { return mphfs.empty(); }

    uint64_t lookup(uint64_t uint64_kmer, uint64_t log2_num_strings_in_bucket) const {
        assert(log2_num_strings_in_bucket >= uint64_t(min_log2 + 1));
        assert(log2_num_strings_in_bucket <= log2_max_num_strings_in_bucket);
        uint64_t partition_id = log2_num_strings_in_bucket - (min_log2 + 1);
        if (log2_num_strings_in_bucket == log2_max_num_strings_in_bucket or
            log2_num_strings_in_bucket > max_log2) {
            partition_id = positions.size() - 1;
        }
        auto const& mphf = mphfs[partition_id];
        auto const& P = positions[partition_id];
        uint64_t position = P.access(mphf(uint64_kmer));
        return position;
    }

    uint64_t num_bits() const {
        uint64_t n = (sizeof(min_log2) + sizeof(max_log2) + sizeof(max_num_strings_in_bucket) +
                      sizeof(log2_max_num_strings_in_bucket)) *
                     8;
        for (uint64_t partition_id = 0; partition_id != mphfs.size(); ++partition_id) {
            auto const& mphf = mphfs[partition_id];
            auto const& P = positions[partition_id];
            n += mphf.num_bits() + P.bytes() * 8;
        }
        return n;
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visitor.visit(min_log2);
        visitor.visit(max_log2);
        visitor.visit(max_num_strings_in_bucket);
        visitor.visit(log2_max_num_strings_in_bucket);
        visitor.visit(mphfs);
        visitor.visit(positions);
    }

    uint16_t min_log2;
    uint16_t max_log2;
    uint32_t max_num_strings_in_bucket;  // useless
    uint32_t log2_max_num_strings_in_bucket;
    std::vector<pthash_mphf_type> mphfs;
    std::vector<pthash::compact_vector> positions;
};

}  // namespace sshash