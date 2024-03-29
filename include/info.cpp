#include "dictionary.hpp"
#include "skew_index.hpp"
#include "spdlog_piscem/spdlog.h"

namespace sshash {

uint64_t skew_index::print_info() const {
    uint64_t num_partitions = mphfs.size();
    uint64_t lower = 1ULL << min_log2;
    uint64_t upper = 2 * lower;
    uint64_t num_kmers_in_skew_index = 0;
    for (uint64_t partition_id = 0; partition_id != num_partitions; ++partition_id) {
        uint64_t n = mphfs[partition_id].num_keys();
        assert(n == positions[partition_id].size());
        spdlog_piscem::info(
            "num_kmers belonging to buckets of size > {} and <= {}: {}; bits/kmer = {} (mphf) + {} "
            "(positions)",
            lower, upper, n, static_cast<double>(mphfs[partition_id].num_bits()) / n,
            (positions[partition_id].bytes() * 8.0) / n);

        num_kmers_in_skew_index += n;
        lower = upper;
        upper = 2 * lower;
    }
    return num_kmers_in_skew_index;
}

void dictionary::print_space_breakdown() const {
    spdlog_piscem::info("total index size: {}[MB]",
                 essentials::convert((num_bits() + 7) / 8, essentials::MB));
    spdlog_piscem::info("SPACE BREAKDOWN:");
    spdlog_piscem::info("  minimizers: {} [bits/kmer]",
                 static_cast<double>(m_minimizers.num_bits()) / size());
    spdlog_piscem::info("  pieces: {} [bits/kmer]",
                 static_cast<double>(m_buckets.pieces.num_bits()) / size());
    spdlog_piscem::info("  num_super_kmers_before_bucket: {} [bits/kmer]",
                 static_cast<double>(m_buckets.num_super_kmers_before_bucket.num_bits()) / size());
    spdlog_piscem::info("  offsets: {} [bits/kmer]",
                 static_cast<double>(8 * m_buckets.offsets.bytes()) / size());
    spdlog_piscem::info("  strings: {} [bits/kmer]",
                 static_cast<double>(8 * m_buckets.strings.bytes()) / size());
    spdlog_piscem::info("  skew_index: {} [bits/kmer]",
                 static_cast<double>(m_skew_index.num_bits()) / size());
    spdlog_piscem::info("  weights: {} [bits/kmer]", static_cast<double>(m_weights.num_bits()) / size());
    m_weights.print_space_breakdown(size());
    spdlog_piscem::info("  --------------");
    spdlog_piscem::info("  total: {} [bits/kmer]", static_cast<double>(num_bits()) / size());
}

void dictionary::print_info() const {
    spdlog_piscem::info("=== dictionary info:");
    spdlog_piscem::info("num_kmers = {}", size());
    spdlog_piscem::info("k = {}", k());
    spdlog_piscem::info("num_minimizers = {}", m_minimizers.size());
    spdlog_piscem::info("m = {}", m());
    spdlog_piscem::info("canonicalized = {}", (canonicalized() ? "true" : "false"));
    spdlog_piscem::info("weighted = {}", (weighted() ? "true" : "false"));

    spdlog_piscem::info("num_super_kmers = {}", m_buckets.offsets.size());
    spdlog_piscem::info("num_pieces = {} (+{} [bits/kmer])", m_buckets.pieces.size(),
                 (2.0 * m_buckets.pieces.size() * (k() - 1)) / size());
    spdlog_piscem::info("bits_per_offset = ceil(log2({})) = {}", m_buckets.strings.size() / 2,
                 std::ceil(std::log2(m_buckets.strings.size() / 2)));
    uint64_t num_kmers_in_skew_index = m_skew_index.print_info();
    spdlog_piscem::info("num_kmers_in_skew_index {} ({}%)", num_kmers_in_skew_index,
                 (num_kmers_in_skew_index * 100.0) / size());
    print_space_breakdown();
}

}  // namespace sshash
