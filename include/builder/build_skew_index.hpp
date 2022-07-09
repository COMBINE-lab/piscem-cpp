#pragma once

#include "../../external/pthash/include/pthash.hpp"
#include "../spdlog/spdlog.h"

namespace sshash {

void build_skew_index(skew_index& m_skew_index, parse_data& data, buckets const& m_buckets,
                      build_configuration const& build_config,
                      buckets_statistics const& buckets_stats) {
    uint64_t min_log2_size = m_skew_index.min_log2;
    uint64_t max_log2_size = m_skew_index.max_log2;

    uint64_t max_num_super_kmers_in_bucket = buckets_stats.max_num_super_kmers_in_bucket();
    m_skew_index.log2_max_num_super_kmers_in_bucket =
        std::ceil(std::log2(buckets_stats.max_num_super_kmers_in_bucket()));

    spdlog::info("max_num_super_kmers_in_bucket {}", max_num_super_kmers_in_bucket);
    spdlog::info("log2_max_num_super_kmers_in_bucket {}",
                 m_skew_index.log2_max_num_super_kmers_in_bucket);

    mm::file_source<minimizer_tuple> input(data.minimizers.get_minimizers_filename(),
                                           mm::advice::sequential);

    uint64_t num_buckets_in_skew_index = 0;
    uint64_t num_super_kmers_in_skew_index = 0;
    for (minimizers_tuples_iterator it(input.data(), input.data() + input.size()); it.has_next();
         it.next()) {
        uint64_t list_size = it.list().size();
        if (list_size > (1ULL << min_log2_size)) {
            num_super_kmers_in_skew_index += list_size;
            ++num_buckets_in_skew_index;
        }
    }
    spdlog::info("num_buckets_in_skew_index {} / {} ({}%)", num_buckets_in_skew_index,
                 buckets_stats.num_buckets(),
                 (num_buckets_in_skew_index * 100.0) / buckets_stats.num_buckets());

    if (num_buckets_in_skew_index == 0) {
        input.close();
        return;
    }

    std::vector<list_type> lists;
    lists.reserve(num_buckets_in_skew_index);
    std::vector<minimizer_tuple> lists_tuples;  // backed memory
    lists_tuples.reserve(num_super_kmers_in_skew_index);
    for (minimizers_tuples_iterator it(input.data(), input.data() + input.size()); it.has_next();
         it.next()) {
        auto list = it.list();
        if (list.size() > (1ULL << min_log2_size)) {
            minimizer_tuple const* begin = lists_tuples.data() + lists_tuples.size();
            std::copy(list.begin_ptr(), list.end_ptr(), std::back_inserter(lists_tuples));
            minimizer_tuple const* end = lists_tuples.data() + lists_tuples.size();
            lists.push_back(list_type(begin, end));
        }
    }
    assert(lists.size() == num_buckets_in_skew_index);
    input.close();

    std::sort(lists.begin(), lists.end(),
              [](list_type const& x, list_type const& y) { return x.size() < y.size(); });

    uint64_t num_partitions = max_log2_size - min_log2_size + 1;
    if (buckets_stats.max_num_super_kmers_in_bucket() < (1ULL << max_log2_size)) {
        num_partitions = m_skew_index.log2_max_num_super_kmers_in_bucket - min_log2_size;
    }
    spdlog::info("num_partitions {}", num_partitions);

    std::vector<uint64_t> num_kmers_in_partition(num_partitions, 0);
    m_skew_index.mphfs.resize(num_partitions);
    m_skew_index.positions.resize(num_partitions);

    {
        spdlog::info("computing partitions...");

        uint64_t partition_id = 0;
        uint64_t lower = 1ULL << min_log2_size;
        uint64_t upper = 2 * lower;
        uint64_t num_kmers_in_skew_index = 0;
        for (uint64_t i = 0; i != lists.size() + 1; ++i) {
            if (i == lists.size() or lists[i].size() > upper) {
                spdlog::info("num_kmers belonging to buckets of size > {} and <= {}: {}", lower,
                             upper, num_kmers_in_partition[partition_id]);

                if (num_kmers_in_partition[partition_id] == 0) {
                    spdlog::critical(
                        "==> Empty bucket detected:\nthere is no k-mer that belongs to a list of "
                        "size > {} and <= {}",
                        lower, upper);
                    throw empty_bucket_runtime_error();
                }
                util::check_hash_collision_probability(num_kmers_in_partition[partition_id]);
                num_kmers_in_skew_index += num_kmers_in_partition[partition_id];
                partition_id += 1;

                if (i == lists.size()) break;

                lower = upper;
                upper = 2 * lower;
                if (partition_id == num_partitions - 1) upper = max_num_super_kmers_in_bucket;

                /*
                    If still larger than upper, then there is an empty bucket
                    and we should try different parameters.
                */
                if (lists[i].size() > upper) {
                    spdlog::critical(
                        "==> Empty bucket detected:\nthere is no k-mer that belongs to a list of "
                        "size > {} and <= {}",
                        lower, upper);
                    throw empty_bucket_runtime_error();
                }
            }

            assert(lists[i].size() > lower and lists[i].size() <= upper);
            for (auto [offset, num_kmers_in_super_kmer] : lists[i]) {
                (void)offset;  // unused
                num_kmers_in_partition[partition_id] += num_kmers_in_super_kmer;
            }
        }
        assert(partition_id == num_partitions);
        spdlog::info("num_kmers_in_skew_index {} ({}%)", num_kmers_in_skew_index,
                     (num_kmers_in_skew_index * 100.0) / buckets_stats.num_kmers());

        assert(num_kmers_in_skew_index == std::accumulate(num_kmers_in_partition.begin(),
                                                          num_kmers_in_partition.end(),
                                                          uint64_t(0)));
    }

    {
        pthash::build_configuration mphf_config;
        mphf_config.c = build_config.c;
        mphf_config.alpha = 0.94;
        mphf_config.seed = 1234567890;  // my favourite seed
        mphf_config.minimal_output = true;
        mphf_config.verbose_output = false;
        mphf_config.num_threads = std::thread::hardware_concurrency() >= 8 ? 8 : 1;

        spdlog::info("building PTHash mphfs (with {} threads) and positions...",
                     mphf_config.num_threads);

        uint64_t partition_id = 0;
        uint64_t lower = 1ULL << min_log2_size;
        uint64_t upper = 2 * lower;
        uint64_t num_bits_per_pos = min_log2_size + 1;

        /* tmp storage for keys and super_kmer_ids ******/
        std::vector<uint64_t> keys_in_partition;
        std::vector<uint32_t> super_kmer_ids_in_partition;
        keys_in_partition.reserve(num_kmers_in_partition[partition_id]);
        super_kmer_ids_in_partition.reserve(num_kmers_in_partition[partition_id]);
        pthash::compact_vector::builder cvb_positions;
        cvb_positions.resize(num_kmers_in_partition[partition_id], num_bits_per_pos);
        /*******/

        for (uint64_t i = 0; i != lists.size() + 1; ++i) {
            if (i == lists.size() or lists[i].size() > upper) {
                spdlog::info("lower {}; upper{}; num_bits_per_pos {}", lower, upper,
                             num_bits_per_pos);

                auto& mphf = m_skew_index.mphfs[partition_id];
                assert(num_kmers_in_partition[partition_id] == keys_in_partition.size());
                assert(num_kmers_in_partition[partition_id] == super_kmer_ids_in_partition.size());
                mphf.build_in_internal_memory(keys_in_partition.begin(), keys_in_partition.size(),
                                              mphf_config);

                spdlog::info("  built mphs[{}] for {} keys; bits/key = {}", partition_id,
                             keys_in_partition.size(),
                             static_cast<double>(mphf.num_bits()) / mphf.num_keys());

                for (uint64_t i = 0; i != keys_in_partition.size(); ++i) {
                    uint64_t kmer = keys_in_partition[i];
                    uint64_t pos = mphf(kmer);
                    uint32_t super_kmer_id = super_kmer_ids_in_partition[i];
                    cvb_positions.set(pos, super_kmer_id);
                }
                auto& positions = m_skew_index.positions[partition_id];
                cvb_positions.build(positions);

                spdlog::info("  built positions[{}] for {} keys; bits/key = {}", partition_id, positions.size(), (positions.bytes() * 8.0) / positions.size());

                partition_id += 1;

                if (i == lists.size()) break;

                lower = upper;
                upper = 2 * lower;
                num_bits_per_pos += 1;
                if (partition_id == num_partitions - 1) {
                    upper = max_num_super_kmers_in_bucket;
                    num_bits_per_pos = m_skew_index.log2_max_num_super_kmers_in_bucket;
                }

                keys_in_partition.clear();
                super_kmer_ids_in_partition.clear();
                keys_in_partition.reserve(num_kmers_in_partition[partition_id]);
                super_kmer_ids_in_partition.reserve(num_kmers_in_partition[partition_id]);
                cvb_positions.resize(num_kmers_in_partition[partition_id], num_bits_per_pos);
            }

            assert(lists[i].size() > lower and lists[i].size() <= upper);
            uint64_t super_kmer_id = 0;
            for (auto [offset, num_kmers_in_super_kmer] : lists[i]) {
                bit_vector_iterator bv_it(m_buckets.strings, 2 * offset);
                for (uint64_t i = 0; i != num_kmers_in_super_kmer; ++i) {
                    uint64_t kmer = bv_it.read(2 * build_config.k);
                    keys_in_partition.push_back(kmer);
                    super_kmer_ids_in_partition.push_back(super_kmer_id);
                    bv_it.eat(2);
                }
                assert(super_kmer_id < (1ULL << cvb_positions.width()));
                ++super_kmer_id;
            }
        }
        assert(partition_id == num_partitions);
    }

    spdlog::info("num_bits_for_skew_index {} ({} [bits/kmer])", m_skew_index.num_bits(),
                 static_cast<double>(m_skew_index.num_bits()) / buckets_stats.num_kmers());
}

}  // namespace sshash
