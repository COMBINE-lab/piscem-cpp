#pragma once

#include "../spdlog_piscem/spdlog.h"
#include "../gz/zip_stream.hpp"

namespace sshash {

struct parse_data {
    parse_data(std::string const& tmp_dirname) : num_kmers(0), minimizers(tmp_dirname) {}
    uint64_t num_kmers;
    minimizers_tuples minimizers;
    compact_string_pool strings;
    weights::builder weights_builder;
};

void parse_file_from_cuttlefish(std::istream& is, parse_data& data,
                                build_configuration const& build_config) {
    uint64_t k = build_config.k;
    uint64_t m = build_config.m;
    uint64_t seed = build_config.seed;
    uint64_t max_num_kmers_in_super_kmer = k - m + 1;
    uint64_t block_size = 2 * k - m;  // max_num_kmers_in_super_kmer + k - 1

    if (max_num_kmers_in_super_kmer >= (1ULL << (sizeof(num_kmers_in_super_kmer_uint_type) * 8))) {
        throw std::runtime_error(
            "max_num_kmers_in_super_kmer " + std::to_string(max_num_kmers_in_super_kmer) +
            " does not fit into " + std::to_string(sizeof(num_kmers_in_super_kmer_uint_type) * 8) +
            " bits");
    }

    /* fit into the wanted number of bits */
    assert(max_num_kmers_in_super_kmer < (1ULL << (sizeof(num_kmers_in_super_kmer_uint_type) * 8)));

    compact_string_pool::builder builder(k);

    std::string sequence;
    uint64_t prev_minimizer = constants::invalid_uint64;

    uint64_t begin = 0;  // begin of parsed super_kmer in sequence
    uint64_t end = 0;    // end of parsed super_kmer in sequence
    uint64_t num_sequences = 0;
    uint64_t num_bases = 0;
    bool glue = false;

    auto append_super_kmer = [&]() {
        if (sequence.empty() or prev_minimizer == constants::invalid_uint64 or begin == end) return;

        assert(end > begin);
        char const* super_kmer = sequence.data() + begin;
        uint64_t size = (end - begin) + k - 1;
        assert(util::is_valid(super_kmer, size));

        /* if num_kmers_in_super_kmer > k - m + 1, then split the super_kmer into blocks */
        uint64_t num_kmers_in_super_kmer = end - begin;
        uint64_t num_blocks = num_kmers_in_super_kmer / max_num_kmers_in_super_kmer +
                              (num_kmers_in_super_kmer % max_num_kmers_in_super_kmer != 0);
        assert(num_blocks > 0);
        for (uint64_t i = 0; i != num_blocks; ++i) {
            uint64_t n = block_size;
            if (i == num_blocks - 1) n = size;
            uint64_t num_kmers_in_block = n - k + 1;
            assert(num_kmers_in_block <= max_num_kmers_in_super_kmer);
            data.minimizers.emplace_back(prev_minimizer, builder.offset, num_kmers_in_block);
            builder.append(super_kmer + i * max_num_kmers_in_super_kmer, n, glue);
            if (glue) {
                assert(data.minimizers.back().offset > k - 1);
                data.minimizers.back().offset -= k - 1;
            }
            size -= max_num_kmers_in_super_kmer;
            glue = true;
        }
    };

    uint64_t seq_len = 0;
    uint64_t sum_of_weights = 0;
    data.weights_builder.init();

    /* intervals of weights */
    uint64_t weight_value = constants::invalid_uint64;
    uint64_t weight_length = 0;

    while (!is.eof()) {
        std::getline(is, sequence);  // header sequence
        auto tsep = sequence.find('\t');
        sequence = sequence.substr(tsep + 1);
        if (sequence.size() < k) continue;

        if (++num_sequences % 100000 == 0) {
            spdlog_piscem::info("read {} sequences, {} bases, {}, kmers", num_sequences, num_bases,
                         data.num_kmers);
        }

        begin = 0;
        end = 0;
        glue = false;  // start a new piece
        prev_minimizer = constants::invalid_uint64;
        num_bases += sequence.size();

        if (build_config.weighted and seq_len != sequence.size()) {
            spdlog_piscem::critical("expected a sequence of length {}, but got one of length {}.", seq_len,
                             sequence.length());
            throw std::runtime_error("file is malformed");
        }

        while (end != sequence.size() - k + 1) {
            char const* kmer = sequence.data() + end;
            assert(util::is_valid(kmer, k));
            uint64_t uint64_kmer = util::string_to_uint64_no_reverse(kmer, k);
            uint64_t minimizer = util::compute_minimizer(uint64_kmer, k, m, seed);

            if (build_config.canonical_parsing) {
                uint64_t uint64_kmer_rc = util::compute_reverse_complement(uint64_kmer, k);
                uint64_t minimizer_rc = util::compute_minimizer(uint64_kmer_rc, k, m, seed);
                minimizer = std::min<uint64_t>(minimizer, minimizer_rc);
            }

            if (prev_minimizer == constants::invalid_uint64) prev_minimizer = minimizer;
            if (minimizer != prev_minimizer) {
                append_super_kmer();
                begin = end;
                prev_minimizer = minimizer;
                glue = true;
            }

            ++data.num_kmers;
            ++end;
        }

        append_super_kmer();
    }

    data.minimizers.finalize();
    builder.finalize();
    builder.build(data.strings);

    spdlog_piscem::info("read {} sequences, {} bases, {} kmers.", num_sequences, num_bases,
                 data.num_kmers);
    spdlog_piscem::info("num_super_kmers {}", data.strings.num_super_kmers());
    spdlog_piscem::info("num_pieces {} (+{} [bits/kmer])", data.strings.pieces.size(),
                 (2.0 * data.strings.pieces.size() * (k - 1)) / data.num_kmers);

    assert(data.strings.pieces.size() == num_sequences + 1);

    if (build_config.weighted) {
        spdlog_piscem::info("sum_of_weights {}", sum_of_weights);
        data.weights_builder.push_weight_interval(weight_value, weight_length);
        data.weights_builder.finalize(data.num_kmers);
    }
}

parse_data parse_file(std::string const& filename, build_configuration const& build_config) {
    std::ifstream is(filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + filename + "'");
    spdlog_piscem::info("reading file '{}'...", filename);
    parse_data data(build_config.tmp_dirname);
    if (util::ends_with(filename, ".gz")) {
        zip_istream zis(is);
        parse_file_from_cuttlefish(zis, data, build_config);
    } else {
        parse_file_from_cuttlefish(is, data, build_config);
    }
    is.close();
    return data;
}

}  // namespace sshash
