#pragma once

#include <vector>
#include <unordered_map>  // count the distinct weights with freq information

#include "ef_sequence.hpp"
#include "spdlog/spdlog.h"

namespace sshash {

struct weights {
    struct builder {
        builder() {}

        void init() { m_weight_interval_lengths.push_back(0); }

        void eat(uint64_t weight) {
            assert(weight > 0);
            auto it = m_weights_map.find(weight);
            if (it != m_weights_map.cend()) {  // found
                (*it).second += 1;
            } else {
                m_weights_map[weight] = 1;
            }
        }

        void push_weight_interval(uint64_t value, uint64_t length) {
            m_weight_interval_values.push_back(value);
            m_weight_interval_lengths.push_back(m_weight_interval_lengths.back() + length);
        }

        uint64_t num_weight_intervals() const { return m_weight_interval_values.size(); }

        void finalize(uint64_t num_kmers) {
            assert(
                std::is_sorted(m_weight_interval_lengths.begin(), m_weight_interval_lengths.end()));

            spdlog::info("num_weight_intervals {}", num_weight_intervals());

            uint64_t num_distinct_weights = m_weights_map.size();

            spdlog::info("found {} distinct weights (ceil(log2({})) = {})", num_distinct_weights, num_distinct_weights, std::ceil(std::log2(num_distinct_weights)));

            m_weights.reserve(num_distinct_weights);
            uint64_t n = 0;
            uint64_t largest_weight = 0;
            for (auto p : m_weights_map) {
                if (p.first > largest_weight) largest_weight = p.first;
                n += p.second;
                m_weights.push_back(p);
            }
            assert(largest_weight > 0);

            spdlog::info("largest_weight+1 = {} (ceil(log2(largest_weight + 1)) = {}", largest_weight+1, std::ceil(std::log2(largest_weight + 1)));

            if (n != num_kmers) {
              spdlog::critical("expected {} kmers but got {}", num_kmers, n);
                throw std::runtime_error("file is malformed");
            }

            std::sort(m_weights.begin(), m_weights.end(), [](auto const& x, auto const& y) {
                if (x.second != y.second) return x.second > y.second;
                return x.first < y.first;
            });

            uint64_t rest = num_kmers - m_weights.front().second;
            spdlog::info("kmers that do not have the most frequent weight: {} ({}%)", rest, (rest * 100.0) / num_kmers); 

            m_weight_dictionary_builder.resize(num_distinct_weights,
                                               std::ceil(std::log2(largest_weight + 1)));
            for (uint64_t id = 0; id != num_distinct_weights; ++id) {
                uint64_t weight_value = m_weights[id].first;
                m_weight_dictionary_builder.set(id, weight_value);
                m_weights_map[weight_value] = id;
            }
        }

        void build(weights& index) {
            uint64_t num_distinct_weights = m_weight_dictionary_builder.size();
            pthash::compact_vector::builder weight_interval_values;
            weight_interval_values.resize(
                m_weight_interval_values.size(),
                num_distinct_weights == 1 ? 1 : std::ceil(std::log2(num_distinct_weights)));
            uint64_t prev_weight = constants::invalid_uint64;
            for (uint64_t i = 0; i != m_weight_interval_values.size(); ++i) {
                uint64_t weight = m_weight_interval_values[i];
                if (i != 0) {
                    if (weight == prev_weight) {
                        std::cerr << "Error at " << i << "/" << m_weight_interval_values.size()
                                  << ": cannot have two consecutive intervals with the same "
                                     "weight value"
                                  << std::endl;
                        throw std::runtime_error("weight intervals are malformed");
                    }
                }
                prev_weight = weight;
                uint64_t id = m_weights_map[weight];
                assert(id < num_distinct_weights);
                weight_interval_values.set(i, id);
            }
            weight_interval_values.build(index.m_weight_interval_values);
            index.m_weight_interval_lengths.encode(m_weight_interval_lengths.begin(),
                                                   m_weight_interval_lengths.size(),
                                                   m_weight_interval_lengths.back());

            m_weight_dictionary_builder.build(index.m_weight_dictionary);
        }

        /* return the average empirical entropy per weight */
        double print_info(uint64_t num_kmers) {
            assert(!m_weights.empty());
            double expected_weight = 0.0;
            double entropy_weights = 0.0;
            uint64_t print = 0;
            for (auto p : m_weights) {
                double prob = static_cast<double>(p.second) / num_kmers;
                expected_weight += p.first * prob;
                entropy_weights += prob * std::log2(1.0 / prob);
                print += 1;
                if (print <= 10) {
                  spdlog::info("weight: {} freq: {} ({}%)", p.first, p.second, (p.second * 100.0) / num_kmers); 
                }
            }
            spdlog::info("expected_weight {}", expected_weight); 
            spdlog::info("entropy_weights {} [bits/kmer]", entropy_weights); 
            return entropy_weights;
        }

    private:
        /* (weight,frequency) pairs during construction, then (weight,id) after sorting */
        std::unordered_map<uint64_t, uint64_t> m_weights_map;
        std::vector<std::pair<uint64_t, uint64_t>> m_weights;  // (weight,frequency)

        std::vector<uint64_t> m_weight_interval_values;
        std::vector<uint64_t> m_weight_interval_lengths;

        pthash::compact_vector::builder m_weight_dictionary_builder;
    };

    bool empty() const { return m_weight_dictionary.size() == 0; }

    uint64_t weight(uint64_t kmer_id) const {
        uint64_t i = m_weight_interval_lengths.prev_leq(kmer_id);
        uint64_t id = m_weight_interval_values.access(i);
        uint64_t weight = m_weight_dictionary.access(id);
        return weight;
    }

    uint64_t num_bits() const {
        return m_weight_interval_values.bytes() * 8 + m_weight_interval_lengths.num_bits() +
               m_weight_dictionary.bytes() * 8;
    }

    void print_space_breakdown(uint64_t num_kmers) const {
      spdlog::info("    weight_interval_values: {} [bits/kmer]", static_cast<double>(m_weight_interval_values.bytes() * 8) / num_kmers);
      spdlog::info("    weight_interval_lengths: {} [bits/kmer]", static_cast<double>(m_weight_interval_lengths.num_bits()) / num_kmers);
      spdlog::info("    weight_dictionary: {} [bits/kmer]", static_cast<double>(m_weight_dictionary.bytes() * 8) / num_kmers);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visitor.visit(m_weight_interval_values);
        visitor.visit(m_weight_interval_lengths);
        visitor.visit(m_weight_dictionary);
    }

private:
    pthash::compact_vector m_weight_interval_values;
    ef_sequence<true> m_weight_interval_lengths;
    pthash::compact_vector m_weight_dictionary;
};

}  // namespace sshash
