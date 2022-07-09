#pragma once

#include "util.hpp"
#include "minimizers.hpp"
#include "buckets.hpp"
#include "skew_index.hpp"
#include "weights.hpp"

namespace mindex {
    class reference_index;
}

namespace sshash {

struct dictionary {
    dictionary() : m_size(0), m_seed(0), m_k(0), m_m(0), m_canonical_parsing(0) {}

    void build(std::string const& filename, build_configuration const& build_config);

    uint64_t size() const { return m_size; }
    uint64_t seed() const { return m_seed; }
    uint64_t k() const { return m_k; }
    uint64_t m() const { return m_m; }
    bool canonicalized() const { return m_canonical_parsing; }
    bool weighted() const { return !m_weights.empty(); }

    uint64_t lookup(char const* string_kmer, bool check_reverse_complement_too = true) const;
    uint64_t lookup_uint64(uint64_t uint64_kmer, bool check_reverse_complement_too = true) const;

    lookup_result lookup_advanced(char const* string_kmer,
                                  bool check_reverse_complement_too = true) const;
    lookup_result lookup_advanced_uint64(uint64_t uint64_kmer,
                                         bool check_reverse_complement_too = true) const;

    /* Return the number of kmers in contig. Since contigs do not have duplicates,
       the length of the contig is always size + k - 1. */
    uint64_t contig_size(uint64_t contig_id) const;
    // std::vector<uint32_t> contig_neighbours(uint64_t contig_id) const;

    uint64_t weight(uint64_t kmer_id) const;

    void access(uint64_t kmer_id, char* string_kmer) const;

    bool is_member(char const* string_kmer, bool check_reverse_complement_too = true) const;
    bool is_member_uint64(uint64_t uint64_kmer, bool check_reverse_complement_too = true) const;

    friend struct streaming_query_canonical_parsing;
    friend struct streaming_query_regular_parsing;
    friend class ::mindex::reference_index;
    
    streaming_query_report streaming_query_from_file(std::string const& filename,
                                                     bool multiline) const;

    struct iterator {
        iterator(dictionary const* ptr, uint64_t kmer_id = 0) {
            it = ptr->m_buckets.at(kmer_id, ptr->m_k, ptr->m_size);
        }

        bool has_next() const { return it.has_next(); }
        std::pair<uint64_t, std::string> next() { return it.next(); }

    private:
        typename buckets::iterator it;
    };

    iterator begin() const { return iterator(this); }

    iterator at(uint64_t kmer_id) const {
        assert(kmer_id < size());
        return iterator(this, kmer_id);
    }

    uint64_t num_bits() const;
    void print_info() const;
    void print_space_breakdown() const;

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visitor.visit(m_size);
        visitor.visit(m_seed);
        visitor.visit(m_k);
        visitor.visit(m_m);
        visitor.visit(m_canonical_parsing);
        visitor.visit(m_minimizers);
        visitor.visit(m_buckets);
        visitor.visit(m_skew_index);
        visitor.visit(m_weights);
    }

private:
    uint64_t m_size;
    uint64_t m_seed;
    uint16_t m_k;
    uint16_t m_m;
    uint16_t m_canonical_parsing;
    minimizers m_minimizers;
    buckets m_buckets;
    skew_index m_skew_index;
    weights m_weights;

    lookup_result lookup_uint64_regular_parsing(uint64_t uint64_kmer) const;
    lookup_result lookup_uint64_canonical_parsing(uint64_t uint64_kmer) const;
};

}  // namespace sshash