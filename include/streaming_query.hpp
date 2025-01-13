#ifndef STREAMING_QUERY_HPP
#define STREAMING_QUERY_HPP

#include "../external/sshash/external/pthash/external/essentials/include/essentials.hpp"
#include "../external/sshash/include/dictionary.hpp"
#include "../external/sshash/include/query/streaming_query_canonical_parsing.hpp"
#include "../external/sshash/include/util.hpp"
#include "CanonicalKmerIterator.hpp"
#include "util_piscem.hpp"
#include <sstream>

namespace mindex {
class reference_index;
}

namespace piscem {

constexpr int32_t invalid_query_offset = std::numeric_limits<int32_t>::lowest();
constexpr uint64_t invalid_contig_id = std::numeric_limits<uint64_t>::max();

class streaming_query {
public:
  inline streaming_query(sshash::dictionary const *d)
    : m_d(d), m_prev_query_offset(invalid_query_offset),
      m_prev_contig_id(invalid_contig_id), m_prev_kmer_id(invalid_query_offset),
      m_ref_contig_it(sshash::bit_vector_iterator(d->strings(), 0)),
      m_k(d->k()) {}

  streaming_query(const streaming_query &other) = delete;
  streaming_query(streaming_query &&other) = default;
  streaming_query &operator=(const streaming_query &other) = delete;
  streaming_query &operator=(streaming_query &&other) = default;

  ~streaming_query() {
    if constexpr (m_print_stats) {
      // print out statistics about the number of
      // extension lookups compared to the number
      // of stateless lookups.
      auto ns = num_searches();   // m_n_search;
      auto ne = num_extensions(); // nm_n_extend;
      if (ns > 0) {
        std::stringstream ss;
        ss << "method : " << search_meth() << "\n";
        ss << "\nnumber of searches = " << ns
           << ", number of extensions = " << ne
           << ", neighbors = " << m_neighbors << ", extension ratio = "
           << static_cast<double>(ne) / static_cast<double>(ns)
           << ", neihbor ratio = "
           << static_cast<double>(m_neighbors) / static_cast<double>(ns)
           << "\n";
        std::cerr << ss.str();
      }
    }
  }

  inline void reset_state() {
    m_prev_query_offset = invalid_query_offset;
    m_prev_kmer_id = invalid_query_offset;
    m_is_present = false;
    m_start = true;
    m_remaining_contig_bases = 0;
  }

  inline void start() { reset_state(); }

  inline void do_stateless_lookup(const char *kmer_s) {
    // lookup directly in the index without assuming any relation
    // to the previously queried k-mer
    m_prev_res = m_d->lookup_advanced(kmer_s);
    m_direction = m_prev_res.kmer_orientation ? -1 : 1;
    m_prev_kmer_id = m_prev_res.kmer_id;
    m_is_present = (m_prev_res.kmer_id != sshash::constants::invalid_uint64);
    m_start = !m_is_present;
    m_n_search += m_is_present ? 1 : 0;
    if (m_is_present) {
      uint64_t kmer_offset =
        2 * (m_prev_res.kmer_id + (m_prev_res.contig_id * (m_k - 1)));
      kmer_offset += (m_direction > 0) ? 0 : (2 * m_k);
      m_ref_contig_it.at(kmer_offset);
      set_remaining_contig_bases();
    } else {
      reset_state();
    }
  }

  inline void set_remaining_contig_bases() {
    // if moving forward, we have (contig-length - (pos + k)) positions left
    // if moving backward, we have (pos) positions left.
    m_remaining_contig_bases =
      (m_direction == 1)
        ? (m_prev_res.contig_size - (m_prev_res.kmer_id_in_contig + m_k))
        : (m_prev_res.kmer_id_in_contig);
  }

  inline sshash::lookup_result
  query_lookup(pufferfish::CanonicalKmerIterator &kmit,
               sshash::ef_sequence<false> &m_ctg_offsets,
               pthash::compact_vector &m_ctg_entries) {

    auto query_offset = kmit->second;

    // if the current query offset position is
    // the next position after the stored query
    // offset position, then we can apply the
    // relevant optimizations.  Otherwise, we
    // should consider this as basically a "new"
    // query
    //

    // if this isn't the first search after a reset
    if (m_prev_query_offset != invalid_query_offset) {
      // if there is no contig left to walk, or if
      // the current k-mer is not the successor of the
      // previous k-mer on the query, then reset the state.
      if (((m_prev_query_offset + 1) != query_offset) ||
          (m_remaining_contig_bases <= 0)) {
        reset_state();
      }
    }

    // get the k-mer string
    const char *kmer_s = kmit.seq().substr(query_offset).data();

    // if we need to do a stateless lookup
    if (m_start) {
      // NOTE: technically, the `CanonicalKmerIterator` should
      // never yield an invalid k-mer, so we shouldn't have to
      // actually check this.
      if (!sshash::util::is_valid(kmer_s, m_k)) {
        return sshash::lookup_result();
      }
      do_stateless_lookup(kmer_s);
    } else {
      // try to get the next k-mer
      uint64_t next_kmer_id = m_prev_kmer_id + m_direction;
      auto ref_kmer =
        (m_direction > 0)
          ? (m_ref_contig_it.eat(2), m_ref_contig_it.read(2 * m_k))
          : (m_ref_contig_it.eat_reverse(2),
             m_ref_contig_it.read_reverse(2 * m_k));
      auto match_type = kmit->first.isEquivalent(ref_kmer);
      m_is_present = (match_type != KmerMatchType::NO_MATCH);

      if (!m_is_present) {
        // the next k-mer was not what was expected
        // we're doing a fresh lookup
        do_stateless_lookup(kmer_s);
      } else {
        // the next k-mer was what was expected
        m_n_extend++;
        m_remaining_contig_bases--;
        m_start = false;
        m_prev_kmer_id = next_kmer_id;
        m_prev_res.kmer_id += m_direction;
        m_prev_res.kmer_id_in_contig += m_direction;
      }
    }

    // the query offset is unconditionally the
    // offset we just searched
    m_prev_query_offset = query_offset;

    // if we found the query, and the contig id is different
    // from that of the last found contig, then we have to refresh the
    // contig spant.
    if (m_is_present && (m_prev_res.contig_id != m_prev_contig_id)) {
      auto start_pos = m_ctg_offsets.access(m_prev_res.contig_id);
      auto end_pos = m_ctg_offsets.access(m_prev_res.contig_id + 1);
      size_t len = end_pos - start_pos;
      m_ctg_span = {m_ctg_entries.at(start_pos),
                    m_ctg_entries.at(start_pos + len), len};
      m_prev_contig_id = m_prev_res.contig_id;
    }
    return m_prev_res;
  }

  uint64_t num_searches() { return m_n_search; }

  uint64_t num_extensions() { return m_n_extend; }

  std::string search_meth() { return "custom"; }

  inline bool is_present() { return m_is_present; }

  inline sshash::util::contig_span contig_span() { return m_ctg_span; }

private:
  sshash::dictionary const *m_d;

  int32_t m_prev_query_offset;
  uint64_t m_prev_contig_id;
  uint64_t m_prev_kmer_id;
  uint64_t m_neighbors{0};

  bool m_start = true;
  bool m_is_present = false;
  int32_t m_direction{0};
  int32_t m_n_search{0};
  int32_t m_n_extend{0};
  sshash::lookup_result m_prev_res;
  sshash::util::contig_span m_ctg_span;
  sshash::bit_vector_iterator m_ref_contig_it;
  int32_t m_remaining_contig_bases{0};
  uint64_t m_k;
  static constexpr bool m_print_stats{false};
};
} // namespace piscem

#endif // STREAMING_QUERY_HPP
