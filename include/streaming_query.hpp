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
  inline streaming_query(
    sshash::dictionary const *d,
    std::unique_ptr<sshash::bit_vector_iterator> &&contig_it)
    : m_d(d), m_q(d), m_prev_query_offset(invalid_query_offset),
      m_prev_contig_id(invalid_contig_id), m_prev_kmer_id(invalid_query_offset),
      m_ref_contig_it(std::move(contig_it)), m_k(d->k()) {}

  streaming_query(const streaming_query &other) = delete;
  streaming_query(streaming_query &&other) = default;
  streaming_query &operator=(const streaming_query &other) = delete;
  streaming_query &operator=(streaming_query &&other) = default;

  ~streaming_query() {
    auto ns = m_n_search;
    auto ne = m_n_extend;
    if (ns > 0) {
      std::stringstream ss;
      ss << "\nnumber of searches = " << ns << ", number of extensions = " << ne
         << ", neighbors = " << m_neighbors << ", extension ratio = "
         << static_cast<double>(ne) / static_cast<double>(ns) << "\n";
      std::cerr << ss.str();
    }
  }

  inline void reset_state() {
    m_prev_query_offset = invalid_query_offset;
    m_is_present = false;
    m_start = true;
    m_remaining_contig_bases = 0;
    // m_q.reset_state();
  }

  inline void start() {
    m_prev_query_offset = invalid_query_offset;
    m_prev_kmer_id = invalid_query_offset;
    m_is_present = false;
    m_start = true;
    m_remaining_contig_bases = 0;
    // m_q.start();
  }

  inline void do_stateless_lookup(const char *kmer_s) {
    // lookup directly in the index without assuming any relation
    // to the previously queried k-mer
    m_n_search++;
    m_prev_res = m_d->lookup_advanced(kmer_s);
    m_direction = m_prev_res.kmer_orientation ? -1 : 1;
    m_prev_kmer_id = m_prev_res.kmer_id;
    m_is_present = (m_prev_res.kmer_id != sshash::constants::invalid_uint64);
    m_start = !m_is_present;
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
      do_stateless_lookup(kmer_s);
      // we found this k-mer via stateless lookup
      if (m_is_present) {
        set_remaining_contig_bases();
      }
    } else {
      // try to get the next k-mer
      uint64_t next_kmer_id = m_prev_kmer_id + m_direction;
      uint64_t kmer_offset = next_kmer_id + (m_prev_res.contig_id * (m_k - 1));
      m_ref_contig_it->at(2 * kmer_offset);
      auto ref_kmer = m_ref_contig_it->read(2 * m_k);
      auto match_type = kmit->first.isEquivalent(ref_kmer);
      m_is_present = (match_type != KmerMatchType::NO_MATCH);

      if (!m_is_present) {
        // the next k-mer was not what was expected
        reset_state();
        // we're doing a fresh lookup
        do_stateless_lookup(kmer_s);
        if (m_is_present) {
          set_remaining_contig_bases();
        } else {
          reset_state();
        }
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

  /*
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
    if (m_prev_query_offset != invalid_query_offset) {
      if ((m_prev_query_offset + 1) != query_offset) {
        reset_state();
      }
    }
    m_prev_query_offset = query_offset;

    const char *kmer = kmit.seq().data() + query_offset;
    sshash::lookup_result res = m_q.lookup_advanced(kmer);
    m_is_present = (res.kmer_id != sshash::constants::invalid_uint64);

    uint64_t neighbor_inc =
      m_is_present && ((res.kmer_id == m_prev_kmer_id + 1) ||
                       (res.kmer_id == m_prev_kmer_id - 1));
    m_neighbors += neighbor_inc;
    m_prev_kmer_id = res.kmer_id;
    if (m_is_present && (res.contig_id != m_prev_contig_id)) {
      auto start_pos = m_ctg_offsets.access(res.contig_id);
      auto end_pos = m_ctg_offsets.access(res.contig_id + 1);
      size_t len = end_pos - start_pos;
      m_ctg_span = {m_ctg_entries.at(start_pos),
                    m_ctg_entries.at(start_pos + len), len};
      m_prev_contig_id = res.contig_id;
    }
    return res;
  }
  */
  inline bool is_present() { return m_is_present; }

  inline sshash::util::contig_span contig_span() { return m_ctg_span; }

private:
  sshash::dictionary const *m_d;
  sshash::streaming_query_canonical_parsing m_q;

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
  std::unique_ptr<sshash::bit_vector_iterator> m_ref_contig_it{nullptr};
  int32_t m_remaining_contig_bases{0};
  uint64_t m_k;

  // ref_contig_it.at(2 * ph.globalPos_);
};
} // namespace piscem

#endif // STREAMING_QUERY_HPP
