#ifndef STREAMING_QUERY_HPP
#define STREAMING_QUERY_HPP

#include "../external/sshash/external/pthash/external/essentials/include/essentials.hpp"
#include "../external/sshash/include/dictionary.hpp"
#include "../external/sshash/include/query/streaming_query_canonical_parsing.hpp"
#include "../external/sshash/include/util.hpp"
#include "../include/unordered_dense.h"
// #include "../include/bcf/cuckoofilter.h"
#include "../include/boost/unordered/concurrent_flat_map.hpp"
#include "CanonicalKmerIterator.hpp"
#include "util_piscem.hpp"
#include <sstream>
#include <limits>

namespace mindex {
class reference_index;
}

namespace piscem {

constexpr int32_t invalid_query_offset = std::numeric_limits<int32_t>::lowest();
constexpr uint64_t invalid_contig_id = std::numeric_limits<uint64_t>::max();

struct empty_map_t {};
struct empty_filter_t {
  empty_filter_t(size_t _s) { (void)_s; }
};

// if `with_cache` is `true`, then this class will 
// use a cache for retaining the end-of-unitig k-mers 
// to speed up lookup. If it is instantiated with `false`, 
// no such cache will be used.
template <bool with_cache>
class streaming_query {
  using cache_t = boost::concurrent_flat_map<uint64_t, sshash::lookup_result>;
  //using cache_t = std::conditional_t<with_cache, ankerl::unordered_dense::map<uint64_t, sshash::lookup_result>, empty_map_t>;
  //using filter_t = std::conditional_t<with_cache, cuckoofilter::CuckooFilter<uint64_t, 12>, empty_filter_t>;
public:
  
  inline streaming_query(sshash::dictionary const *d, piscem::unitig_end_cache_t* unitig_end_cache = nullptr)
    : m_d(d), m_prev_query_offset(invalid_query_offset),
      m_prev_contig_id(invalid_contig_id), m_prev_kmer_id(sshash::constants::invalid_uint64),
      /*m_unitig_ends_filter(m_max_cache_size), */
      m_unitig_ends((unitig_end_cache == nullptr) ? nullptr : unitig_end_cache->get_map()),
      m_max_cache_size((unitig_end_cache == nullptr) ? 0 : unitig_end_cache->get_capacity()),
      m_ref_contig_it(sshash::bit_vector_iterator(d->strings(), 0)),
      m_k(d->k()) {}

  streaming_query(const streaming_query &other) = delete;
  streaming_query(streaming_query &&other) = default;
  streaming_query &operator=(const streaming_query &other) = delete;
  streaming_query &operator=(streaming_query &&other) = default;

  ~streaming_query() {
    if constexpr (with_cache && m_print_stats) {
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
           << ", extension ratio = "
           << static_cast<double>(ne) / static_cast<double>(ns) << "\n";
        if constexpr (with_cache) {
           ss << ", cache size = " << m_unitig_ends->size() << "\n"
              << ", cache hits = " << m_num_cache_hits << "\n";
        }

        std::cerr << ss.str();
      }
    }
  }

  inline void reset_state() {
    m_prev_query_offset = invalid_query_offset;
    m_prev_kmer_id = sshash::constants::invalid_uint64;
    m_is_present = false;
    m_start = true;
    m_remaining_contig_bases = 0;
  }

  inline void start() { reset_state(); }

  inline void do_stateless_lookup(const char *kmer_s, CanonicalKmer& kmer) {
    // lookup directly in the index without assuming any relation
    // to the previously queried k-mer
    bool was_cached = false;
    uint64_t canon_kmer = kmer.getCanonicalWord();
    bool fw_is_canonical = kmer.isFwCanonical();
   
    if constexpr (with_cache) {
      // first check the end cache if we are looking in a relevant place
      if (m_cache_end) {
        // ankerl map
        /*
        auto cache_it = m_unitig_ends.find(kmer.getCanonicalWord());
        if (cache_it != m_unitig_ends.end()) {
          m_prev_res = cache_it->second;
          m_num_cache_hits++;

          // marked as 1 if, when we looked up in the actual index, 
          // the forward k-mer was the canonical k-mer, and 0 otherwise; 
          uint64_t inserted_fw = (m_prev_res.kmer_orientation & 0x3) >> 1;
          m_prev_res.kmer_orientation = (m_prev_res.kmer_orientation & 0x1);
          if (inserted_fw != fw_is_canonical) {
            m_prev_res.kmer_orientation = 1 - m_prev_res.kmer_orientation; 
          }
          was_cached = true;
        }
        */

        // boost concurrent
        // num_visited is 1 if we found the value matching the key in the 
        // map and 0 otherwise.
        size_t num_visited = m_unitig_ends->cvisit(canon_kmer,
          [this](const auto& x) { this->m_prev_res = x.second; });
        if (num_visited == 1) {
          m_num_cache_hits++;
          // marked as 1 if, when we looked up in the actual index, 
          // the forward k-mer was the canonical k-mer, and 0 otherwise; 
          bool fw_was_canonical = (((m_prev_res.kmer_orientation & 0x2) >> 1) == 1) ? true : false;
          m_prev_res.kmer_orientation = (m_prev_res.kmer_orientation & 0x1);
          if (fw_was_canonical != fw_is_canonical) {
            m_prev_res.kmer_orientation = 1 - m_prev_res.kmer_orientation; 
          }
          was_cached = true;
        }
      }

      if (!was_cached) { 
        m_prev_res = m_d->lookup_advanced(kmer_s); 
      }
    } else {
      m_prev_res = m_d->lookup_advanced(kmer_s); 
    }

    m_direction = m_prev_res.kmer_orientation ? -1 : 1;
    m_prev_kmer_id = m_prev_res.kmer_id;
    m_is_present = (m_prev_res.kmer_id != sshash::constants::invalid_uint64);
    m_start = !m_is_present;
    m_n_search += m_is_present ? 1 : 0;
    if (m_is_present) {
      if constexpr(with_cache) {
        /*
        // starts off as false until the filter is full,
        // then remains true (but our cache is bounded size so
        // this is OK).
        bool do_insert = !m_use_filter;
        if (!was_cached && m_use_filter) {
          if (m_unitig_ends_filter.Contain(canon_kmer) == cuckoofilter::Ok) {
            // key was in the filter already so now add to the cache
            do_insert = true;
          } else if (m_unitig_ends_filter.Add(canon_kmer) != cuckoofilter::Ok) {
            // insert failed so stop trying to use the filter
            m_use_filter = false;
          }
          if (m_unitig_ends_filter.Size() >= m_max_cache_size) {
            m_use_filter = false;
          }
        }
        */
        if (!was_cached && m_cache_end && m_unitig_ends->size() < m_max_cache_size) {
          auto res_copy = m_prev_res;
          res_copy.kmer_orientation |= fw_is_canonical ? 0x2 : 0x0;
          // boost concurrent
          m_unitig_ends->try_emplace_or_cvisit(canon_kmer, std::move(res_copy), [](const auto& x) { (void)x; });
          // ankerl hash
          //m_unitig_ends[kmer.getCanonicalWord()] = res_copy;
        }
      }
      uint64_t kmer_offset =
        2 * (m_prev_res.kmer_id + (m_prev_res.contig_id * (m_k - 1)));
      kmer_offset += (m_direction > 0) ? 0 : (2 * m_k);
      m_ref_contig_it.at(kmer_offset);
      set_remaining_contig_bases();
    } else {
      reset_state();
    }
    m_cache_end = false;
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
    int32_t query_advance = (query_offset > m_prev_query_offset) ?
      (query_offset - m_prev_query_offset) : (reset_state(), invalid_query_offset);

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
      if (m_remaining_contig_bases < query_advance) {
        reset_state();
        m_cache_end = true;
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
      do_stateless_lookup(kmer_s, kmit->first);
    } else {
      // try to get the next k-mer
      uint64_t next_kmer_id = m_prev_kmer_id + (m_direction * query_advance);
      // can't "eat" more than 63 bits, so make sure we don't attempt to do so
      auto ref_kmer =
        (query_advance <= 31) ? 
        ((m_direction > 0)
          ? (m_ref_contig_it.eat(2 * query_advance), m_ref_contig_it.read(2 * m_k))
          : (m_ref_contig_it.eat_reverse(2 * query_advance),
             m_ref_contig_it.read_reverse(2 * m_k))) : 
        // if we can't just shift to the new k-mer, then use the `at` method;
        (m_ref_contig_it.at(next_kmer_id), m_ref_contig_it.read(2 * m_k));
      
      auto match_type = kmit->first.isEquivalent(ref_kmer);
      m_is_present = (match_type != KmerMatchType::NO_MATCH);
      if (!m_is_present) {
        // the next k-mer was not what was expected
        // we're doing a fresh lookup
        do_stateless_lookup(kmer_s, kmit->first);
      } else {
        // the next k-mer was what was expected
        m_n_extend++;
        m_remaining_contig_bases -= query_advance;
        m_start = false;
        m_prev_kmer_id = next_kmer_id;
        m_prev_res.kmer_id += (m_direction * query_advance);
        m_prev_res.kmer_id_in_contig += (m_direction * query_advance);

        // record the orientation of the previous match and look at the orientation
        // of the current match.
        auto prev_orientation = m_prev_res.kmer_orientation;
        m_direction = (match_type == KmerMatchType::IDENTITY_MATCH) ? 
          (m_prev_res.kmer_orientation = sshash::constants::forward_orientation, 1) : 
          (m_prev_res.kmer_orientation = sshash::constants::backward_orientation, -1);
        // if the orientation switched (should be exceedingly rare), recompute the 
        // distance to the end of the unitig because we have changed traversal direction.
        if (prev_orientation != m_prev_res.kmer_orientation) {
          // NOTE: This is a really strange case. We are walking along a contig in some direction
          // and we find the k-mer where we expect it, but the orientation is the opposite of what
          // we expect.  How, exactly does this happen?
          // In any case, it seems the right thing to do is to re-compute the distance to 
          // the end of the contig in the new orientation
          set_remaining_contig_bases();
        }
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
 
  
  //filter_t m_unitig_ends_filter;
  cache_t *m_unitig_ends;
  size_t m_max_cache_size{0};
  //cache_t m_unitig_ends;
  bool m_use_filter = false;//with_cache;
  bool m_cache_end = false;
  bool m_start = true;
  bool m_is_present = false;
  int32_t m_direction{0};
  int32_t m_n_search{0};
  int32_t m_n_extend{0};
  uint64_t m_num_cache_hits{0};
  sshash::lookup_result m_prev_res;
  sshash::util::contig_span m_ctg_span;
  sshash::bit_vector_iterator m_ref_contig_it;
  int32_t m_remaining_contig_bases{0};
  uint64_t m_k;
  static constexpr bool m_print_stats{false};
  //static constexpr uint64_t m_max_cache_size{5000000};
};
} // namespace piscem

#endif // STREAMING_QUERY_HPP
