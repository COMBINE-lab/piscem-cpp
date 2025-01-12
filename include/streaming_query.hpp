#ifndef STREAMING_QUERY_HPP
#define STREAMING_QUERY_HPP

#include "../external/sshash/include/query/streaming_query_canonical_parsing.hpp"
#include "../external/sshash/include/dictionary.hpp"
#include "../external/sshash/include/util.hpp"
#include "CanonicalKmerIterator.hpp"

namespace piscem {

constexpr int32_t invalid_query_offset = std::numeric_limits<int32_t>::lowest();

class streaming_query {
public:
  streaming_query(sshash::dictionary const* dict) : m_q(dict), m_prev_query_offset(invalid_query_offset) {}
  
  void reset_state() {
    m_prev_query_offset = invalid_query_offset;
    m_q.reset_state();
  }

  void start() {
    m_prev_query_offset = invalid_query_offset;
    m_q.start();
  }

  sshash::lookup_result query_lookup(pufferfish::CanonicalKmerIterator& kmit) {
        auto query_offset = kmit->second;
        // if the current query offset position is
        // the next position after the stored query
        // offset position, then we can apply the
        // relevant optimizations.  Otherwise, we
        // should consider this as basically a "new"
        // query
        if (m_prev_query_offset != std::numeric_limits<int32_t>::lowest()) {
          if ((m_prev_query_offset + 1) != query_offset) { reset_state(); }
        }
        m_prev_query_offset = query_offset;
        const char* kmer = kmit.seq().substr(query_offset).data();
        return m_q.lookup_advanced(kmer);
  }


private:
  sshash::streaming_query_canonical_parsing m_q;
  int32_t m_prev_query_offset;
};
}

#endif // STREAMING_QUERY_HPP

