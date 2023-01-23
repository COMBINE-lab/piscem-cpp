#pragma once

#include <vector>
#include "util.hpp"
#include "ef_sequence.hpp"

namespace sshash {

class basic_contig_table {
public:
    template <typename Visitor>
    void visit(Visitor& visitor) {
        visitor.visit(m_ref_len_bits);
        visitor.visit(m_ctg_offsets);
        visitor.visit(m_ctg_entries);
    }
  
    sshash::util::contig_span contig_entries(uint64_t contig_id) {
      auto start_pos = m_ctg_offsets.access(contig_id);
      auto end_pos = m_ctg_offsets.access(contig_id + 1);
      size_t len = end_pos - start_pos;
      return {m_ctg_entries.at(start_pos), m_ctg_entries.at(start_pos + len), len};
    }

    uint64_t m_ref_len_bits;
    pthash::compact_vector m_ctg_entries;
    sshash::ef_sequence<false> m_ctg_offsets;
};

}  // namespace sshash
