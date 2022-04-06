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

    uint64_t m_ref_len_bits;
    pthash::compact_vector m_ctg_entries;
    sshash::ef_sequence<false> m_ctg_offsets;
};

}  // namespace sshash