#pragma once

#include <vector>
#include "util.hpp"
#include "ef_sequence.hpp"

namespace sshash {

class equivalence_class_map {
public:
    template <typename Visitor>
    void visit(Visitor& visitor) {
        visitor.visit(m_tile_ec_ids);
        visitor.visit(m_label_list_offsets);
        visitor.visit(m_label_entries);
    }

    inline sshash::util::ec_span entries_for_ec(uint64_t ec_id) const {
        auto start_pos = m_label_list_offsets.access(ec_id);
        auto end_pos = m_label_list_offsets.access(ec_id + 1);
        size_t len = end_pos - start_pos;
        return {m_label_entries.at(start_pos), m_label_entries.at(start_pos + len), len};
    }

    inline sshash::util::ec_span entries_for_tile(uint64_t tile_id) const {
        return entries_for_ec(m_tile_ec_ids[tile_id]);
    }

    pthash::compact_vector m_tile_ec_ids;
    pthash::compact_vector m_label_entries;
    sshash::ef_sequence<false> m_label_list_offsets;
};

}  // namespace sshash
