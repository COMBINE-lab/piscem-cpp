#pragma once

#include <vector>
#include "util_piscem.hpp"
#include "../external/sshash/external/pthash/external/bits/include/elias_fano.hpp"

namespace sshash {

class equivalence_class_map {
public:
        template <typename Visitor>
    void visit(Visitor& visitor) const {
      visit_impl(visitor, *this);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
      visit_impl(visitor, *this);
    }
    template <typename Visitor, typename T>
    static void visit_impl(Visitor& visitor, T&& t) {
        visitor.visit(t.m_tile_ec_ids);
        visitor.visit(t.m_label_list_offsets);
        visitor.visit(t.m_label_entries);
    }

    

    inline sshash::util::ec_span entries_for_ec(uint64_t ec_id) const {
        auto start_pos = m_label_list_offsets.access(ec_id);
        auto end_pos = m_label_list_offsets.access(ec_id + 1);
        size_t len = end_pos - start_pos;
        return {m_label_entries.get_iterator_at(start_pos), m_label_entries.get_iterator_at(start_pos + len), len};
    }

    inline sshash::util::ec_span entries_for_tile(uint64_t tile_id) const {
        return entries_for_ec(m_tile_ec_ids[tile_id]);
    }

    inline uint64_t ec_for_tile(uint64_t tile_id) const {
        return m_tile_ec_ids[tile_id];
    }

    bits::compact_vector m_tile_ec_ids;
    bits::compact_vector m_label_entries;
    bits::elias_fano<false,false> m_label_list_offsets;
};

}  // namespace sshash
