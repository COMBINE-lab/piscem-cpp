#pragma once

#include "dictionary.hpp"
#include "basic_contig_table.hpp"
#include "../external/pthash/external/essentials/include/essentials.hpp"
#include "query/contig_info_query_canonical_parsing.cpp"
#include "bit_vector_iterator.hpp"
#include "CanonicalKmerIterator.hpp"
#include "projected_hits.hpp"

/*namespace sshash {
    class contig_info_query_canonical_parsing;
    class dictionary;
}*/
namespace mindex {
class reference_index {
public:
    reference_index(const std::string& basename) {
        std::string dict_name = basename+".sshash";
        essentials::load(m_dict, dict_name.c_str());
        std::string ctg_name = basename+".ctab";
        essentials::load(bct, ctg_name.c_str());
    }

    projected_hits query(pufferfish::CanonicalKmerIterator kmit, sshash::contig_info_query_canonical_parsing& q) {
        auto qres = q.get_contig_pos(kmit->first.fwWord(), kmit->first.rcWord(), kmit->second);

        constexpr uint64_t invalid_u64 = std::numeric_limits<uint64_t>::max();
        constexpr uint32_t invalid_u32 = std::numeric_limits<uint32_t>::max();

        if (qres.is_member) {
            auto start_pos = bct.m_ctg_offsets.access(qres.contig_id);
            auto end_pos = bct.m_ctg_offsets.access(qres.contig_id+1);
            size_t len = end_pos - start_pos;
            nonstd::span s(bct.m_ctg_entries.begin() + start_pos, len);

            uint32_t contig_id = (qres.contig_id > invalid_u32) ? invalid_u32 : static_cast<uint32_t>(qres.contig_id);
            uint32_t contig_offset = (qres.contig_offset > invalid_u32) ? invalid_u32 : static_cast<uint32_t>(qres.contig_offset);
            uint32_t contig_length = (qres.contig_length> invalid_u32) ? invalid_u32 : static_cast<uint32_t>(qres.contig_length);

            return projected_hits {
                contig_id,
                contig_offset,
                qres.is_forward,
                contig_length,
                qres.global_pos,
                static_cast<uint32_t>(m_dict.k()),
                s
            };
        } else {
            return {
                invalid_u32, invalid_u32, false, invalid_u32, invalid_u64, static_cast<uint32_t>(m_dict.k()), 
                nonstd::span<sshash::util::Position>()
            };
        }
    }   

    uint64_t k() const { return m_dict.k(); }
    const sshash::dictionary* get_dict() const { return &m_dict; }
    pthash::bit_vector& contigs() { return m_dict.m_buckets.strings; }
private:
    sshash::dictionary m_dict;
    sshash::basic_contig_table bct;
};
}