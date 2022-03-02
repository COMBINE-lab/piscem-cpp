#pragma once

#include "dictionary.hpp"
#include "basic_contig_table.hpp"
#include "../external/pthash/external/essentials/include/essentials.hpp"
#include "query/contig_info_query_canonical_parsing.cpp"
#include "CanonicalKmerIterator.hpp"
#include "projected_hits.hpp"
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

        if (qres.is_member) {
            auto start_pos = bct.m_ctg_offsets.access(qres.contig_id);
            auto end_pos = bct.m_ctg_offsets.access(qres.contig_id+1);
            size_t len = end_pos - start_pos;
            nonstd::span s(bct.m_ctg_entries.begin() + start_pos, len);

            return projected_hits {
                static_cast<uint32_t>(qres.contig_id),
                static_cast<uint32_t>(qres.contig_offset),
                qres.is_forward,
                static_cast<uint32_t>(qres.contig_length),
                static_cast<uint32_t>(m_dict.k()),
                s
            };
        } else {
            return {
                0, 0, false, 0, static_cast<uint32_t>(m_dict.k()), 
                nonstd::span<sshash::util::Position>()
            };
        }
    }   

    uint64_t k() const { return m_dict.k(); }
    const sshash::dictionary* get_dict() const { return &m_dict; }
private:
    sshash::dictionary m_dict;
    sshash::basic_contig_table bct;
};