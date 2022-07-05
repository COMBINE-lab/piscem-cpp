#pragma once

#include <fstream>

#include "dictionary.hpp"
#include "basic_contig_table.hpp"
#include "../external/pthash/external/essentials/include/essentials.hpp"
#include "../include/bitsery/bitsery.h"
#include "../include/bitsery/adapter/stream.h"
#include "../include/bitsery/brief_syntax/vector.h"
#include "../include/bitsery/brief_syntax/string.h"
//#include "query/contig_info_query_canonical_parsing.cpp"
#include "query/streaming_query_canonical_parsing.hpp"
#include "bit_vector_iterator.hpp"
#include "CanonicalKmerIterator.hpp"
#include "projected_hits.hpp"
#include "util.hpp"
#include "spdlog/spdlog.h"

namespace mindex {
class reference_index {
public:
    reference_index(const std::string& basename) {
        spdlog::info("loading index from {}", basename);
        std::string dict_name = basename + ".sshash";
        essentials::load(m_dict, dict_name.c_str());
        std::string ctg_name = basename + ".ctab";
        essentials::load(m_bct, ctg_name.c_str());
        // based on the number of bits used to encode reference positions
        // read from the file, set the shift we have to perform on a
        // contig table entry to read off the reference id (= m_ref_len_bits + 1)
        // where the +1 is for the orientation bit.
        sshash::util::_ref_shift = (m_bct.m_ref_len_bits + 1);
        // based on the value of m_ref_len_bits, select the appropriate mask to use
        // when decoding a reference position from a contig table entry.
        sshash::util::_pos_mask = sshash::util::pos_masks[m_bct.m_ref_len_bits];

        std::string ref_info = basename + ".refinfo";

        std::fstream s{ref_info.c_str(), s.binary | s.in};
        auto state = bitsery::quickDeserialization<bitsery::InputStreamAdapter>(s, m_ref_names);
        state = bitsery::quickDeserialization<bitsery::InputStreamAdapter>(s, m_ref_lens);
        spdlog::info("done loading index");
    }

    projected_hits query(pufferfish::CanonicalKmerIterator kmit,
                         sshash::streaming_query_canonical_parsing& q) {
        auto qres = q.get_contig_pos(kmit->first.fwWord(), kmit->first.rcWord(), kmit->second);

        constexpr uint64_t invalid_u64 = std::numeric_limits<uint64_t>::max();
        constexpr uint32_t invalid_u32 = std::numeric_limits<uint32_t>::max();

        bool is_member = (qres.kmer_id != sshash::constants::invalid_uint64);

        // std::cout << "== ANSWER\n";
        // std::cout << "read position = " << kmit->second << "\n";
        // std::cout << "kmer_id " << qres.kmer_id << '\n';
        // std::cout << "kmer_id_in_contig " << qres.kmer_id_in_contig << '\n';
        // std::cout << "kmer_orientation " << qres.kmer_orientation << '\n';
        // std::cout << "contig_id " << qres.contig_id << '\n';
        // std::cout << "contig_size " << qres.contig_size << '\n';

        if (is_member) {
            qres.contig_size += m_dict.k() - 1;
            auto ctsize = m_bct.m_ctg_offsets.size();
            auto start_pos = m_bct.m_ctg_offsets.access(qres.contig_id);
            auto end_pos = m_bct.m_ctg_offsets.access(qres.contig_id + 1);
            size_t len = end_pos - start_pos;
            sshash::util::contig_span s{m_bct.m_ctg_entries.at(start_pos),
                                        m_bct.m_ctg_entries.at(start_pos + len), len};

            uint32_t contig_id = (qres.contig_id > invalid_u32)
                                     ? invalid_u32
                                     : static_cast<uint32_t>(qres.contig_id);
            uint32_t contig_offset = (qres.kmer_id_in_contig > invalid_u32)
                                         ? invalid_u32
                                         : static_cast<uint32_t>(qres.kmer_id_in_contig);
            uint32_t contig_length = (qres.contig_size > invalid_u32)
                                         ? invalid_u32
                                         : static_cast<uint32_t>(qres.contig_size);

            bool is_forward = (qres.kmer_orientation == sshash::constants::forward_orientation);

            return projected_hits{contig_id,
                                  contig_offset,
                                  is_forward,
                                  contig_length,
                                  qres.kmer_id,
                                  static_cast<uint32_t>(m_dict.k()),
                                  s};
        } else {
            return {invalid_u32, invalid_u32, false,
                    invalid_u32, invalid_u64, static_cast<uint32_t>(m_dict.k()),
                    {}};
        }
    }

    uint64_t k() const { return m_dict.k(); }
    const sshash::dictionary* get_dict() const { return &m_dict; }
    pthash::bit_vector& contigs() { return m_dict.m_buckets.strings; }
    const std::string& ref_name(size_t i) const { return m_ref_names[i]; }
    uint64_t ref_len(size_t i) const { return m_ref_lens[i]; }
    uint64_t num_refs() const { return m_ref_names.size(); }
    const sshash::basic_contig_table& get_contig_table() const { return m_bct; }

private:
    sshash::dictionary m_dict;
    sshash::basic_contig_table m_bct;
    std::vector<std::string> m_ref_names;
    std::vector<uint64_t> m_ref_lens;
};
}  // namespace mindex
