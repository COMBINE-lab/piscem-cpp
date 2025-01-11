#pragma once

#include <fstream>
#include <optional>

#include "../external/pthash/external/essentials/include/essentials.hpp"
#include "../external/sshash/include/dictionary.hpp"
#include "../external/sshash/include/util.hpp"
#include "../include/streaming_query.hpp"
#include "CanonicalKmerIterator.hpp"
#include "basic_contig_table.hpp"
#include "bit_vector_iterator.hpp"
#include "bitsery/adapter/stream.h"
#include "bitsery/bitsery.h"
#include "bitsery/brief_syntax/string.h"
#include "bitsery/brief_syntax/vector.h"
#include "equivalence_class_map.hpp"
#include "ghc/filesystem.hpp"
#include "json.hpp"
#include "projected_hits.hpp"
#include "ref_sig_info.hpp"
#include "spdlog_piscem/spdlog.h"
#include "util.hpp"

namespace mindex {
class reference_index {
public:
  reference_index(const std::string &basename,
                  bool attempt_load_ec_map = false) {
    spdlog_piscem::info("loading index from {}", basename);
    std::string sigfile_name = basename + ".sigs.json";

    sig_info = ref_sig_info_t::from_path(sigfile_name);

    std::string dict_name = basename + ".sshash";

    spdlog_piscem::info("sshash loaded");
    essentials::load(m_dict, dict_name.c_str());
    spdlog_piscem::info("ctab loaded");
    std::string ctg_name = basename + ".ctab";
    essentials::load(m_bct, ctg_name.c_str());

    if (attempt_load_ec_map) {
      std::string ectab_name = basename + ".ectab";
      if (ghc::filesystem::exists(ectab_name)) {
        m_has_ec_tab = true;
        essentials::load(m_ec_tab, ectab_name.c_str());
      } else {
        spdlog_piscem::warn("user requested an option that required loading "
                            "the ec map, but that was "
                            "not built for this index. The ec map will not be "
                            "loaded, and any feature "
                            "requiring it "
                            "will be disabled.");
        m_has_ec_tab = false;
      }
    }

    // based on the number of bits used to encode reference positions
    // read from the file, set the shift we have to perform on a
    // contig table entry to read off the reference id (= m_ref_len_bits + 1)
    // where the +1 is for the orientation bit.
    sshash::util::PiscemIndexUtils::ref_shift(m_bct.m_ref_len_bits + 1);
    // based on the value of m_ref_len_bits, select the appropriate mask to use
    // when decoding a reference position from a contig table entry.
    sshash::util::PiscemIndexUtils::pos_mask(
      sshash::util::pos_masks[m_bct.m_ref_len_bits]);

    std::string ref_info = basename + ".refinfo";

    std::fstream s{ref_info.c_str(), s.binary | s.in};
    auto state = bitsery::quickDeserialization<bitsery::InputStreamAdapter>(
      s, m_ref_names);
    state =
      bitsery::quickDeserialization<bitsery::InputStreamAdapter>(s, m_ref_lens);
    spdlog_piscem::info("done loading index");
  }

  projected_hits query(pufferfish::CanonicalKmerIterator &kmit,
                       piscem::streaming_query &q) {
    auto qres = q.query_lookup(kmit, m_bct.m_ctg_offsets, m_bct.m_ctg_entries);
    constexpr uint64_t invalid_u64 = std::numeric_limits<uint64_t>::max();
    constexpr uint32_t invalid_u32 = std::numeric_limits<uint32_t>::max();

    if (q.is_present()) {
      const auto kval = m_dict.k();
      qres.contig_size += kval - 1;
      sshash::util::contig_span s = q.contig_span();
      uint32_t contig_id = (qres.contig_id > invalid_u32)
                             ? invalid_u32
                             : static_cast<uint32_t>(qres.contig_id);
      uint32_t contig_offset =
        (qres.kmer_id_in_contig > invalid_u32)
          ? invalid_u32
          : static_cast<uint32_t>(qres.kmer_id_in_contig);
      uint32_t contig_length = (qres.contig_size > invalid_u32)
                                 ? invalid_u32
                                 : static_cast<uint32_t>(qres.contig_size);

      bool is_forward =
        (qres.kmer_orientation == sshash::constants::forward_orientation);

      // because the query gives us a global
      // ID and not a global offset, we have to
      // adjust it here.
      uint64_t global_offset = qres.kmer_id + (contig_id * (kval - 1));
      return projected_hits{
        contig_id,
        contig_offset,
        is_forward,
        contig_length,
        global_offset, // qres.kmer_id (<- id, not contig offset!)
        static_cast<uint32_t>(kval),
        s};
    } else {
      return {invalid_u32, invalid_u32, false,
              invalid_u32, invalid_u64, static_cast<uint32_t>(m_dict.k()),
              {}};
    }
  }

  uint64_t k() const { return m_dict.k(); }
  const sshash::dictionary *get_dict() const { return &m_dict; }
  pthash::bit_vector const &contigs() { return m_dict.strings(); }
  const std::string &ref_name(size_t i) const { return m_ref_names[i]; }
  uint64_t ref_len(size_t i) const { return m_ref_lens[i]; }
  uint64_t num_refs() const { return m_ref_names.size(); }
  const sshash::basic_contig_table &get_contig_table() const { return m_bct; }

  bool has_ec_table() const { return m_has_ec_tab; }
  const sshash::equivalence_class_map &get_ec_table() { return m_ec_tab; }

  std::optional<ref_sig_info_t> ref_sig_info() const { return sig_info; }

private:
  sshash::dictionary m_dict;
  sshash::basic_contig_table m_bct;
  sshash::equivalence_class_map m_ec_tab;
  std::vector<std::string> m_ref_names;
  std::vector<uint64_t> m_ref_lens;
  std::optional<ref_sig_info_t> sig_info;
  // will be set to true if we have & load
  // and equivalence class table.
  bool m_has_ec_tab{false};
};
} // namespace mindex
