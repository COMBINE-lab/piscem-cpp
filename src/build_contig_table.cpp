#include <iostream>
#include <vector>
#include "../include/basic_contig_table.hpp"
#include "../external/pthash/external/essentials/include/essentials.hpp"
#include "../include/ef_sequence.hpp"
#include "../include/equivalence_class_map.hpp"
#include "../external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "../include/util.hpp"
#include "../include/parallel_hashmap/phmap.h"
#include "../include/bitsery/bitsery.h"
#include "../include/bitsery/brief_syntax.h"
#include "../include/bitsery/adapter/stream.h"
#include "../include/bitsery/brief_syntax/vector.h"
#include "../include/bitsery/brief_syntax/string.h"
#include "../include/spdlog_piscem/spdlog.h"
#include "../include/json.hpp"
#include "../external/pthash/include/utils/hasher.hpp"

using namespace sshash;
using phmap::flat_hash_map;

struct rank_count {
    uint64_t rank;
    uint32_t len;
    uint64_t count;
};

struct rank_offset {
  uint32_t rank;
  uint64_t offset;
};

enum class dir_status : uint8_t {
  FW=0, RC=1, BOTH=2
};


namespace std {

  template <>
  struct hash<std::vector<std::tuple<uint32_t, dir_status>>>
  {
    std::size_t operator()(const std::vector<std::tuple<uint32_t, dir_status>>& k) const
    {
       const void* data = reinterpret_cast<const void*>(&k[0]);
       uint64_t len = k.size() * sizeof(std::tuple<uint32_t, dir_status>);
       return pthash::MurmurHash2_64(data, len, 0);
    }
  };

}

bool build_contig_table(const std::string& input_filename, uint64_t k,
                        bool build_eq_table,
                        const std::string& output_filename) {
    flat_hash_map<uint64_t, rank_count> id_to_rank;
    const std::string refstr = "Reference";
    const auto hlen = refstr.length();

    // where we will write the reference info
    std::string out_refinfo = output_filename + ".refinfo";
    std::fstream s{out_refinfo.c_str(), s.binary | s.trunc | s.out};
    bitsery::Serializer<bitsery::OutputBufferedStreamAdapter> ser{s};

    std::vector<uint64_t> segment_order;
    {
        // First, we will pass over the segment file to collect the
        // identifier and length of each segment.
        std::ifstream seg_file(input_filename + ".cf_seg");
        uint64_t idx = 0;
        while (!seg_file.eof()) {
            uint64_t seg_id;
            uint32_t seg_len;
            std::string seg;
            while (seg_file >> seg_id >> seg) {
                segment_order.push_back(seg_id);
                seg_len = seg.length();
                id_to_rank[seg_id] = {idx, seg_len, 0};
                ++idx;
            }
        }
        spdlog_piscem::info("computed all segment lengts");
    }

    size_t num_refs = 0;
    size_t max_ref_len = 0;
    {
        // In the first pass over the cf_seq file we
        // will assign each segment an ID based on the
        // order of its first appearance (rank),
        // will count how many times each segment occurs, and
        // will compute the lengths of all reference
        // sequences.
        std::ifstream ifile(input_filename + ".cf_seq");

        uint64_t refctr = 0;
        bool first = true;
        uint64_t current_offset = 0;
        std::vector<std::string> ref_names;
        std::vector<uint64_t> ref_lens;

        while (!ifile.eof()) {
            std::string tok;
            while (ifile >> tok) {
                // this is a new reference
                if (tok.compare(0, hlen, refstr) == 0) {
                    auto sp = tok.find("Sequence:");
                    auto ep = sp + 9;

                    if (!first) {
                        if (ref_names.empty()) {
                            spdlog_piscem::error(
                                "ref_names is empty, but first is false; should not happen!");
                        }
                        auto rn = ref_names.back();
                        uint64_t len = current_offset + (k - 1);
                        ref_lens.push_back(len);
                        max_ref_len = std::max(static_cast<uint64_t>(max_ref_len), len);
                        if (refctr % 10000 == 0) {
                            spdlog_piscem::info("finished processing reference #{} : {}, len : {}", refctr,
                                         rn, len);
                        }
                        ++refctr;
                    }

                    std::string refname = tok.substr(ep);
                    ref_names.push_back(refname);

                    current_offset = 0;
                    first = false;
                } else {  // this should be a segment entry
                    bool is_n_tile = false;
                    if (!((tok.back() == '-') or (tok.back() == '+'))) {
                        // in this case, the first character must be an 'N'
                        if ( tok.front() == 'N' ) { 
                          is_n_tile = true;
                        } else {
                          spdlog_piscem::critical("Unless a tiling entry is an 'N' entry, it must end with '+' or '-'. "
                                           "Found unexpected last character [{}] of tiling entry.",
                                           tok.back());
                          std::exit(1);
                        }
                    }

                    if (is_n_tile) {
                      // skip the number of 'N's and the overlapping (k-1)-mer
                      tok.erase(0,1); // remove the actual 'N'
                      uint64_t num_ns = std::stoul(tok, nullptr, 10);
                      // if this was the first tile it needs special handling.
                      // specifically, we *shouldn't* skip the k-1 overlap we should 
                      // just skip the leading 'N's.
                      if (current_offset > 0) {
                        current_offset += (k-1);
                      }
                      current_offset += num_ns;
                    } else {
                      tok.pop_back();
                      uint64_t id = std::stoul(tok, nullptr, 10);

                      auto rit = id_to_rank.find(id);
                      if (rit == id_to_rank.end()) {
                        spdlog_piscem::critical(
                            "encountered segment {} that was not found in the id_to_rank "
                            "dictionary!",
                            id);
                        std::exit(1);
                      } else {
                        rit->second.count += 1;
                        // then we increment the current offset
                        current_offset += rit->second.len - (k - 1);
                      }
                    }
                }
            }
        }
        if (ref_names.empty()) {
            spdlog_piscem::critical("ref_names is empty, but first is false; should not happen!");
            std::exit(1);
        }
        auto rn = ref_names.back();
        // for the last reference
        uint64_t len = current_offset + (k - 1);
        ref_lens.push_back(len);
        max_ref_len = std::max(static_cast<uint64_t>(max_ref_len), len);

        // from the json file we will get info about any
        // short references
        {
            using short_refs_t = std::vector<std::pair<std::string, size_t>>;
            nlohmann::json dbg_info;
            std::ifstream json_file(input_filename + ".json");
            json_file >> dbg_info;
            if (dbg_info.contains("short refs")) {
                short_refs_t short_refs_info = dbg_info["short refs"].get<short_refs_t>();
                for (auto& p : short_refs_info) {
                    ref_names.push_back(p.first);
                    ref_lens.push_back(static_cast<uint64_t>(p.second));
                    max_ref_len = std::max(static_cast<uint64_t>(max_ref_len), len);
                }
            }
        }

        if (refctr % 10000 == 0) {
            spdlog_piscem::info("finished processing reference #{} : {}, len : {}", refctr, rn, len);
        }

        num_refs = ref_lens.size();

        {
            ser(ref_names);
            // flush to writer
            ser.adapter().flush();
            // s.close();
        }

        {
            // serialize this part of the segment table.
            ser(ref_lens);
            // flush to writer
            ser.adapter().flush();
        }
    }
    // close the buffer that we are serializing to
    s.close();

    uint64_t ref_len_bits = std::ceil(std::log2(max_ref_len + 1));
    uint64_t num_ref_bits = std::ceil(std::log2(num_refs + 1));
    uint64_t total_ctg_bits = ref_len_bits + num_ref_bits + 1;

    // to get to the ref we shift ref_len_bits + 1 (orientation bit)
    sshash::util::_ref_shift = ref_len_bits + 1;
    sshash::util::_pos_mask = sshash::util::pos_masks[ref_len_bits];

    spdlog_piscem::info("completed first pass over paths.");
    spdlog_piscem::info("there were {} segments.", id_to_rank.size());
    spdlog_piscem::info("max ref len = {}, requires {} bits.", max_ref_len, ref_len_bits);
    spdlog_piscem::info("max refs = {}, requires {} bits.", num_refs, num_ref_bits);

    uint64_t tot_seg_occ = 0;
    for (auto& kv : id_to_rank) { tot_seg_occ += kv.second.count; }

    spdlog_piscem::info("there were {} total segment occurrences", tot_seg_occ);
    spdlog_piscem::info("computing cumulative offset vector.");

    basic_contig_table bct;
    bct.m_ref_len_bits = ref_len_bits;
    {
        // next we compute an offest vector that will point
        // to where, in the concatenated contig occurrence
        // table, the sublist for each contig starts.
        // this looks like:
        // [0, #occ(ctg_0), #occ(ctg_0)+#occ(ctg_1), ...]
        // in other words, it is a cumulative sum, padded with 0
        // at the start.
        std::vector<uint64_t> contig_offsets;
        contig_offsets.reserve(id_to_rank.size() + 1);
        contig_offsets.push_back(0);
        uint64_t total_occ = 0;
        for (auto seg_id : segment_order) {
            total_occ += id_to_rank[seg_id].count;
            contig_offsets.push_back(total_occ);
        }

        spdlog_piscem::info("converting segment coutns to offsets.");
        // now convert each `count` entry for each contig to
        // the current offset where its next entry will be written
        for (size_t i = 0; i < segment_order.size(); ++i) {
            auto seg_id = segment_order[i];
            id_to_rank[seg_id].count = contig_offsets[i];
            if (id_to_rank[seg_id].rank != i) {
                spdlog_piscem::critical("expected segment {} to have rank {}, but it had rank {}", seg_id,
                                 i, id_to_rank[seg_id].rank);
                std::exit(1);
            }
        }
        // since the contig offset vector is a monotonic sequence
        // it is amenable to Elias-Fano compression, so compress it
        // as such and write it.
        // ef_sequence efo;
        bct.m_ctg_offsets.encode(contig_offsets.begin(), contig_offsets.size(),
                                 contig_offsets.back());
        // std::string cto_fname = output_filename+"_coff.bin";
        // essentials::save(efo, cto_fname.c_str());
    }

    spdlog_piscem::info("second pass over seq file to fill in contig entries.");
    {
        // Finally, we'll go over the sequences of segments again
        // and build the final table.
        auto seg_table_builder = pthash::compact_vector::builder(tot_seg_occ, total_ctg_bits);
        // auto& seg_table = bct.m_ctg_entries;
        // seg_table.resize(tot_seg_occ);
        std::ifstream ifile(input_filename + ".cf_seq");

        uint64_t refctr = 0;
        bool first = true;
        uint64_t current_offset = 0;

        while (!ifile.eof()) {
            std::string tok;
            while (ifile >> tok) {
                // this is a new reference
                if (tok.compare(0, hlen, refstr) == 0) {
                    if (!first) { ++refctr; }
                    if (refctr % 10000 == 0) { spdlog_piscem::info("processing reference #{}", refctr); }
                    first = false;
                    current_offset = 0;
                } else {  // this should be a segment entry
                    bool is_n_tile = false;
                    bool is_fw = true;
                    if (tok.back() == '-') {
                        is_fw = false;
                    } else if (tok.back() == '+') {
                        is_fw = true;
                    } else if (tok.front() == 'N') {
                        is_n_tile = true;
                    } else {
                        spdlog_piscem::critical("Unless a tiling entry is an 'N' entry, it must end with '+' or '-'. "
                                         "Found unexpected last character [{}] of tiling entry.",
                                          tok.back());
                        std::exit(1);
                    }
                    
                    if (is_n_tile) {
                      tok.erase(0,1); // remove the leading 'N'
                      uint64_t num_ns = std::stoul(tok, nullptr, 10);
                      // if this was the first tile it needs special handling.
                      // specifically, we *shouldn't* skip the k-1 overlap we should 
                      // just skip the leading 'N's.
                      if (current_offset > 0) {
                        current_offset += (k-1);
                      }
                      current_offset += num_ns;
                    } else {
                      tok.pop_back();
                      uint64_t id = std::stoul(tok, nullptr, 10);
                      // get the entry for this segment
                      auto& v = id_to_rank[id];

                      // insert the next entry for this segment
                      // at the index given by v.count.
                      auto entry_idx = v.count;
                      uint64_t encoded_entry =
                        sshash::util::encode_contig_entry(refctr, current_offset, is_fw);
                      //[entry_idx].update(refctr, current_offset, is_fw);
                      seg_table_builder.set(entry_idx, encoded_entry);
                      // then we increment entry_idx for next time
                      v.count += 1;
                      // then we increment the current offset
                      current_offset += v.len - (k - 1);
                    }
                }
            }
        }
        seg_table_builder.build(bct.m_ctg_entries);
    }

    std::string out_ctab = output_filename + ".ctab";
    essentials::save(bct, out_ctab.c_str());

    if (build_eq_table) {
      // map equivalence class content to id
      flat_hash_map<std::vector<std::tuple<uint32_t, dir_status>>, rank_offset>  ec_id_map;
      std::vector<std::tuple<uint32_t, dir_status>> label;
      uint64_t largest_label = 0;

      // the number of tiles
      size_t num_tiles = id_to_rank.size();

      // the ec table
      equivalence_class_map ect;

      // will hold the total length of concatenated labels of all equivalence 
      // classes
      size_t total_label_length = 0;
      {
        // will hold the starting position in the globally concatenated 
        // list, for the sublist corresponding to each equivalence class
        std::vector<uint64_t> label_list_offsets;

        // holds the corresponding ec id for each contig
        std::vector<uint32_t> tile_ec_ids;
        tile_ec_ids.reserve(num_tiles);

        for (size_t tile_idx = 0; tile_idx < num_tiles; ++tile_idx) {
          label.clear();
          uint32_t prev_tid = 0;
          dir_status prev_dir = dir_status::FW;
          bool first = true;
          sshash::util::contig_span ctg_entry_span = bct.contig_entries(tile_idx);
          for (auto ce : ctg_entry_span) {
            // note, these sshash::util:: functions should be safe to call because 
            // we have set the relevant static members above.
            uint32_t tid = sshash::util::transcript_id(ce);
            dir_status dir = sshash::util::orientation(ce) ? dir_status::FW : dir_status::RC;
            largest_label = std::max(static_cast<uint64_t>(tid), largest_label);

            // skip adjacent duplicates (we can still get dups because of orientation 
            // switching), but if duplicate target / ori pairs are adjacent, don't 
            // add them to avoid the vector growing unnecessarily.
            if (first or tid != prev_tid) {
              label.push_back({tid, dir});
            } else if (tid == prev_tid and dir != prev_dir) {
              // if dir != prev_dir, then 
              // dir == FW and prev_dir == RC
              // or dir == RC and prev_dir == FW
              // or prev_dir == BOTH and dir == FW | RC 
              // in any such case, the right thing to do is 
              // to set (or keep) the dir as BOTH.
              std::get<1>(label.back()) = dir_status::BOTH;
            }

            prev_tid = tid;
            prev_dir = dir;
            first = false;
          }
          // remove any duplicate entries --- shouldn't be any!
          std::sort(label.begin(), label.end());
          auto last_valid_it = std::unique(label.begin(), label.end());
          if (last_valid_it != label.end()) {
            std::cerr << "sort | unique should not be necessary on ec label!\n";
            label.erase(last_valid_it, label.end());
          }

          // see if we know about this class already
          uint32_t next_ec_rank = static_cast<uint32_t>(ec_id_map.size());
          auto ec_it = ec_id_map.insert({label, {next_ec_rank, total_label_length}});
          // if this is newly inserted
          if (ec_it.second) {
            label_list_offsets.push_back(total_label_length);
            total_label_length += label.size();
          }
          // may have been set / updated in the if above
          tile_ec_ids.push_back((ec_it.first)->second.rank);

        } // created the hash map, now pack it into an efficient structure 
          // last label list ending position
        label_list_offsets.push_back(total_label_length);
        // since the offset vector is a monotonic sequence
        // it is amenable to Elias-Fano compression, so compress it
        // as such and write it.
        ect.m_label_list_offsets.encode(label_list_offsets.begin(), label_list_offsets.size(),
            label_list_offsets.back());
        
        uint64_t tile_id_width = std::ceil(std::log2(ec_id_map.size() + 1));
        pthash::compact_vector::builder ec_id_builder(tile_ec_ids.begin(), tile_ec_ids.size(), tile_id_width);
        ec_id_builder.build(ect.m_tile_ec_ids);
      } // end scope to free labal_list_offsets 

      // the +2 is for the orientation bits;
      uint64_t label_width = std::ceil(std::log2(largest_label + 1)) + 2;
      pthash::compact_vector::builder label_builder(total_label_length, label_width);

      // now, iterate over the equivalence class map and pack the 
      // label information into the concatenated vector
      for (auto& kv : ec_id_map) {
        // where we start writing 
        uint64_t offset_pos = kv.second.offset;
        for (auto ref_dir : kv.first) {
          auto txpid = std::get<0>(ref_dir);
          auto dir = std::get<1>(ref_dir);
          // pack this into the compact_vector 
          uint64_t packed_entry = (txpid << 2);
          switch (dir) {
            case dir_status::FW:
              break;
            case dir_status::RC:
              packed_entry |= 0x1;
              break;
            case dir_status::BOTH:
              packed_entry |= 0x2;
          }

          label_builder.set(offset_pos, packed_entry);
          ++offset_pos;
        }
      }

      label_builder.build(ect.m_label_entries);

      std::string out_ectab = output_filename + ".ectab";
      essentials::save(ect, out_ectab.c_str());
    }

    /*
    std::cerr << "verifying contig table invariants.\n";

    for (size_t i = 0; i < bct.m_ctg_offsets.size() - 1; ++i) {
        auto curr_start = bct.m_ctg_offsets.access(i);
        auto curr_end = bct.m_ctg_offsets.access(i + 1);

        uint32_t prior_ref = 0;
        uint32_t prior_pos = 0;
        for (size_t j = curr_start; j < curr_end; ++j) {
            auto p = bct.m_ctg_entries[j];

            if (p.transcript_id() < prior_ref) {
                std::cerr << "m_ctg_offsets[" << i << "] = " << curr_start << "\n";
                std::cerr << "m_ctg_offsets[" << i + 1 << "] = " << curr_end << "\n";
                std::cerr << "entry " << j << " had ref id = " << p.transcript_id() << " < "
                          << prior_ref << "\n";
            } else if ((p.transcript_id() == prior_ref) && (p.pos() < prior_pos)) {
                std::cerr << "m_ctg_offsets[" << i << "] = " << curr_start << "\n";
                std::cerr << "m_ctg_offsets[" << i + 1 << "] = " << curr_end << "\n";
                std::cerr << "entry " << j << " had ref id = " << prior_ref << " and pos "
                          << p.pos() << " > " << prior_pos << "\n";
            }
            prior_ref = p.transcript_id();
            prior_pos = p.pos();
        }
    }
    */
    return true;
}

int build_contig_table_main(const std::string& input_filename, uint64_t k,
                            bool build_eq_table,
                            const std::string& output_filename) {
    bool success = build_contig_table(input_filename, k, build_eq_table, output_filename);
    if (!success) {
        spdlog_piscem::critical("failed to build contig table.");
        return 1;
    }
    return 0;
}
