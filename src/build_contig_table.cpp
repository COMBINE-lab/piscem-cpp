#include <iostream>
#include <vector>
#include "../include/basic_contig_table.hpp"
#include "../external/pthash/external/essentials/include/essentials.hpp"
#include "../include/ef_sequence.hpp"
#include "../external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "../include/util.hpp"
#include "../include/parallel_hashmap/phmap.h"
#include "../include/bitsery/bitsery.h"
#include "../include/bitsery/brief_syntax.h"
#include "../include/bitsery/adapter/stream.h"
#include "../include/bitsery/brief_syntax/vector.h"
#include "../include/bitsery/brief_syntax/string.h"
#include "../include/spdlog/spdlog.h"
#include "../include/json.hpp"

using namespace sshash;
using phmap::flat_hash_map;

struct rank_count {
    uint64_t rank;
    uint32_t len;
    uint64_t count;
};

bool build_contig_table(const std::string& input_filename, uint64_t k,
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
        spdlog::info("computed all segment lengts");
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
                            spdlog::error(
                                "ref_names is empty, but first is false; should not happen!");
                        }
                        auto rn = ref_names.back();
                        uint64_t len = current_offset + (k - 1);
                        ref_lens.push_back(len);
                        max_ref_len = std::max(static_cast<uint64_t>(max_ref_len), len);
                        if (refctr % 10000 == 0) {
                            spdlog::info("finished processing reference #{} : {}, len : {}", refctr,
                                         rn, len);
                        }
                        ++refctr;
                    }

                    std::string refname = tok.substr(ep);
                    ref_names.push_back(refname);

                    current_offset = 0;
                    first = false;
                } else {  // this should be a segment entry

                    if (!((tok.back() == '-') or (tok.back() == '+'))) {
                        spdlog::critical("unexpected last character of tiling entry [{}]",
                                         tok.back());
                        std::exit(1);
                    }

                    tok.pop_back();
                    uint64_t id = std::stoul(tok, nullptr, 0);

                    auto rit = id_to_rank.find(id);
                    if (rit == id_to_rank.end()) {
                        spdlog::critical(
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
        if (ref_names.empty()) {
            spdlog::critical("ref_names is empty, but first is false; should not happen!");
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
            spdlog::info("finished processing reference #{} : {}, len : {}", refctr, rn, len);
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

    spdlog::info("completed first pass over paths.");
    spdlog::info("there were {} segments.", id_to_rank.size());
    spdlog::info("max ref len = {}, requires {} bits.", max_ref_len, ref_len_bits);
    spdlog::info("max refs = {}, requires {} bits.", num_refs, num_ref_bits);

    uint64_t tot_seg_occ = 0;
    for (auto& kv : id_to_rank) { tot_seg_occ += kv.second.count; }

    spdlog::info("there were {} total segment occurrences", tot_seg_occ);
    spdlog::info("computing cumulative offset vector.");

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

        spdlog::info("converting segment coutns to offsets.");
        // now convert each `count` entry for each contig to
        // the current offset where its next entry will be written
        for (size_t i = 0; i < segment_order.size(); ++i) {
            auto seg_id = segment_order[i];
            id_to_rank[seg_id].count = contig_offsets[i];
            if (id_to_rank[seg_id].rank != i) {
                spdlog::critical("expected segment {} to have rank {}, but it had rank {}", seg_id,
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

    spdlog::info("second pass over seq file to fill in contig entries.");
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
                    if (refctr % 10000 == 0) { spdlog::info("processing reference #{}", refctr); }
                    first = false;
                    current_offset = 0;
                } else {  // this should be a segment entry
                    bool is_fw = true;
                    if (tok.back() == '-') {
                        is_fw = false;
                    } else if (tok.back() == '+') {
                        is_fw = true;
                    } else {
                        spdlog::critical("unexpected last character of tiling entry [{}]",
                                         tok.back());
                        std::exit(1);
                    }
                    tok.pop_back();
                    uint64_t id = std::stoul(tok, nullptr, 0);
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
        seg_table_builder.build(bct.m_ctg_entries);
    }

    std::string out_ctab = output_filename + ".ctab";
    essentials::save(bct, out_ctab.c_str());

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
                            const std::string& output_filename) {
    bool success = build_contig_table(input_filename, k, output_filename);
    if (!success) {
        spdlog::critical("failed to build contig table.");
        return 1;
    }
    return 0;
}
