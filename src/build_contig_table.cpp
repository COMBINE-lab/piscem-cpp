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

    {
        // In the first pass over the cf_seq file we
        // will assign each segment an ID based on the
        // order of its first appearance (rank), and
        // will count how many times each segment occurs.
        std::ifstream ifile(input_filename + ".cf_seq");

        uint64_t refctr = 0;
        bool first = true;
        uint64_t next_rank = 0;
        std::vector<std::string> ref_names;

        while (!ifile.eof()) {
            std::string tok;
            while (ifile >> tok) {
                // this is a new reference
                if (tok.compare(0, hlen, refstr) == 0) {

                    auto sp = tok.find("Sequence:");
                    auto ep = sp+9;
                    std::string refname = tok.substr(ep);
                    ref_names.push_back(refname);

                    if (!first) {
                        ++refctr;
                        std::cerr << "processing reference #" << refctr << " : " << refname << "\n";
                    }
                    first = false;
                } else {  // this should be a segemnt entry

                    if (!((tok.back() == '-') or (tok.back() == '+'))) {
                        std::cerr << "unexpected last character of tiling entry [" << tok.back()
                                  << "\n";
                        std::exit(1);
                    }

                    tok.pop_back();
                    uint64_t id = std::stoul(tok, nullptr, 0);

                    auto rit = id_to_rank.find(id);
                    if (rit == id_to_rank.end()) {
                        id_to_rank[id] = {next_rank, 0, 1};
                        ++next_rank;
                    } else {
                        rit->second.count += 1;
                    }
                }
            }
        }
        
        ser(ref_names);
        // flush to writer
        ser.adapter().flush();
        //s.close();
    }

    std::cerr << "completed first pass over paths.\n";
    std::cerr << "there were " << id_to_rank.size() << " segments\n";

    uint64_t tot_seg_occ = 0;
    for (auto& kv : id_to_rank) { tot_seg_occ += kv.second.count; }

    std::cerr << "there were " << tot_seg_occ << " total segement occurrences\n";

    std::vector<uint64_t> segment_order;
    segment_order.reserve(id_to_rank.size() + 1);
    {
        // now we go over the file defining the segments
        // to determine the length of each segment.
        std::ifstream seg_file(input_filename + ".cf_seg");
        while (!seg_file.eof()) {
            uint64_t seg_id;
            uint32_t seg_len;
            std::string seg;
            while (seg_file >> seg_id >> seg) {
                seg_len = seg.length();
                id_to_rank[seg_id].len = seg_len;
                segment_order.push_back(seg_id);
            }
        }
    }
    std::cerr << "computed all segment lengths.\n";

    std::cerr << "computing cumulative offset vector.\n";

    basic_contig_table bct;
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

        std::cerr << "converting segment counts to offsets.\n";
        // now convert each `count` entry for each contig to
        // the current offset where it's next entry will be written
        for (size_t i = 0; i < segment_order.size(); ++i) {
            auto seg_id = segment_order[i];
            id_to_rank[seg_id].count = contig_offsets[i];
            id_to_rank[seg_id].rank = i;
        }
        // since the contig offset vector is a monotonic sequence
        // it is amenable to Elias-Fano compression, so compress it
        // as such and write it.
        // ef_sequence efo;
        bct.m_ctg_offsets.encode(contig_offsets.begin(), contig_offsets.size());
        // std::string cto_fname = output_filename+"_coff.bin";
        // essentials::save(efo, cto_fname.c_str());
    }

    std::cerr << "second pass over seq file to fill in contig entires.\n";
    {
        // Finally, we'll go over the sequences of segments again
        // and build the final table.
        auto& seg_table = bct.m_ctg_entries;
        seg_table.resize(tot_seg_occ);
        std::ifstream ifile(input_filename + ".cf_seq");

        uint64_t refctr = 0;
        bool first = true;
        uint64_t current_offset = 0;
        std::vector<uint64_t> ref_lens;

        while (!ifile.eof()) {
            std::string tok;
            uint64_t tctr = 0;
            while (ifile >> tok) {
                // this is a new reference
                if (tok.compare(0, hlen, refstr) == 0) {
                    if (!first) {
                        ++refctr;
                        uint64_t len = current_offset + (k - 1);
                        ref_lens.push_back(len);
                    }
                    std::cerr << "processing reference #" << refctr << "\n";
                    first = false;
                    current_offset = 0;
                    tctr = 0;
                } else {  // this should be a segemnt entry

                    bool is_fw = true;
                    if (tok.back() == '-') {
                        is_fw = false;
                    } else if (tok.back() == '+') {
                        is_fw = true;
                    } else {
                        std::cerr << "unexpected last character of tiling entry [" << tok.back()
                                  << "\n";
                        std::exit(1);
                    }
                    tok.pop_back();
                    uint64_t id = std::stoul(tok, nullptr, 0);
                    // get the entry for this segment
                    auto& v = id_to_rank[id];
                    ++tctr;

                    // insert the next entry for this segment
                    // at the index given by v.count.
                    auto entry_idx = v.count;
                    seg_table[entry_idx].update(refctr, current_offset, is_fw);
                    // then we increment entry_idx for next time
                    v.count += 1;
                    // then we increment the current offset
                    current_offset += v.len - (k - 1);
                }
            }
        }
        // for the last reference
        uint64_t len = current_offset + (k - 1);
        ref_lens.push_back(len);
        // serialize this part of the segment table.
        ser(ref_lens);
        // flush to writer
        ser.adapter().flush();
    }
    // close the buffer that we are serializing to
    s.close();

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

int build_contig_table_main(const std::string& input_filename, uint64_t k, const std::string& output_filename) {
    bool success = build_contig_table(input_filename, k, output_filename);
    if (!success) { std::cerr << "failed to build contig table.\n"; }
    return 0;
}
