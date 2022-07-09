
#include "../external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "../include/reference_index.hpp"
#include "../include/kseq++.hpp"
#include "../include/CanonicalKmerIterator.hpp"
//#include "../include/query/contig_info_query_canonical_parsing.cpp"
#include "../include/query/streaming_query_canonical_parsing.hpp"
#include "../include/projected_hits.hpp"
#include "zlib.h"

using namespace klibpp;

void check_index(mindex::reference_index& ri, const std::string& ref_fname ){

    CanonicalKmer::k(ri.k());

    KSeq record;
    gzFile fp = gzopen(ref_fname.c_str(), "r");
    auto ks = make_kstream(fp, gzread, mode::in);

    pufferfish::CanonicalKmerIterator kend;
    sshash::streaming_query_canonical_parsing q(ri.get_dict());
    uint64_t refnum = 0;
    uint64_t global_idx = 0;

    while (ks >> record) {
        q.start();
        pufferfish::CanonicalKmerIterator kit(record.seq);
        while (kit != kend) {
            auto ref_hits = ri.query(kit, q);
            bool found_ref_pos = false;

            if ((global_idx % 1000000 == 0) and global_idx > 0) {
                std::cerr << "queried " << global_idx << " k-mers\n";
            }
            ++global_idx;
            if (ref_hits.empty()) {
                std::cerr << "querying reference : " << refnum << " at position " << kit->second << "\n";
                std::cerr << "didn't find k-mer " << kit->first.fwMer().toStr() << " at all in the index!\n";
            }
            
            for (auto h : ref_hits.refRange) {
                auto txp_id = sshash::util::transcript_id(h);
                auto rp = ref_hits.decode_hit(h);
                if ((txp_id == refnum) and 
                    (static_cast<int>(rp.pos) == kit->second) and 
                    (rp.isFW == true)) {
                    found_ref_pos = true;
                    break;
                }
            }
            if (!found_ref_pos) {
                std::cerr << "querying reference : " << refnum << " at position " << kit->second << "\n";
                std::cerr << "did NOT find reference position in hit list!\n";
                for (auto h : ref_hits.refRange) {

                    auto rp = ref_hits.decode_hit(h);
                    
                    std::cerr << "\t[ kmer: " << kit->first.fwMer().toStr()
                              //<< ", rawpos: " << h.pos_
                              << ", contig: " << ref_hits.contig_id() 
                              << ", clen: " << ref_hits.contigLen_ 
                              << ", cstart: " << sshash::util::pos(h) 
                              << ", coff: " << ref_hits.contigPos_
                              << ", cori: " << (ref_hits.contigOrientation_ ? "fw" : "rc")
                              << ", ref: " << sshash::util::transcript_id(h) 
                              << ", pos: " << rp.pos 
                              << ", ori: " << (rp.isFW ? "fw" : "rc") << "]\n";
                    std::exit(1);
                }
            }
            ++kit;
        }
        ++refnum;
    }

}

int main(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);

    /* mandatory arguments */
    parser.add("input_filename",
               "input index prefix.");
    parser.add("ref_filename",
               "reference filename.");

    if (!parser.parse()) return 1;

    auto input_filename = parser.get<std::string>("input_filename");
    auto ref_filename = parser.get<std::string>("ref_filename");

    mindex::reference_index ri(input_filename);

    check_index(ri, ref_filename);
    return 0;
}
