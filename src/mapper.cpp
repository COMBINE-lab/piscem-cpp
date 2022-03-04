
#include "../external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "../include/reference_index.hpp"
#include "../include/kseq++.hpp"
#include "../include/CanonicalKmerIterator.hpp"
//#include "../include/query/contig_info_query_canonical_parsing.cpp"
#include "../include/projected_hits.hpp"
#include "hit_searcher.cpp"
#include "zlib.h"
#include <iostream>

using namespace klibpp;

void do_map(mindex::reference_index& ri, const std::string& reads_filename) {
    CanonicalKmer::k(ri.k());

    KSeq record;
    gzFile fp = gzopen(reads_filename.c_str(), "r");
    auto ks = make_kstream(fp, gzread, mode::in);

    sshash::contig_info_query_canonical_parsing q(ri.get_dict());
    mindex::hit_searcher hs(&ri);
    uint64_t read_num = 0;
    while (ks >> record) {
        ++read_num;
        std::cout << "readnum : " << read_num << "\n";
        q.start();
        hs.clear();
        bool had_left_hit = hs.get_raw_hits_sketch(record.seq, q, true, false);
        auto& left_hits = hs.get_left_hits();

        if (!left_hits.empty()) {
            for (auto& lh : left_hits) {
                auto& ph = lh.second;
                std::cout << "qp : " << lh.first 
                          << ", hit : [ ctg:" << ph.contigIdx_ << " @ " 
                          << ph.contigPos_ << " {" << (ph.contigOrientation_ ? "fw" : "rc") << "}]\n";
            }
        }

    }
}

int main(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);

    /* mandatory arguments */
    parser.add("input_filename",
               "input index prefix.");
    parser.add("reads",
               "read filename.");

    if (!parser.parse()) return 1;

    auto input_filename = parser.get<std::string>("input_filename");
    auto read_filename = parser.get<std::string>("reads");

    mindex::reference_index ri(input_filename);
    
    do_map(ri, read_filename);
    return 0;
}