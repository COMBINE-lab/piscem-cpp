#include "../include/CanonicalKmerIterator.hpp"
#include "../include/Kmer.hpp"
#include "../include/hit_searcher.hpp"
#include <iostream>
using namespace std;

// All temp code will be created and tested here

bool get_raw_hits_everykmer(mindex::hit_searcher &hs,
                  std::string &read,
                  sshash::streaming_query_canonical_parsing& qc,
                  bool isLeft,
                  bool verbose) {
    (void) verbose;
    projected_hits phits;
    auto& raw_hits = isLeft ? hs.get_left_hits() : hs.get_right_hits();
    CanonicalKmer::k(hs.get_k());
    pufferfish::CanonicalKmerIterator kit(read), kit_end;
    while(kit != kit_end) {
        phits = hs.get_index()->query(kit, qc);
        raw_hits.push_back(std::make_pair(kit->second, phits));
        kit++;
    }
    return raw_hits.size() != 0;
}

int main() {
    mindex::reference_index ri("/fs/cbcb-lab/rob/students/noor/Atacseq/piscem_analysis/hg38_ind_k23/hg38_ind_k23");
    mindex::hit_searcher hs(&ri);
    sshash::streaming_query_canonical_parsing q(ri.get_dict());
    std::string seq1 = "CCTTCTCCTAATTCCGCAAATGTGAAGGGTAGGGGGACGTTAAGGACCTG";
    get_raw_hits_everykmer(hs, 
                        seq1,
                        q, true, true);
    std::cout << hs.get_left_hits().size() << "\n";
    return 0;
}