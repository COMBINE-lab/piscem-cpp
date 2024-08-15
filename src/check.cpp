#include <iostream>

#include "../external/sshash/external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "../include/reference_index.hpp"
#include "../include/streaming_query.hpp"
#include "../external/sshash/src/common.hpp"
#include "../external/sshash/src/bench_utils.hpp"
#include "../external/sshash/src/check_utils.hpp"
#include "../include/FastxParser.hpp"
#include "../include/CanonicalKmerIterator.hpp"

using namespace sshash;

int main(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "Must be a file generated with src/build.cpp", "-i", true);
    parser.add("ref_filename", "Reference from which the index was built", "-r", true);
    if (!parser.parse()) return 1;

    auto index_filename = parser.get<std::string>("index_filename");
    auto ref_filename = parser.get<std::string>("ref_filename");

  	mindex::reference_index ri(index_filename);
    piscem::streaming_query q(ri.get_dict());

    // set the canonical k-mer size globally
    CanonicalKmer::k(ri.k());

    fastx_parser::FastxParser<klibpp::KSeq> rparser({ref_filename}, 1); 
    rparser.start();
    auto rg = rparser.getReadGroup();

    while (rparser.refill(rg)) {
        // Here, rg will contain a chunk of references 
        // we can process.
        for (auto& record : rg) {

          pufferfish::CanonicalKmerIterator end;
          pufferfish::CanonicalKmerIterator kit(record.seq);
          while (kit != end) {
            //auto km = kit->first;
            //auto pos = kit->second;

            bool found = false;
            auto& read_pos = kit->second;
            auto proj_hits = ri.query(kit, q);
            if (!proj_hits.empty()) {
              auto& refs = proj_hits.refRange;

              for (auto v : refs) {
                const auto& ref_pos_ori = proj_hits.decode_hit(v);
                uint32_t tid = sshash::util::transcript_id(v);
                auto& ref_name = ri.ref_name(tid);
                int32_t pos = static_cast<int32_t>(ref_pos_ori.pos);
                //bool ori = ref_pos_ori.isFW;
                if ((pos == read_pos) and (record.name == ref_name)) {
                  found = true;
                  break;
                }
              }
            }

            if (!found) {
              std::cerr << "couldn't find k-mer " << kit->first.to_str() << " at position " << read_pos << "\n";
            }

            ++kit;
          }
          
        }
    }


    rparser.stop();
    //check_dictionary(dict);


    return 0;
}
