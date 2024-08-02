#include <iostream>

#include "../external/sshash/external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "../external/sshash/include/dictionary.hpp"
#include "../external/sshash/include/query/membership_query.hpp"
#include "../include/print_info.cpp"

using namespace sshash;

int main(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "Must be a file generated with src/build.cpp");
    if (!parser.parse()) return 1;

    auto index_filename = parser.get<std::string>("index_filename");
    dictionary dict;
    uint64_t num_bytes_read = essentials::load(dict, index_filename.c_str());
    size_t np = dict.get_buckets().pieces.size(); 
    for (size_t i = 0; i < np - 1; ++i) {
        auto c = dict.get_buckets().pieces.access(i);
        auto n = dict.get_buckets().pieces.access(i+1);

        std::cout << (n-c) << "\n";
    }
}