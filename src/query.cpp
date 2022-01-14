#include <iostream>

#include "../external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "../include/dictionary.hpp"
#include "../include/query/membership_query.hpp"
#include "../include/print_info.cpp"

using namespace sshash;

int main(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "Must be a file generated with src/build.cpp");
    parser.add("query_filename",
               "Must be a FASTA/FASTQ file (.fa/fasta or .fq/fastq extension) compressed with gzip "
               "or not.");
    parser.add("multiline",
               "Use this option if more the one DNA line must be parsed after each header."
               " Only valid for FASTA files (not FASTQ).",
               "--multiline", true);
    parser.add("print_index_info", "Print index information.", "--print-index-info", true);
    if (!parser.parse()) return 1;

    auto index_filename = parser.get<std::string>("index_filename");
    auto query_filename = parser.get<std::string>("query_filename");

    dictionary dict;
    essentials::logger("loading index from file '" + index_filename + "'...");
    uint64_t num_bytes_read = essentials::load(dict, index_filename.c_str());
    std::cout << "index size: " << essentials::convert(num_bytes_read, essentials::MB) << " [MB] ("
              << (num_bytes_read * 8.0) / dict.size() << " [bits/kmer])" << std::endl;
    if (parser.get<bool>("print_index_info")) dict.print_info();

    bool multiline = parser.get<bool>("multiline");

    essentials::logger("performing queries from file '" + query_filename + "'...");
    essentials::timer<std::chrono::high_resolution_clock, std::chrono::microseconds> t;
    t.start();
    auto result = dict.membership_query_from_file(query_filename, multiline);
    t.stop();
    essentials::logger("DONE");

    std::cout << "==== query report:\n";
    std::cout << "num_kmers = " << result.num_kmers << std::endl;
    std::cout << "num_valid_kmers = " << result.num_valid_kmers << " ("
              << (result.num_valid_kmers * 100.0) / result.num_kmers << "% of kmers)" << std::endl;
    std::cout << "num_positive_kmers = " << result.num_positive_kmers << " ("
              << (result.num_positive_kmers * 100.0) / result.num_valid_kmers << "% of valid kmers)"
              << std::endl;
    std::cout << "num_searches = " << result.num_searches << "/" << result.num_positive_kmers
              << " (" << (result.num_searches * 100.0) / result.num_positive_kmers << "%)"
              << std::endl;
    std::cout << "num_extensions = " << result.num_extensions << "/" << result.num_positive_kmers
              << " (" << (result.num_extensions * 100.0) / result.num_positive_kmers << "%)"
              << std::endl;
    std::cout << "elapsed = " << t.elapsed() / 1000 << " millisec / ";
    std::cout << t.elapsed() / 1000000 << " sec / ";
    std::cout << t.elapsed() / 1000000 / 60 << " min / ";
    std::cout << (t.elapsed() * 1000) / result.num_kmers << " ns/kmer" << std::endl;

    return 0;
}
