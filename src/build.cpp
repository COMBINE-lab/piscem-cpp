#include <iostream>

#include "../external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "../include/dictionary.hpp"
#include "../include/builder/build.cpp"
#include "bench_utils.hpp"
#include "check_utils.hpp"

using namespace sshash;

int main(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);

    /* mandatory arguments */
    parser.add("input_filename",
               "Must be a FASTA file (.fa/fasta extension) compressed with gzip (.gz) or not:\n"
               "\t- without duplicate nor invalid kmers\n"
               "\t- one DNA sequence per line.\n"
               "\tFor example, it could be the de Bruijn graph topology output by BCALM.");
    parser.add("k", "K-mer length (must be <= " + std::to_string(constants::max_k) + ").");
    parser.add("m", "Minimizer length (must be < k).");

    /* optional arguments */
    parser.add("seed",
               "Seed for construction (default is " + std::to_string(constants::seed) + ").", "-s",
               false);
    parser.add("l",
               "A (integer) constant that controls the space/time trade-off of the dictionary. "
               "A reasonable values lies between 2 and 12 (default is " +
                   std::to_string(constants::min_l) + ").",
               "-l", false);
    parser.add("c",
               "A (floating point) constant that trades construction speed for space effectiveness "
               "of minimal perfect hashing. "
               "A reasonable value lies between 3.0 and 10.0 (default is " +
                   std::to_string(constants::c) + ").",
               "-c", false);
    parser.add("canonical_parsing",
               "Canonical parsing of k-mers. This option changes the parsing and results in a "
               "trade-off between index space and lookup time.",
               "--canonical-parsing", true);
    parser.add("store_abundances", "Also store the abundances in compressed format.",
               "--abundances", true);
    parser.add("output_filename", "Output file name where the data structure will be serialized.",
               "-o", false);
    parser.add("check", "Check correctness after construction.", "--check", true);
    parser.add("bench", "Run benchmark after construction.", "--bench", true);
    parser.add("verbose", "Verbose output during construction.", "--verbose", true);

    if (!parser.parse()) return 1;

    auto input_filename = parser.get<std::string>("input_filename");
    auto k = parser.get<uint64_t>("k");
    auto m = parser.get<uint64_t>("m");

    dictionary dict;

    build_configuration build_config;
    build_config.k = k;
    build_config.m = m;

    if (parser.parsed("seed")) build_config.seed = parser.get<uint64_t>("seed");
    if (parser.parsed("l")) build_config.l = parser.get<double>("l");
    if (parser.parsed("c")) build_config.c = parser.get<double>("c");
    build_config.canonical_parsing = parser.get<bool>("canonical_parsing");
    build_config.store_abundances = parser.get<bool>("store_abundances");
    build_config.verbose = parser.get<bool>("verbose");
    build_config.print();

    dict.build(input_filename, build_config);
    assert(dict.k() == k);

    bool check = parser.get<bool>("check");
    if (check) {
        check_correctness_lookup_access(dict, input_filename);
        if (build_config.store_abundances) check_correctness_abundances(dict, input_filename);
        check_correctness_iterator(dict);
    }
    bool bench = parser.get<bool>("bench");
    if (bench) {
        perf_test_lookup_access(dict);
        perf_test_iterator(dict);
    }
    if (parser.parsed("output_filename")) {
        auto output_filename = parser.get<std::string>("output_filename");
        essentials::logger("saving data structure to disk...");
        essentials::save(dict, output_filename.c_str());
        essentials::logger("DONE");
    }

    return 0;
}
