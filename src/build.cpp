#include <iostream>
#include <memory>
#include <thread>

#include "../external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "../include/dictionary.hpp"
#include "../include/spdlog_piscem/spdlog.h"
#include "../include/spdlog_piscem/sinks/stdout_color_sinks.h"
#include "../include/cli11/CLI11.hpp"
#include "bench_utils.hpp"
#include "check_utils.hpp"
#include "build_contig_table.cpp"
#include "bench_utils.hpp"
#include "check_utils.hpp"

using namespace sshash;

#ifdef __cplusplus
extern "C" {
#endif
int run_build(int argc, char** argv);
#ifdef __cplusplus
}
#endif

int run_build(int argc, char** argv) {
    constexpr uint32_t min_threads = 1;
    constexpr uint32_t target_threads = 16;
    uint32_t default_num_threads = std::max(
        min_threads,
        std::min(target_threads, static_cast<uint32_t>(std::thread::hardware_concurrency())));

    bool quiet = false;
    bool build_ec_table = false;
    bool check = false;
    bool bench = false;

    std::string input_files_basename;
    std::string output_filename;
    std::string tmp_dirname;
    build_configuration build_config;

    CLI::App app{"piscem index builder"};
    /* Required arguments. */
    app.add_option(
           "-i,--input", input_files_basename,
           "Must be the basename of input cuttlefish files (expected suffixes are .cf_seq and "
           ".cf_seg, possibly ending with '.gz'.)")
        ->required();
    app.add_option("-k,--klen", build_config.k,
                   "K-mer length (must be <= " + std::to_string(constants::max_k) + ")")
        ->required();
    app.add_option("-m,--minimizer-len", build_config.m, "Minimizer length (must be < k).")->required();
    app.add_option("-o,--output", output_filename,
                   "path to output file prefix where the data strucutre will be serialized.")
        ->required();

    /* optional arguments */
    app.add_flag("--quiet", quiet, "Only write errors or critical messages to the log");
    app.add_option("-s,--seed", build_config.seed, "Seed for construction")
        ->default_val(constants::seed);
    app.add_option("-l,--load", build_config.l,
                   "A (integer) constant that controls the space/time trade-off of the dictionary. "
                   "A reasonable values lies between 2 and 12")
        ->default_val(constants::min_l);
    app.add_option(
           "-c,--cscale", build_config.c,
           "A (floating point) constant that trades construction speed for space effectiveness "
           "of minimal perfect hashing. "
           "A reasonable value lies between 3.0 and 10.0")
        ->default_val(constants::c);

    CLI::Option* tmpdir_opt =
        app.add_option("-d,--tempdir", tmp_dirname,
                       "Temporary directory used for construction in external memory.")
            ->default_val(constants::default_tmp_dirname);

    app.add_option("-t,--num-threads", build_config.num_threads,
                   "Number of threads to use for hash construction (much of the other index "
                   "building is currently single-threaded.")
        ->default_val(default_num_threads);
    app.add_flag("--canonical-parsing", build_config.canonical_parsing,
                 "Canonical parsing of k-mers. This option changes the parsing and results in a "
                 "trade-off between index space and lookup time.");
    app.add_flag("--build-ec-table", build_ec_table,
                 "build orientation-aware equivalence class table an include it in the index.");
    app.add_flag("--weighted", build_config.weighted,
                 "Also store the weights in compressed format.");
    app.add_flag("--verbose", build_config.verbose, "Verbose output during construction.");
    app.add_flag("--check", check, "Check correctness after construction.");
    app.add_flag("--bench", bench, "Run benchmark after construction.");

    CLI11_PARSE(app, argc, argv);

    spdlog_piscem::drop_all();
    auto logger = spdlog_piscem::create<spdlog_piscem::sinks::stderr_color_sink_mt>("");
    logger->set_pattern("%+");
    spdlog_piscem::set_default_logger(logger);

    if (quiet) {
      spdlog_piscem::set_level(spdlog_piscem::level::warn);
    }

    // make sure the number of requested threads is OK
    if (build_config.num_threads == 0) {
        spdlog_piscem::warn("specified 0 threads, defaulting to 1");
        build_config.num_threads = 1;
    }
    uint64_t max_num_threads = std::thread::hardware_concurrency();
    if (build_config.num_threads > max_num_threads) {
        build_config.num_threads = max_num_threads;
        spdlog_piscem::warn("too many threads specified, defaulting to {}", build_config.num_threads);
    }

    // if it was passed in
    if (!tmpdir_opt->empty()) {
        build_config.tmp_dirname = tmp_dirname;
        essentials::create_directory(build_config.tmp_dirname);
    }
    // if (!quiet) { build_config.print(); }

    {
        // make this scope here and put dict inside of it to
        // ensure it goes out of scope before we build the
        // contig table
        auto input_seq = input_files_basename + ".cf_seg";
        dictionary dict;
        dict.build(input_seq, build_config);
        assert(dict.k() == build_config.k);
        auto output_seqidx = output_filename + ".sshash";
        spdlog_piscem::info("saving data structure to disk...");
        essentials::save(dict, output_seqidx.c_str());
        spdlog_piscem::info("DONE");

        if (check) {
            check_correctness_lookup_access(dict, input_seq);
            if (build_config.weighted) check_correctness_weights(dict, input_seq);
            check_correctness_iterator(dict);
        }
        if (bench) {
            perf_test_lookup_access(dict);
            if (dict.weighted()) perf_test_lookup_weight(dict);
            perf_test_iterator(dict);
        }
    }

    // now build the contig table
    bool ctab_ok =
        build_contig_table_main(input_files_basename, build_config.k, build_ec_table, output_filename);
    spdlog_piscem::drop_all();
    return ctab_ok;
}
