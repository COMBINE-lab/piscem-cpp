#include <chrono>
#include <iostream>
#include <string>
#include <vector>

#include "../include/reference_index.hpp"
#include "../include/streaming_query.hpp"
#include "../include/CanonicalKmerIterator.hpp"
#include "../include/FastxParser.hpp"
#include "../include/cli11/CLI11.hpp"
#include "../include/spdlog_piscem/spdlog.h"
#include "../include/spdlog_piscem/sinks/stdout_color_sinks.h"

int main(int argc, char** argv) {
    std::string index_filename;
    std::string query_filename;
    bool locate = false;

    CLI::App app{"sshash streaming lookup benchmark"};
    app.add_option("-i,--index", index_filename, "Index prefix")->required();
    app.add_option("-q,--query", query_filename, "Query FASTA/FASTQ file")->required();
    app.add_flag("--locate", locate, "Also perform locate (contig table lookup) for each hit");
    CLI11_PARSE(app, argc, argv);

    spdlog_piscem::drop_all();
    auto logger = spdlog_piscem::create<spdlog_piscem::sinks::stderr_color_sink_mt>("");
    logger->set_pattern("%+");
    spdlog_piscem::set_default_logger(logger);

    spdlog_piscem::info("loading index from {}", index_filename);
    mindex::reference_index ri(index_filename);
    spdlog_piscem::info("index loaded");

    CanonicalKmer::k(ri.k());
    const uint64_t k = ri.k();

    // Load all query sequences into memory
    spdlog_piscem::info("loading queries from {}", query_filename);
    std::vector<std::string> sequences;
    sequences.reserve(300000);

    {
        fastx_parser::ParserConfig pc;
        std::vector<std::string> rfiles{query_filename};
        fastx_parser::FastxParser<fastx_parser::ReadSeq> rparser(pc, rfiles);
        rparser.start();
        auto rg = rparser.getReadGroup();
        while (rparser.refill(rg)) {
            for (auto& record : rg) {
                sequences.push_back(std::string(record.first().seq));
            }
        }
        rparser.stop();
    }

    uint64_t total_kmers = 0;
    for (auto& s : sequences) {
        if (s.size() >= k) total_kmers += s.size() - k + 1;
    }
    spdlog_piscem::info("loaded {} sequences, {} k-mer positions", sequences.size(), total_kmers);

    // Set up streaming query
    piscem::unitig_end_cache_t unitig_end_cache(5000000);
    piscem::streaming_query<true> q(ri.get_dict(), &unitig_end_cache);

    uint64_t found = 0;
    uint64_t num_kmers = 0;

    // Time ONLY the query loop
    spdlog_piscem::info("starting benchmark (locate={})", locate);
    auto t_start = std::chrono::high_resolution_clock::now();

    for (auto& seq : sequences) {
        if (seq.size() < k) continue;

        pufferfish::CanonicalKmerIterator kit(seq);
        pufferfish::CanonicalKmerIterator end;

        while (kit != end) {
            auto proj_hits = ri.query(kit, q);
            num_kmers++;
            if (!proj_hits.empty()) {
                found++;
                if (locate) {
                    for (auto v : proj_hits.refRange) {
                        auto ref_pos_ori = proj_hits.decode_hit(v);
                        (void)ref_pos_ori;
                    }
                }
            }
            ++kit;
        }
    }

    auto t_stop = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_stop - t_start);

    double ns_per_kmer = static_cast<double>(elapsed.count()) / num_kmers;
    uint64_t extensions = q.num_extensions();
    uint64_t searches = q.num_searches();

    std::cout << "==== streaming lookup report:\n";
    std::cout << "num_kmers = " << num_kmers << "\n";
    std::cout << "found_kmers = " << found << " ("
              << (num_kmers > 0 ? static_cast<double>(found) / num_kmers * 100.0 : 0)
              << "%)\n";
    std::cout << "searches = " << searches << "\n";
    std::cout << "extensions = " << extensions << "\n";
    std::cout << "extension_ratio = "
              << (searches > 0 ? static_cast<double>(extensions) / searches : 0) << "\n";
    std::cout << "time_per_kmer = " << ns_per_kmer << " ns\n";
    std::cout << "total_time = " << elapsed.count() / 1e9 << " s\n";

    spdlog_piscem::drop_all();
    return 0;
}
