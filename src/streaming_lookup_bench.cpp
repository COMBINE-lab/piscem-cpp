#include <chrono>
#include <iostream>
#include <string>
#include <vector>

#include "../include/reference_index.hpp"
#include "../include/streaming_query.hpp"
#include "../include/lean_streaming_query.hpp"
#include "../include/CanonicalKmerIterator.hpp"
#include "../include/FastxParser.hpp"
#include "../include/cli11/CLI11.hpp"
#include "../include/spdlog_piscem/spdlog.h"
#include "../include/spdlog_piscem/sinks/stdout_color_sinks.h"

int main(int argc, char** argv) {
    std::string index_filename;
    std::string query_filename;
    bool locate = false;
    bool sshash_native = false;
    bool lean = false;
    bool point_lookup = false;

    bool validate = false;

    CLI::App app{"sshash streaming lookup benchmark"};
    app.add_option("-i,--index", index_filename, "Index prefix")->required();
    app.add_option("-q,--query", query_filename, "Query FASTA/FASTQ file")->required();
    app.add_flag("--locate", locate, "Also perform locate (contig table lookup) for each hit");
    app.add_flag("--sshash-native", sshash_native,
                 "Use sshash's built-in streaming query instead of piscem-cpp's wrapper");
    app.add_flag("--lean", lean,
                 "Use the lean streaming query (sshash engine + contig table lookup)");
    app.add_flag("--point", point_lookup, "Non-streaming point lookup (independent per-kmer queries)");
    app.add_flag("--validate", validate,
                 "Validate lean iterator kmer words match CanonicalKmerIterator");
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

    if (validate) {
        spdlog_piscem::info("validating lean iterator kmer words against CanonicalKmerIterator...");
        piscem::lean_read_iterator lit(ri.get_dict(), ri.get_contig_table());
        uint64_t checked = 0, mismatches = 0;
        uint64_t max_check = std::numeric_limits<uint64_t>::max();

        for (auto& seq : sequences) {
            if (checked >= max_check) break;
            if (seq.size() < k) continue;

            pufferfish::CanonicalKmerIterator kit(seq);
            pufferfish::CanonicalKmerIterator end;
            lit.start(seq.data(), static_cast<int32_t>(seq.size()));

            while (kit != end && !lit.is_exhausted() && checked < max_check) {
                if (kit->second != lit.pos()) {
                    std::cerr << "POSITION MISMATCH: kit=" << kit->second
                              << " lit=" << lit.pos() << "\n";
                    mismatches++;
                    break;
                }

                uint64_t kit_fw = kit->first.fwWord();
                uint64_t kit_rc = kit->first.rcWord();
                uint64_t lit_fw = lit.fw_word();
                uint64_t lit_rc = lit.rc_word();

                if (kit_fw != lit_fw || kit_rc != lit_rc) {
                    if (mismatches < 10) {
                        std::cerr << "KMER MISMATCH at pos " << kit->second
                                  << " in seq of len " << seq.size() << ":\n"
                                  << "  kit fw=0x" << std::hex << kit_fw
                                  << " rc=0x" << kit_rc << std::dec << "\n"
                                  << "  lit fw=0x" << std::hex << lit_fw
                                  << " rc=0x" << lit_rc << std::dec << "\n";
                    }
                    mismatches++;
                }
                checked++;
                ++kit;
                ++lit;
            }
        }

        std::cout << "Validated " << checked << " kmer positions, "
                  << mismatches << " mismatches.\n";
        if (mismatches == 0) {
            std::cout << "PASS: lean iterator kmer words match CanonicalKmerIterator.\n";
        } else {
            std::cout << "FAIL: " << mismatches << " mismatches found.\n";
        }
        spdlog_piscem::drop_all();
        return mismatches > 0 ? 1 : 0;
    }

    uint64_t found = 0;
    uint64_t num_kmers = 0;
    uint64_t extensions = 0;
    uint64_t searches = 0;

    auto run_sshash_native = [&]<bool canonical>() {
        sshash::streaming_query<piscem::piscem_dictionary, canonical> sq(ri.get_dict());

        spdlog_piscem::info("starting benchmark (sshash-native, canonical={}, locate={})",
                            canonical, locate);
        auto t_start = std::chrono::high_resolution_clock::now();

        for (auto& seq : sequences) {
            if (seq.size() < k) continue;
            sq.reset();
            const char* data = seq.data();
            uint64_t n_kmers = seq.size() - k + 1;
            for (uint64_t i = 0; i < n_kmers; ++i) {
                auto res = sq.lookup(data + i);
                num_kmers++;
                if (res.kmer_id != sshash::constants::invalid_uint64) {
                    found++;
                }
            }
        }

        auto t_stop = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_stop - t_start);
        double ns_per_kmer = static_cast<double>(elapsed.count()) / num_kmers;
        extensions = sq.num_extensions();
        searches = sq.num_searches();

        std::cout << "==== streaming lookup report (sshash-native):\n";
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
    };

    if (point_lookup) {
        found = 0;
        num_kmers = 0;

        spdlog_piscem::info("starting benchmark (point-lookup)");
        auto t_start = std::chrono::high_resolution_clock::now();

        for (auto& seq : sequences) {
            if (seq.size() < k) continue;
            const char* data = seq.data();
            uint64_t n_kmers = seq.size() - k + 1;
            for (uint64_t i = 0; i < n_kmers; ++i) {
                auto res = ri.get_dict()->lookup(data + i, true);
                num_kmers++;
                if (res.kmer_id != sshash::constants::invalid_uint64) {
                    found++;
                }
            }
        }

        auto t_stop = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_stop - t_start);
        double ns_per_kmer = static_cast<double>(elapsed.count()) / num_kmers;

        std::cout << "==== streaming lookup report (point-lookup):\n";
        std::cout << "num_kmers = " << num_kmers << "\n";
        std::cout << "found_kmers = " << found << " ("
                  << (num_kmers > 0 ? static_cast<double>(found) / num_kmers * 100.0 : 0)
                  << "%)\n";
        std::cout << "time_per_kmer = " << ns_per_kmer << " ns\n";
        std::cout << "total_time = " << elapsed.count() / 1e9 << " s\n";
    } else if (lean) {
        piscem::lean_read_iterator lit(ri.get_dict(), ri.get_contig_table());

        spdlog_piscem::info("starting benchmark (lean-iterator, locate={})", locate);
        auto t_start = std::chrono::high_resolution_clock::now();

        for (auto& seq : sequences) {
            if (seq.size() < k) continue;
            lit.start(seq.data(), static_cast<int32_t>(seq.size()));

            while (!lit.is_exhausted()) {
                auto res = lit.streaming_lookup();
                num_kmers++;
                if (res.kmer_id != sshash::constants::invalid_uint64) {
                    found++;
                    if (locate) {
                        for (auto v : lit.contig_span()) {
                            auto pos = sshash::util::pos(v);
                            auto ori = sshash::util::orientation(v);
                            (void)pos;
                            (void)ori;
                        }
                    }
                }
                ++lit;
            }
        }

        auto t_stop = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_stop - t_start);
        double ns_per_kmer = static_cast<double>(elapsed.count()) / num_kmers;
        extensions = lit.num_extensions();
        searches = lit.num_searches();

        std::cout << "==== streaming lookup report (lean-iterator):\n";
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
    } else if (sshash_native) {
        if (ri.get_dict()->canonical()) {
            run_sshash_native.operator()<true>();
        } else {
            run_sshash_native.operator()<false>();
        }
    } else {
        // Use piscem-cpp's streaming query wrapper (no unitig_end_cache)
        piscem::streaming_query<false> q(ri.get_dict());

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
        extensions = q.num_extensions();
        searches = q.num_searches();

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
    }

    spdlog_piscem::drop_all();
    return 0;
}
