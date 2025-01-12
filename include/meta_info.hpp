#ifndef __PISCEM_META_INFO__
#define __PISCEM_META_INFO__

#include <optional>

#include "ref_sig_info.hpp"
#include "json.hpp"
#include "../include/ghc/filesystem.hpp"

namespace piscem {

enum RunMode {
    bulk,
    scrna,
    scatac,
};

namespace meta_info {
class run_stats {
  using ParamMapT = std::unordered_map<std::string, std::string>;
public:
    run_stats() = default;
    inline void cmd_line(std::string& cmd_line_in) { cmd_line_ = cmd_line_in; }
    inline void mode(RunMode mode_in) { mode_ = mode_in; }
    inline void num_reads(uint64_t num_reads_in) { num_reads_ = num_reads_in; }
    inline void num_hits(uint64_t num_hits_in) { num_hits_ = num_hits_in; }
    inline void num_multihits(uint64_t num_multi_hits_in) { num_multi_ = num_multi_hits_in; }
    inline void num_poisoned(uint64_t num_poisoned_in) { num_poisoned_ = num_poisoned_in; }
    inline void num_kmatch(uint64_t num_kmatch_in) { num_k_match_ = num_kmatch_in; }
    inline void num_lmatch(uint64_t num_lmatch_in) { num_l_match_ = num_lmatch_in; }
    inline void num_rmatch(uint64_t num_rmatch_in) { num_r_match_ = num_rmatch_in; }
    inline void num_dovematch(uint64_t num_dovematch_in) { num_dove_match_ = num_dovematch_in; }
    inline void num_dovenum(uint64_t num_dovenum_in) { num_dove_num_ = num_dovenum_in; }
    inline void num_ovmatch(uint64_t num_ovmatch_in) { num_ov_match_ = num_ovmatch_in; }
    inline void num_ovnum(uint64_t num_ovnum_in) { num_ov_num_ = num_ovnum_in; }
    inline void num_rorphan(uint64_t num_rorphan_in) { num_r_orphan_ = num_rorphan_in; }
    inline void num_lorphan(uint64_t num_lorphan_in) { num_l_orphan_ = num_lorphan_in; }
    inline void num_negkmers(uint64_t num_negkmers_in) { num_neg_kmers_ = num_negkmers_in; }
    inline void num_seconds(double num_sec) { num_seconds_ = num_sec; }
    inline void important_params(ParamMapT param_map) { important_params_ = param_map; };
    inline void ref_sig_info(std::optional<ref_sig_info_t> ref_sig_info) { ref_sig_info_ = ref_sig_info; };


    inline std::string cmd_line() const { return cmd_line_; }
    inline RunMode mode () const { return mode_; }
    inline std::string mode_str () const { 
            switch(mode_){
                case RunMode::scatac:
                    return std::string{"sc-atac"};
                case RunMode::scrna:
                    return std::string{"sc-rna"};
                case RunMode::bulk:
                    return std::string{"bulk"};
            }
        // should not reach here, but the C++ compiler
        // isn't smart enough to know the switch above
        // is exhaustive.
        return std::string{"unknown"};
    }
    inline uint64_t num_reads() const { return num_reads_; }
    inline uint64_t num_hits() const { return num_hits_; }
    inline uint64_t num_multihits() const { return num_multi_; }
    inline uint64_t num_poisoned() const { return num_poisoned_; }
    inline uint64_t num_kmatch() const { return num_k_match_; }
    inline uint64_t num_lmatch() const { return num_l_match_; }
    inline uint64_t num_rmatch() const { return num_r_match_; }
    inline uint64_t num_dovematch() const { return num_dove_match_; }
    inline uint64_t num_dovenum() const { return num_dove_num_; }
    inline uint64_t num_ovmatch() const { return num_ov_match_; }
    inline uint64_t num_ovnum() const { return num_ov_num_; }
    inline uint64_t num_rorphan() const { return num_r_orphan_; }
    inline uint64_t num_lorphan() const { return num_l_orphan_; }
    inline uint64_t num_negkmers() const { return num_neg_kmers_; }
    inline double num_seconds() const { return num_seconds_; }
    inline const ParamMapT& important_params() const { return important_params_; }
    inline const std::optional<ref_sig_info_t>& ref_sig_info() const { return ref_sig_info_; }

private:
    std::string cmd_line_{""};
    RunMode mode_{RunMode::bulk};
    uint64_t num_reads_{0};
    uint64_t num_hits_{0};
    uint64_t num_multi_{0};
    uint64_t num_poisoned_{0};
    uint64_t num_k_match_{0};
    uint64_t num_l_match_{0};
    uint64_t num_r_match_{0};
    uint64_t num_dove_match_{0};
    uint64_t num_dove_num_{0};
    uint64_t num_ov_match_{0};
    uint64_t num_ov_num_{0};
    uint64_t num_r_orphan_{0};
    uint64_t num_l_orphan_{0};
    uint64_t num_neg_kmers_{0};

    double num_seconds_{0};
    ParamMapT important_params_;
    std::optional<ref_sig_info_t> ref_sig_info_{std::nullopt};
};

inline bool write_map_info(run_stats& rs, ghc::filesystem::path& map_info_file_path) {
    using json = nlohmann::json;

    json j;
    j["cmdline"] = rs.cmd_line();
    j["mode"] = rs.mode_str();
    j["num_reads"] = rs.num_reads();
    j["num_mapped"] = rs.num_hits();
    j["num_poisoned"] = rs.num_poisoned();
    double percent_mapped = (100.0 * static_cast<double>(rs.num_hits())) / rs.num_reads();
    j["percent_mapped"] = percent_mapped;

    // only ouput stats tracked in sc-atac processing
    // if we are in sc-atac mode.
    if (rs.mode() == RunMode::scatac) {
        j["num_multihits"] = rs.num_multihits();
        j["ks_matched"] = rs.num_kmatch();
        j["l_matched"] = rs.num_lmatch();
        j["r_matched"] = rs.num_rmatch();
        j["dove_matched"] = rs.num_dovematch();
        j["dove_num"] = rs.num_dovenum();
        j["ov_matched"] = rs.num_ovmatch();
        j["ov_num"] = rs.num_ovnum();
        j["l_orphan"] = rs.num_lorphan();
        j["r_orphan"] = rs.num_rorphan();
        j["num_negkmers"] = rs.num_negkmers();
    }

    j["runtime_seconds"] = rs.num_seconds();

    if (rs.ref_sig_info()) {
      rs.ref_sig_info()->add_to_json(j); 
    } 

    for (auto& [k, v] : rs.important_params()) {
      j[k] = v;
    }
    // write prettified JSON to another file
    std::ofstream o(map_info_file_path.string());

    if (!o.good()) { return false; }

    o << std::setw(4) << j << std::endl;

    if (!o) { return false; }
    return true;
}

}  // namespace meta_info
}  // namespace piscem
#endif  // __PISCEM_META_INFO__
