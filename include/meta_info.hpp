#ifndef __PISCEM_META_INFO__
#define __PISCEM_META_INFO__

#include <optional>

#include "../include/ref_sig_info.hpp"
#include "../include/json.hpp"
#include "../include/ghc/filesystem.hpp"

namespace piscem {

namespace meta_info {
class run_stats {
  using ParamMapT = std::unordered_map<std::string, std::string>;
public:
    run_stats() = default;
    inline void cmd_line(std::string& cmd_line_in) { cmd_line_ = cmd_line_in; }
    inline void num_reads(uint64_t num_reads_in) { num_reads_ = num_reads_in; }
    inline void num_hits(uint64_t num_hits_in) { num_hits_ = num_hits_in; }
    inline void num_poisoned(uint64_t num_poisoned_in) { num_poisoned_ = num_poisoned_in; }
    inline void num_kmatch(uint64_t num_kmatch_in) { num_k_match_ = num_kmatch_in; }
    inline void num_seconds(double num_sec) { num_seconds_ = num_sec; }
    inline void important_params(ParamMapT param_map) { important_params_ = param_map; };
    inline void ref_sig_info(std::optional<ref_sig_info_t> ref_sig_info) { ref_sig_info_ = ref_sig_info; };


    inline std::string cmd_line() const { return cmd_line_; }
    inline uint64_t num_reads() const { return num_reads_; }
    inline uint64_t num_hits() const { return num_hits_; }
    inline uint64_t num_poisoned() const { return num_poisoned_; }
    inline uint64_t num_kmatch() const { return num_k_match_; }    
    inline double num_seconds() const { return num_seconds_; }
    inline const ParamMapT& important_params() const { return important_params_; }
    inline const std::optional<ref_sig_info_t>& ref_sig_info() const { return ref_sig_info_; }

private:
    std::string cmd_line_{""};
    uint64_t num_reads_{0};
    uint64_t num_hits_{0};
    uint64_t num_poisoned_{0};
    uint64_t num_k_match_{0};

    double num_seconds_{0};
    ParamMapT important_params_;
    std::optional<ref_sig_info_t> ref_sig_info_{std::nullopt};
};

inline bool write_map_info(run_stats& rs, ghc::filesystem::path& map_info_file_path) {
    using json = nlohmann::json;

    json j;
    j["cmdline"] = rs.cmd_line();
    j["num_reads"] = rs.num_reads();
    j["num_mapped"] = rs.num_hits();
    j["num_poisoned"] = rs.num_poisoned();
    double percent_mapped = (100.0 * static_cast<double>(rs.num_hits())) / rs.num_reads();
    j["percent_mapped"] = percent_mapped;
    j["ks_matched"] = rs.num_kmatch();
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
