#ifndef __PISCEM_META_INFO__
#define __PISCEM_META_INFO__

#include "../include/json.hpp"
#include "../include/ghc/filesystem.hpp"

namespace piscem {

namespace meta_info {
class run_stats {
public:
    run_stats() = default;
    inline void cmd_line(std::string& cmd_line_in) { cmd_line_ = cmd_line_in; }
    inline void num_reads(uint64_t num_reads_in) { num_reads_ = num_reads_in; }
    inline void num_hits(uint64_t num_hits_in) { num_hits_ = num_hits_in; }
    inline void num_seconds(double num_sec) { num_seconds_ = num_sec; }

    inline std::string cmd_line() const { return cmd_line_; }
    inline uint64_t num_reads() const { return num_reads_; }
    inline uint64_t num_hits() const { return num_hits_; }
    inline double num_seconds() const { return num_seconds_; }

private:
    std::string cmd_line_{""};
    uint64_t num_reads_{0};
    uint64_t num_hits_{0};
    double num_seconds_{0};
};

inline bool write_map_info(run_stats& rs, ghc::filesystem::path& map_info_file_path) {
    using json = nlohmann::json;

    json j;
    j["cmdline"] = rs.cmd_line();
    j["num_reads"] = rs.num_reads();
    j["num_mapped"] = rs.num_hits();
    double percent_mapped = (100.0 * static_cast<double>(rs.num_hits())) / rs.num_reads();
    j["percent_mapped"] = percent_mapped;
    j["runtime_seconds"] = rs.num_seconds();
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
