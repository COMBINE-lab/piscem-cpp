#ifndef __REF_SIG_INFO__
#define __REF_SIG_INFO__

#include <optional>
#include "../include/json.hpp"
#include "../include/ghc/filesystem.hpp"

struct ref_sig_info_t {
  std::string sha256_names;
  std::string sha256_seqs;
  std::string sha512_names;
  std::string sha512_seqs;

  void add_to_json(nlohmann::json& json_to_augment) const {
    json_to_augment["signatures"] = { 
        {"sha256_names", sha256_names},
        {"sha512_names", sha512_names},
        {"sha256_seqs", sha256_seqs},
        {"sha512_seqs", sha512_seqs}
    };
  }

  static std::optional<ref_sig_info_t> from_path(const std::string& fpath) {
    using json = nlohmann::json;
    if (ghc::filesystem::exists(fpath)) {
      std::ifstream f(fpath);
      if (!f.good()) {
        spdlog_piscem::warn("The signature file {} exists, but couldn't be opened properly!", fpath);
      } else {
        json data = json::parse(f);
        ref_sig_info_t si;
        si.sha256_names = data["sha256_names"];
        si.sha512_names = data["sha512_names"];
        si.sha256_seqs = data["sha256_seqs"];
        si.sha512_seqs = data["sha512_seqs"];
        return std::optional<ref_sig_info_t>(si);
      }
    }

    return std::nullopt;
  }
};

#endif //__REF_SIG_INFO__
