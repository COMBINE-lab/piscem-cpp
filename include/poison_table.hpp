#pragma once

#include <fstream>

#include "spdlog_piscem/spdlog.h"
#include "../include/ghc/filesystem.hpp"
#include "../include/mapping/utils.hpp"
#include "../include/util.hpp"
#include "../include/parallel_hashmap/phmap.h"
#include "../include/parallel_hashmap/phmap_dump.h"

using poison_map_t = phmap::flat_hash_map<uint64_t, uint64_t, sshash::RobinHoodHash>;
using sshash::poison_occ_t;

class poison_table {
  public:
  static auto exists(const std::string& basename) -> bool {
    std::string pmap_name = basename + ".poison";
    std::string pocc_name = basename + ".poison_occs";

    return (ghc::filesystem::exists(pmap_name) and
            ghc::filesystem::exists(pocc_name));
  }

  poison_table() = default;
  poison_table(const poison_table& other) = delete;
  poison_table& operator=(const poison_table& other) = delete;

  poison_table(poison_table&& other) : 
    poison_map_(std::move(other.poison_map_)),
    offsets_(std::move(other.offsets_)),
    poison_occs_(std::move(other.poison_occs_)) {}
  
  poison_table& operator=(poison_table&& other) { 
    poison_map_ = std::move(other.poison_map_);
    offsets_ = std::move(other.offsets_);
    poison_occs_ = std::move(other.poison_occs_);
    return *this; 
  }

  poison_table(const std::string& basename) {
    std::string pmap_name = basename + ".poison";
    std::string pocc_name = basename + ".poison_occs";
    if (!ghc::filesystem::exists(pmap_name)) {
      spdlog_piscem::critical("Can not load poison map as {} does not exist!", pmap_name);
      std::exit(1);
    }
    if (!ghc::filesystem::exists(pmap_name)) {
      spdlog_piscem::critical("Can not load poison occ table as {} does not exist!", pocc_name);
      std::exit(1);
    }

    {
      phmap::BinaryInputArchive ar_in(pmap_name.c_str());
      spdlog_piscem::info("Loading poison map...");
      poison_map_.phmap_load(ar_in);
      spdlog_piscem::info("done");
    }

    {
      spdlog_piscem::info("Loading poison occ table...");
      std::ifstream poc_file(pocc_name, std::ios::binary);
      if (!poc_file.good()) {
        spdlog_piscem::critical("Error opening poison occ table {}.", pocc_name);
        std::exit(1);
      }
      size_t s{0};
      poc_file.read(reinterpret_cast<char*>(&s), sizeof(s));
      offsets_.resize(s);
      poc_file.read(reinterpret_cast<char*>(offsets_.data()), s * sizeof(uint64_t));

      poc_file.read(reinterpret_cast<char*>(&s), sizeof(s));
      poison_occs_.resize(s);
      poc_file.read(reinterpret_cast<char*>(poison_occs_.data()), s * sizeof(poison_occ_t));
    }
  }

  bool empty() const { return poison_map_.empty(); }
  
  inline bool key_exists(uint64_t km) const { return poison_map_.find(km) != poison_map_.end(); }

  inline bool key_occurs_in_unitigs(uint64_t km, uint32_t u1, uint32_t u2) const { 
    auto key_it = poison_map_.find(km);
    if (key_it == poison_map_.end()) { return false; }
    auto occ_start = offsets_[key_it->second];
    auto occ_end = offsets_[key_it->second+1];
    auto it_start = poison_occs_.begin() + occ_start;
    auto it_end = poison_occs_.begin() + occ_end;
    for (; it_start != it_end; ++it_start) {
      bool found = (it_start->unitig_id == u1) or (it_start->unitig_id == u2);
      if (found) { return true; }
    }
    return false;
  }

  // temporary until we define a proper iterator
  inline poison_map_t::iterator find_kmer(uint64_t km) { return poison_map_.find(km); }
  inline poison_map_t::iterator kmer_end() { return poison_map_.end(); }

private:
  poison_map_t poison_map_;
  std::vector<uint64_t> offsets_;
  std::vector<poison_occ_t> poison_occs_;
};
