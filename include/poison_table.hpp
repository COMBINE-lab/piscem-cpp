#pragma once

#include <fstream>

#include "spdlog_piscem/spdlog.h"
#include "../include/json.hpp"
#include "../include/ghc/filesystem.hpp"
#include "../include/util.hpp"
#include "../include/parallel_hashmap/phmap.h"
#include "../include/parallel_hashmap/phmap_dump.h"
#include "../include/itlib/small_vector.hpp"

using poison_map_t = phmap::flat_hash_map<uint64_t, uint64_t, sshash::RobinHoodHash>;
using sshash::labeled_poison_occ_t;
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
      poc_file.read(reinterpret_cast<char*>(offsets_.data()), s * sizeof(decltype(offsets_)::value_type));

      poc_file.read(reinterpret_cast<char*>(&s), sizeof(s));
      poison_occs_.resize(s);
      poc_file.read(reinterpret_cast<char*>(poison_occs_.data()), s * sizeof(decltype(poison_occs_)::value_type));
    }
  }
  
  poison_map_t& poison_map() { return poison_map_; }

  bool build_from_occs(std::vector<labeled_poison_occ_t>& poison_occs) {
    // we will build a map from each k-mer to the list of 
    // unitigs and positions where it occurs.  Therefore, here 
    // we want to sort by k-mer, then by unitig, then position
    std::sort(poison_occs.begin(), poison_occs.end(),
              [](const labeled_poison_occ_t& a, const labeled_poison_occ_t& b) -> bool {
              if (a.canonical_kmer == b.canonical_kmer) {
              if (a.unitig_id == b.unitig_id) {
              return a.unitig_pos < b.unitig_pos;
              } else {
              return a.unitig_id < b.unitig_id;
              }
              } else {
              return a.canonical_kmer < b.canonical_kmer;
              }
              });
    // remove duplicates
    poison_occs.erase( std::unique(poison_occs.begin(), poison_occs.end()), poison_occs.end());
    spdlog_piscem::info("Total number of distinct poison k-mer occs: {}", poison_occs.size());
    using offset_t = decltype(offsets_)::value_type;
    auto max_offset = std::numeric_limits<offset_t>::max();
    if (poison_occs.size() > max_offset) {
        spdlog_piscem::critical("There were {} total poison k-mer occurrences, but the offset vector can only hold values up to {}",
                                poison_occs.size(),
                                max_offset);
        return false;
    }

    // build the overall map
    offsets_.reserve(poison_occs.size()+1);

    size_t max_range = 0;
    auto occ_it = poison_occs.begin();
    auto range_start_it = occ_it;
    size_t max_range_offset = 0;
    poison_map_[occ_it->canonical_kmer] = static_cast<uint64_t>(offsets_.size());
    offsets_.push_back(0);
    while (occ_it != poison_occs.end()) {
      // we started a new range, push back the starting point of 
      // the this range.
      if (occ_it->canonical_kmer != range_start_it->canonical_kmer) {
        size_t range_len = std::distance(range_start_it, occ_it);
        if (range_len > max_range) {
          max_range = range_len;
          max_range_offset = offsets_.size() - 1;
        }

        offset_t dist_from_start = static_cast<offset_t>(std::distance(poison_occs.begin(), occ_it));
        poison_map_[occ_it->canonical_kmer] = static_cast<uint64_t>(offsets_.size());
        offsets_.push_back(dist_from_start);
        range_start_it = occ_it;
      }
      occ_it++;
    }
    // don't forget the last one
    size_t range_len = std::distance(range_start_it, occ_it);
    if (range_len > max_range) {
      max_range = range_len;
      max_range_offset = offsets_.size() - 1;
    }
    offset_t dist_from_start = static_cast<offset_t>(std::distance(poison_occs.begin(), occ_it));
    if (!offsets_.empty() and (dist_from_start != offsets_.back())) {
      offsets_.push_back(dist_from_start); 
    }
    offsets_.shrink_to_fit();

    spdlog_piscem::info("[poison_table]: The most frequently occuring poison k-mer appeared in {} distinct unitig positions.", max_range);

    // copy without the canonical k-mer.
    poison_occs_.reserve(poison_occs.size());
    std::copy(poison_occs.begin(), poison_occs.end(), std::back_inserter(poison_occs_));

    // print out the occurrences of the most frequent poison k-mer
    auto start = static_cast<int64_t>(offsets_[max_range_offset]);
    auto end = static_cast<int64_t>(offsets_[max_range_offset+1]);
    for (int64_t i = start; i < end; ++i) {
      std::cerr << "occ: " << i-start << ", " << poison_occs[i];
    } 
    return true;
  }

  bool save_to_file(const std::string& output_file, uint64_t global_nk) {
    {
      std::string poc_filename = output_file + "_occs";
      std::ofstream poc_file(poc_filename, std::ios::binary);
      if (!poc_file.good()) {
        spdlog_piscem::critical("could not open occ output file {}", poc_filename);
        return false;
      }

      size_t s = offsets_.size();
      poc_file.write(reinterpret_cast<char*>(&s), sizeof(s));
      poc_file.write(reinterpret_cast<char*>(offsets_.data()), s * sizeof(decltype(offsets_)::value_type));

      s = poison_occs_.size();
      poc_file.write(reinterpret_cast<char*>(&s), sizeof(s));
      poc_file.write(reinterpret_cast<char*>(poison_occs_.data()), s * sizeof(decltype(poison_occs_)::value_type));
    }

    int32_t max_range = 0;
    auto oprev = offsets_.begin();
    auto ocurr = oprev + 1;
    while (ocurr != offsets_.end()) {
      int32_t d = static_cast<int32_t>(*ocurr - *oprev);
      max_range = std::max(max_range, d);
      ++ocurr;
      ++oprev;
    }

    spdlog_piscem::info("[poison_table]: Examined {} total decoy k-mers, recorded {} poison k-mers.", global_nk, poison_map_.size()); 

    phmap::BinaryOutputArchive ar_out(output_file.c_str());
    poison_map_.phmap_dump(ar_out);

    {
      std::string pinf_filename = output_file + ".json";
      using json = nlohmann::json;
      json j;
      j["num_poison_kmers"] = poison_map_.size();
      j["num_poison_occs"] = poison_occs_.size();
      j["num_observed_decoy_kmers"] = global_nk;
      j["max_poison_occ"] = max_range;
      // write prettified JSON to another file
      std::ofstream o(pinf_filename);
      if (!o.good()) { return 1; }
      o << std::setw(4) << j << std::endl;
      if (!o) { return 1; }
    }

    return true;
  }

  bool empty() const { return poison_map_.empty(); }
  
  inline bool key_exists(uint64_t km) const { return poison_map_.find(km) != poison_map_.end(); }

  inline void print_occs(uint64_t km) const { 
    auto key_it = poison_map_.find(km);
    if (key_it == poison_map_.end()) { return; }
    auto occ_start = offsets_[key_it->second];
    auto occ_end = offsets_[key_it->second+1];
    auto it_start = poison_occs_.begin() + occ_start;
    auto it_end = poison_occs_.begin() + occ_end;
    for (; it_start != it_end; ++it_start) {
      std::cerr << "occ: " << *it_start;
    }
  }

  inline bool key_occurs_in_unitig_between(uint64_t km, uint32_t u1, uint32_t lb, uint32_t ub) const { 
    auto key_it = poison_map_.find(km);
    if (key_it == poison_map_.end()) { return false; }
    auto occ_start = offsets_[key_it->second];
    auto occ_end = offsets_[key_it->second+1];
    auto it_start = poison_occs_.begin() + occ_start;
    auto it_end = poison_occs_.begin() + occ_end;

    for (; it_start != it_end; ++it_start) {
      bool found = (it_start->unitig_id == u1) and 
        (it_start->unitig_pos >= lb) and 
        (it_start->unitig_pos <= ub);
      if (found) { return true; }
    }
    return false;
  }

  inline bool key_occurs_in_unitig(uint64_t km, uint32_t u1) const { 
    auto key_it = poison_map_.find(km);
    if (key_it == poison_map_.end()) { return false; }
    auto occ_start = offsets_[key_it->second];
    auto occ_end = offsets_[key_it->second+1];
    auto it_start = poison_occs_.begin() + occ_start;
    auto it_end = poison_occs_.begin() + occ_end;
    for (; it_start != it_end; ++it_start) {
      bool found = (it_start->unitig_id == u1);
      if (found) { return true; }
    }
    return false;
  }

 
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

  inline bool key_occurs_in_unitigs(uint64_t km, 
                                    itlib::small_vector<uint32_t>& unitigs) {
    auto key_it = poison_map_.find(km);
    if (key_it == poison_map_.end()) { return false; }
    auto occ_start = offsets_[key_it->second];
    auto occ_end = offsets_[key_it->second+1];
    auto it_start = poison_occs_.begin() + occ_start;
    auto it_end = poison_occs_.begin() + occ_end;

    auto u_start = unitigs.begin();
    auto u_end = unitigs.end();
    while (it_start != it_end && u_start != u_end ) {
        if (it_start->unitig_id < *u_start) {
            ++it_start;
        } else {
           if (!(*u_start < it_start->unitig_id)) {
               return true;
           }
           ++u_start;
        }
    }
    return false;
  }

  // temporary until we define a proper iterator
  inline poison_map_t::iterator find_kmer(uint64_t km) { return poison_map_.find(km); }
  inline poison_map_t::iterator kmer_end() { return poison_map_.end(); }

private:
  poison_map_t poison_map_;
  std::vector<uint32_t> offsets_;
  std::vector<poison_occ_t> poison_occs_;
};
