#pragma once

#include <fstream>

#include "ghc/filesystem.hpp"
#include "itlib/small_vector.hpp"
#include "json.hpp"
#include "parallel_hashmap/phmap.h"
#include "parallel_hashmap/phmap_dump.h"
#include "../external/sshash/include/util.hpp"
#include "util_piscem.hpp"
#include "spdlog_piscem/spdlog.h"

#include "bitsery/adapter/stream.h"
#include "bitsery/bitsery.h"
#include "bitsery/brief_syntax.h"
#include "bitsery/brief_syntax/string.h"
#include "bitsery/brief_syntax/vector.h"

using poison_map_t =
  phmap::flat_hash_map<uint64_t, uint64_t, sshash::RobinHoodHash>;
using sshash::labeled_poison_occ_t;
using sshash::poison_occ_t;

// This class provides an implementation of "poison" k-mers (aka
// "distinguishing flanking k-mers" or DFKs)[1] to aid in improving
// the specificity of pseudoalignment *without alignment scoring*.
// When alignment scoring is available, a standard decoy-based approach
// can be used[2].  However, when only pseudoalignment or variants of
// pseudoalignemnt are being performed, and no reliable alignment score
// is computed, this approach can be used to improve mapping specificity.
// It relies on identifying "key" k-mers from the decoy sequence that
// would otherwise impinge on the unitigs of the reference de Bruijn
// graph.  If the resulting pseudoalignment contains such a "poison" k-mer
// then the pseudoalignment itself can be considered as "poisoned" and
// the mappings for the read discarded.
//
// The implmentation here keeps the poison k-mer information separate
// from the main index. The hits (correspondences between the reads and
// positions on unitigs of the refernece de Bruijn graph) can then be
// inspected for inclusion of relevant poison k-mers.  The implementation
// here works differently from that of [1], but the core idea and the goal
// is the same.  However, this approach does _not_ attempt to insert the
// poisoned k-mers into the primary index itself, but rather keeps a separate
// poison k-mer table (what this class implements).  This design provides
// certain tradeoffs, but works well with an otherwise "static" compacted
// reference de Bruijn graph based index.
//
// [1] Hjörleifsson, K. E., Sullivan, D. K., Holley, G., Melsted, P. & Pachter,
// L. Accurate quantification of single-nucleus and single-cell RNA-seq
// transcripts. bioRxiv.
// https://www.biorxiv.org/content/early/2022/12/02/2022. 12.02.518832 (2022)
//
// [2] Srivastava, A. et al. Alignment and mapping methodology influence tran-
// script abundance estimation. Genome Biology 21.
// https://doi.org/10.1186/s13059-020-02151-8 (Sept. 2020)
class poison_table {
public:
  // Check if the necessary files exist (prefixed with the provided
  // `basename`) from  which to load the poison table.
  static auto exists(const std::string &basename) -> bool {
    std::string pmap_name = basename + ".poison";

    return ghc::filesystem::exists(pmap_name);
  }

  // construct, copy, or assign the poison table.
  poison_table() = default;
  poison_table(const poison_table &other) = delete;
  poison_table &operator=(const poison_table &other) = delete;

  // move into this poison table from another.
  poison_table(poison_table &&other)
    : poison_map_(std::move(other.poison_map_)),
      offsets_(std::move(other.offsets_)),
      poison_occs_(std::move(other.poison_occs_)) {}

  // move assign into this poison table from another.
  poison_table &operator=(poison_table &&other) {
    poison_map_ = std::move(other.poison_map_);
    offsets_ = std::move(other.offsets_);
    poison_occs_ = std::move(other.poison_occs_);
    return *this;
  }

  // load the relevant poison table (with filed prefixed
  // with the provided `basename`) from disk.
  poison_table(const std::string &basename) {
    std::string pmap_name = basename + ".poison";
    if (!ghc::filesystem::exists(pmap_name)) {
      spdlog_piscem::critical("Can not load poison map as {} does not exist!",
                              pmap_name);
      std::exit(1);
    }
  
    {
      spdlog_piscem::info("Loading poison occ table...");
      std::ifstream poc_file{pmap_name.c_str(), std::ios::binary | std::ios::in};
      if (!poc_file.good()) {
        spdlog_piscem::critical("Error opening poison table {}.",
                                pmap_name);
        std::exit(1);
      }
      auto state = bitsery::quickDeserialization<bitsery::InputStreamAdapter>(
        poc_file, offsets_);
      state = bitsery::quickDeserialization<bitsery::InputStreamAdapter>(
        poc_file, poison_occs_);

      // now that we have read the occ_table, we can safely move 
      // the underlying input stream object into the BinaryInputArchive
      // to load the hash_map.
      phmap::BinaryInputArchive ar_in(std::move(poc_file));
      spdlog_piscem::info("Loading poison map...");
      poison_map_.phmap_load(ar_in);
      spdlog_piscem::info("done");
    }
  }

  // return a reference to the poison "map" used by this table.
  // the poison map is simply the hash table from poison k-mers
  // to their offset in the poison table.
  poison_map_t &poison_map() { return poison_map_; }

  // build the poison table data structure from a collection of
  // (possibly redundant) labeled poision k-mer occurrences.
  bool build_from_occs(std::vector<labeled_poison_occ_t> &poison_occs) {
    // we will build a map from each k-mer to the list of
    // unitigs and positions where it occurs.  Therefore, here
    // we want to sort by k-mer, then by unitig, then position
    std::sort(
      poison_occs.begin(), poison_occs.end(),
      [](const labeled_poison_occ_t &a, const labeled_poison_occ_t &b) -> bool {
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
    poison_occs.erase(std::unique(poison_occs.begin(), poison_occs.end()),
                      poison_occs.end());
    spdlog_piscem::info("Total number of distinct poison k-mer occs: {}",
                        poison_occs.size());
    using offset_t = decltype(offsets_)::value_type;
    auto max_offset = std::numeric_limits<offset_t>::max();
    if (poison_occs.size() > max_offset) {
      spdlog_piscem::critical(
        "There were {} total poison k-mer occurrences, but the offset vector "
        "can only hold values up to {}",
        poison_occs.size(), max_offset);
      return false;
    }

    // build the overall map
    offsets_.reserve(poison_occs.size() + 1);

    size_t max_range = 0;
    auto occ_it = poison_occs.begin();
    auto range_start_it = occ_it;
    size_t max_range_offset = 0;
    poison_map_[occ_it->canonical_kmer] =
      static_cast<uint64_t>(offsets_.size());
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

        offset_t dist_from_start =
          static_cast<offset_t>(std::distance(poison_occs.begin(), occ_it));
        poison_map_[occ_it->canonical_kmer] =
          static_cast<uint64_t>(offsets_.size());
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
    offset_t dist_from_start =
      static_cast<offset_t>(std::distance(poison_occs.begin(), occ_it));
    if (!offsets_.empty() and (dist_from_start != offsets_.back())) {
      offsets_.push_back(dist_from_start);
    }
    offsets_.shrink_to_fit();

    spdlog_piscem::info("[poison_table]: The most frequently occuring poison "
                        "k-mer appeared in {} distinct unitig positions.",
                        max_range);

    // copy without the canonical k-mer.
    poison_occs_.reserve(poison_occs.size());
    std::copy(poison_occs.begin(), poison_occs.end(),
              std::back_inserter(poison_occs_));

    // print out the occurrences of the most frequent poison k-mer
    auto start = static_cast<int64_t>(offsets_[max_range_offset]);
    auto end = static_cast<int64_t>(offsets_[max_range_offset + 1]);
    for (int64_t i = start; i < end; ++i) {
      std::cerr << "occ: " << i - start << ", " << poison_occs[i];
    }
    return true;
  }

  // save this poison table to disk in files prefixed by the provided
  // `output_file` prefix. The `global_nk` variable is the total number
  // of k-mers inspected when building this poison table and is passed in
  // only because it is stored in the serialized JSON file that provides
  // high-level stats about this poison table. It is not technically
  // required for the correctness or use of the table later.
  bool save_to_file(const std::string &output_file, uint64_t global_nk) {
    {
      std::string poc_filename = output_file; 
      // first serialize the offsets and occ table to the 
      // bitsery stream. Be sure to flush the adapter so the 
      // internal stream is at the right point when we go to 
      // write the hash map.
      std::ofstream poc_file{poc_filename.c_str(),
                            std::ios::binary | std::ios::trunc | std::ios::out};
      bitsery::Serializer<bitsery::OutputBufferedStreamAdapter> ser{poc_file};
      if (!poc_file.good()) {
        spdlog_piscem::critical("could not open occ output file {}",
                                poc_filename);
        return false;
      }
      ser(offsets_);
      ser.adapter().flush();
      ser(poison_occs_);
      ser.adapter().flush();

      // now that we wrote the occ table, we can safely move the ofstream 
      // object to which we are writing to the hash_map BinaryOutputArchive
      // and save the hash table itself at the end of the archive.
      phmap::BinaryOutputArchive ar_out(std::move(poc_file));
      poison_map_.phmap_dump(ar_out);
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

    spdlog_piscem::info("[poison_table]: Examined {} total decoy k-mers, "
                        "recorded {} poison k-mers.",
                        global_nk, poison_map_.size());

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
      if (!o.good()) {
        return 1;
      }
      o << std::setw(4) << j << std::endl;
      if (!o) {
        return 1;
      }
    }

    return true;
  }

  // return true if this poison table is empty (contains no
  // poisoned k-mers) and false otherwise.
  bool empty() const { return poison_map_.empty(); }

  // returns true if the key `km` exists in this table
  // (i.e. if it is a poisoned k-mer) and false otherwise.
  inline bool key_exists(uint64_t km) const {
    return poison_map_.find(km) != poison_map_.end();
  }

  // print the occurrences (unitigs and positions) for the
  // poison k-mer `km` to standard error.
  inline void print_occs(uint64_t km) const {
    auto key_it = poison_map_.find(km);
    if (key_it == poison_map_.end()) {
      return;
    }
    auto occ_start = offsets_[key_it->second];
    auto occ_end = offsets_[key_it->second + 1];
    auto it_start = poison_occs_.begin() + occ_start;
    auto it_end = poison_occs_.begin() + occ_end;
    for (; it_start != it_end; ++it_start) {
      std::cerr << "occ: " << *it_start;
    }
  }

  // returns true if the k-mer `km` is a poison k-mer that appears
  // on unitig `u1` between positions `lb` and `ub` (inclusive), and
  // false otherwise.
  inline bool key_occurs_in_unitig_between(uint64_t km, uint32_t u1,
                                           uint32_t lb, uint32_t ub) const {
    auto key_it = poison_map_.find(km);
    if (key_it == poison_map_.end()) {
      return false;
    }
    auto occ_start = offsets_[key_it->second];
    auto occ_end = offsets_[key_it->second + 1];
    auto it_start = poison_occs_.begin() + occ_start;
    auto it_end = poison_occs_.begin() + occ_end;

    for (; it_start != it_end; ++it_start) {
      bool found = (it_start->unitig_id == u1) and
                   (it_start->unitig_pos >= lb) and
                   (it_start->unitig_pos <= ub);
      if (found) {
        return true;
      }
    }
    return false;
  }

  // returns true if k-mer `km` occurs attached to unitig u1.
  inline bool key_occurs_in_unitig(uint64_t km, uint32_t u1) const {
    auto key_it = poison_map_.find(km);
    if (key_it == poison_map_.end()) {
      return false;
    }
    auto occ_start = offsets_[key_it->second];
    auto occ_end = offsets_[key_it->second + 1];
    auto it_start = poison_occs_.begin() + occ_start;
    auto it_end = poison_occs_.begin() + occ_end;
    for (; it_start != it_end; ++it_start) {
      bool found = (it_start->unitig_id == u1);
      if (found) {
        return true;
      }
    }
    return false;
  }

  // temporary until we define a proper iterator
  // returns an iterator in the underlying map corresponding to the
  // k-mer `km`.
  inline poison_map_t::iterator find_kmer(uint64_t km) {
    return poison_map_.find(km);
  }
  // returns the end iterator of the underlying map.
  inline poison_map_t::iterator kmer_end() { return poison_map_.end(); }

private:
  poison_map_t poison_map_;
  std::vector<uint32_t> offsets_;
  std::vector<poison_occ_t> poison_occs_;
};
