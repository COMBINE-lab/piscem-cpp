#ifndef HIT_SEARCHER_HPP
#define HIT_SEARCHER_HPP

#include "CanonicalKmer.hpp"
#include "CanonicalKmerIterator.hpp"
//#include "Util.hpp"

#include <algorithm>
#include <iostream>
#include <iterator>
class hit_searcher {
enum class ExpansionTerminationType : uint8_t { MISMATCH = 0, CONTIG_END, READ_END };  

public:
  explicit hit_searcher(PufferfishIndexT* pfi) : pfi_(pfi) { 
    k = pfi_->k(); 
    setChainSubOptThresh(pre_merge_chain_sub_thresh_);
  }
  
  bool get_raw_hits_sketch(std::string &read,
                  pufferfish::util::QueryCache& qc,
                  bool isLeft=false,
                  bool verbose=false);

void clear();

void setAltSkip(uint32_t altSkip);

pufferfish::util::HitFilterPolicy getHitFilterPolicy() const;

inline std::vector<std::pair<int, pufferfish::util::ProjectedHits>>& get_left_hits() { 
  return left_rawHits;
}
inline std::vector<std::pair<int, pufferfish::util::ProjectedHits>>& get_right_hits() {
  return right_rawHits;
}

private:
  reference_index* pfi_;
  size_t k;
  uint32_t altSkip{3};

  bool isSingleEnd = false;
  std::vector<std::pair<int, pufferfish::util::ProjectedHits>> left_rawHits;
  std::vector<std::pair<int, pufferfish::util::ProjectedHits>> right_rawHits;
};

#endif // HIT_SEARCHER