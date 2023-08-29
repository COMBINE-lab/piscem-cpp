#ifndef HIT_SEARCHER_HPP
#define HIT_SEARCHER_HPP

#include "CanonicalKmer.hpp"
#include "CanonicalKmerIterator.hpp"
#include "projected_hits.hpp"
#include "reference_index.hpp"
//#include "query/streaming_query_canonical_parsing.cpp"
#include "query/streaming_query_canonical_parsing.hpp"
//#include "Util.hpp"
//#include "dictionary.hpp"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <optional>

namespace mindex {

// "smart enum idea from"
// https://stackoverflow.com/questions/21295935/can-a-c-enum-class-have-methods
class SkippingStrategy {
//enum class SkippingStrategy : uint8_t { STRICT = 0, PERMISSIVE };
public:
  enum Value : uint8_t {
    STRICT = 0,
    PERMISSIVE
  };

  SkippingStrategy() = default;
  static std::optional<SkippingStrategy> from_string(const std::string& s) {
    if (s == "strict") {
      return STRICT;
    } else if (s == "permissive") {
      return PERMISSIVE;
    } else {
      return std::nullopt;
    }
  }
  constexpr SkippingStrategy(Value val) : value_(val) {}
  constexpr bool operator==(SkippingStrategy a) const { return value_ == a.value_; }
  constexpr bool operator!=(SkippingStrategy a) const { return value_ != a.value_; }
   
private:
  Value value_;
};

class hit_searcher {
enum class ExpansionTerminationType : uint8_t { MISMATCH = 0, CONTIG_END, READ_END };  

public:
  explicit hit_searcher(reference_index* pfi) : pfi_(pfi) { 
    k = static_cast<size_t>(pfi_->k()); 
  }
  
  bool get_raw_hits_sketch(std::string &read,
                  sshash::streaming_query_canonical_parsing& qc,
                  SkippingStrategy strat = SkippingStrategy::STRICT,
                  bool isLeft=false,
                  bool verbose=false);
  
  bool get_raw_hits_sketch_orig(std::string &read,
                  sshash::streaming_query_canonical_parsing& qc,
                  SkippingStrategy strat = SkippingStrategy::STRICT,
                  bool isLeft=false,
                  bool verbose=false);

void clear();

void setAltSkip(uint32_t altSkip);

inline std::vector<std::pair<int, projected_hits>>& get_left_hits() { 
  return left_rawHits;
}
inline std::vector<std::pair<int, projected_hits>>& get_right_hits() {
  return right_rawHits;
}

inline reference_index* get_index() const { return pfi_; }

private:
  reference_index* pfi_;
  size_t k;
  uint32_t altSkip{3};

  bool isSingleEnd = false;
  std::vector<std::pair<int, projected_hits>> left_rawHits;
  std::vector<std::pair<int, projected_hits>> right_rawHits;
};
}
#endif // HIT_SEARCHER
