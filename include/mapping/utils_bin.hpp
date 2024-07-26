// #pragma once

// #include "../include/parallel_hashmap/phmap.h"
// #include "../include/reference_index.hpp"
// #include "../include/CanonicalKmerIterator.hpp"
// #include "../include/query/streaming_query_canonical_parsing.hpp"
// #include "../include/projected_hits.hpp"
// #include "../include/FastxParser.hpp"
// #include "../include/itlib/small_vector.hpp"

// #include "../include/hit_searcher.hpp"

// #include <algorithm>
// #include <limits>
// #include <vector>
// #include <cassert>
// #include <fstream>
// #include <cmath>  // for std::ceil on linux
// #include <numeric>
// #include <type_traits>



// namespace mapping {

// namespace util_bin {

// constexpr int32_t invalid_frag_len = std::numeric_limits<int32_t>::min();
// constexpr int32_t invalid_mate_pos = std::numeric_limits<int32_t>::min();

// enum class orientation_filter : uint8_t { NONE, FORWARD_ONLY, RC_ONLY };

// class bin_pos {
//     public:
//         explicit bin_pos(mindex::reference_index* pfi,
//                         float thr = 0.7,
//                         uint64_t bin_size = 2000,
//                         uint64_t overlap = 300
//                         )  { 
//             pfi_ = pfi,
//             thr = thr,
//             bin_size = bin_size,
//             overlap = overlap,
//             compute_cum_rank();
//         };

//         uint64_t get_cum_len(size_t i) {
//             return cum_ref_lens[i];
//         }

//         float get_thr() {
//             return thr;
//         }
        
//         std::pair<uint64_t, uint64_t> get_bin_id(uint64_t tid, uint64_t pos) {
//             uint64_t cum_len = get_cum_len(tid);
//             assert(bin_size > overlap);
//             uint64_t bin1 = (cum_len + pos + 1)/bin_size; // 1 added since 0 based
//             uint64_t bin2 = (cum_len + pos + 1) > (bin1+1)*bin_size-overlap ? (bin1+1) :
//                 std::numeric_limits<uint64_t>::max(); // std::numeric_limits<uint32_t>::max() indicates that the kmer does not belong to the overlapping region
//             return {bin1, bin2};
//         }
//         mindex::reference_index* get_ref() { return pfi_;}


//     private:
//         mindex::reference_index* pfi_;
//         std::vector<uint64_t> cum_ref_lens;
//         float thr;
//         uint64_t bin_size;
//         uint64_t overlap;
//         void compute_cum_rank() {
//             int32_t n_refs = static_cast<int32_t>(pfi_->num_refs());
//             cum_ref_lens.reserve(n_refs);
//             cum_ref_lens[0] = 0;
//             for(int32_t i = 1; i < n_refs; i++) {
//                 cum_ref_lens[i] = cum_ref_lens[i-1] + pfi_->ref_len(i-1);
//             }
//         }
// };

// // struct simple_hit {
// //     bool is_fw{false};
// //     bool mate_is_fw{false};
// //     int32_t pos{-1};
// //     float score{0.0};
// //     uint32_t num_hits{0};
// //     uint32_t tid{std::numeric_limits<uint32_t>::max()};
// //     uint64_t bin_id{std::numeric_limits<uint64_t>::max()};
// //     int32_t mate_pos{std::numeric_limits<int32_t>::max()};
// //     int32_t fragment_length{std::numeric_limits<int32_t>::max()};
// //     inline bool valid_pos(int32_t read_len, uint32_t txp_len, int32_t max_over) {
// //         int32_t signed_txp_len = static_cast<int32_t>(txp_len);
// //         return (pos > -max_over) and ((pos + read_len) < (signed_txp_len + max_over));
// //     }
// //     inline bool has_mate() const { return mate_pos != invalid_mate_pos; }
// //     inline bool mate_is_mapped() const { return mate_pos != invalid_mate_pos; }
// //     inline int32_t frag_len() const {
// //         return (fragment_length != invalid_frag_len) ? fragment_length : 0;
// //     }
// // };

// // enum class MappingType : uint8_t {
// //     UNMAPPED = 0,
// //     SINGLE_MAPPED = 1,
// //     MAPPED_FIRST_ORPHAN = 2,
// //     MAPPED_SECOND_ORPHAN = 3,
// //     MAPPED_PAIR = 4
// // };

// void print_hits(const std::vector<mapping::util::simple_hit> &hits ) {
//     for (const auto& hit : hits) {
//         std::cout << "isFw: " << hit.is_fw << std::endl;
//         std::cout << "pos: " << hit.pos << std::endl;
//         std::cout << "num hits: " << hit.num_hits << std::endl;
//         std::cout << "tid: " << hit.tid << std::endl;
//         std::cout << "bin_id: " << hit.bin_id << std::endl;
//         std::cout << "------------------------" << std::endl;
//     }
// }
 
// template <typename mapping_cache_info_t>
// inline void merge_se_mappings(mapping_cache_info_t& map_cache_left,
//                               mapping_cache_info_t& map_cache_right, int32_t left_len,
//                               int32_t right_len, mapping_cache_info_t& map_cache_out 
//                               ) {
//     map_cache_out.clear();
//     auto& accepted_left = map_cache_left.accepted_hits;
//     auto& accepted_right = map_cache_right.accepted_hits;

//     size_t had_matching_kmers_left = map_cache_left.has_matching_kmers;
//     size_t had_matching_kmers_right = map_cache_right.has_matching_kmers;

//     size_t num_accepted_left = accepted_left.size();
//     size_t num_accepted_right = accepted_right.size();
//     std::unordered_map<int32_t, int8_t> hit_pos; // A kmer can map to two bins for the same position, we only want 1 entry
//     // std::cout << "num hits " << num_accepted_left << " " << num_accepted_right << "\n";
//     // std::cout << "matching kmers " << had_matching_kmers_left << " " << had_matching_kmers_right << "\n";
//     if ((num_accepted_left > 0) and (num_accepted_right > 0)) {
//         // std::cout << "entered both\n";
//         // print_hits(accepted_left);
//         // std::cout << "left right\n";
//         // print_hits(accepted_right);
//         // look for paired end mappings
//         // so we have to sort our accepted hits
//         struct {
//             // sort first by orientation, then by transcript id, and finally by position
//             bool operator()(const mapping::util::simple_hit& a,
//                             const mapping::util::simple_hit& b) {
//                 if (a.is_fw != b.is_fw) { return a.is_fw > b.is_fw; }
//                 // orientations are the same
//                 if (a.bin_id != b.bin_id) { return a.bin_id < b.bin_id; }
//                 return a.pos < b.pos;
//             }
//         } simple_hit_less;
//         std::sort(accepted_left.begin(), accepted_left.end(), simple_hit_less);
//         std::sort(accepted_right.begin(), accepted_right.end(), simple_hit_less);

//         const mapping::util::simple_hit smallest_rc_hit = {false, false, -1, 0.0, static_cast<uint32_t>(0), static_cast<uint32_t>(0), 
//             static_cast<uint64_t>(0)};
//         // start of forward sub-list
//         auto first_fw1 = accepted_left.begin();
//         // end of forward sub-list is first non-forward hit
//         auto last_fw1 = std::lower_bound(accepted_left.begin(), accepted_left.end(),
//                                          smallest_rc_hit, simple_hit_less);
//         // start of rc list
//         auto first_rc1 = last_fw1;
//         // end of rc list
//         auto last_rc1 = accepted_left.end();

//         // start of forward sub-list
//         auto first_fw2 = accepted_right.begin();
//         // end of forward sub-list is first non-forward hit
//         auto last_fw2 = std::lower_bound(accepted_right.begin(), accepted_right.end(),
//                                          smallest_rc_hit, simple_hit_less);
//         // start of rc list
//         auto first_rc2 = last_fw2;
//         // end of rc list
//         auto last_rc2 = accepted_right.end();

//         auto back_inserter = std::back_inserter(map_cache_out.accepted_hits);
//         using iter_t = decltype(first_fw1);
//         using out_iter_t = decltype(back_inserter);
//         auto merge_lists = [left_len, right_len, &hit_pos](iter_t first1, iter_t last1, iter_t first2,
//                                                  iter_t last2, out_iter_t out) -> out_iter_t {
//         // auto merge_lists = [left_len, right_len, ri](iter_t first1, iter_t last1, iter_t first2,
//         //                                          iter_t last2, out_iter_t out) -> out_iter_t {
//             // https://en.cppreference.com/w/cpp/algorithm/set_intersection
//             while (first1 != last1 && first2 != last2) {
//                 if (first1->bin_id < first2->bin_id) {
//                     ++first1;
//                 } else {
//                     if (!(first2->bin_id < first1->bin_id)) {
//                         // first1->tid == first2->tid have the same transcript.
//                         int32_t pos_fw = first1->is_fw ? first1->pos : first2->pos;
//                         int32_t pos_rc = first1->is_fw ? first2->pos : first1->pos;
//                         int32_t frag_len = (pos_rc - pos_fw);
//                         // std::cout << frag_len << " fragment length\n";
//                         if (frag_len == 0) {
//                             // std::cout << "0 fragment length";
//                         }
//                         // std::cout << ri.ref_name(first1->bin_id) << " tids " << ri.ref_name(first2->bin_id) << "pos " << pos_fw << " " << pos_rc << " frag len " << frag_len << "\n";
//                         if ((-20 < frag_len) and (frag_len < 1000)) {
//                             // if left is fw and right is rc then
//                             // fragment length is (right_pos + right_len - left_pos) + 1
//                             // otherwise it is (left_pos + left_len - right_pos) + 1
//                             bool right_is_rc = !first2->is_fw;
//                             int32_t tlen = right_is_rc
//                                                ? ((first2->pos + right_len - first1->pos) + 1)
//                                                : ((first1->pos + left_len - first2->pos) + 1);
//                             if (hit_pos.find(first1->pos)==hit_pos.end()) {
//                                hit_pos[first1->pos] = 1;

//                                 *out++ = {first1->is_fw, first2->is_fw, first1->pos, 0.0, std::min(first1->num_hits, first2->num_hits),
//                                       first1->tid, first1->bin_id,   first2->pos,   tlen};
//                             }
//                             ++first1;
//                         }
//                     }
//                     ++first2;
//                 }
//             }
//             return out;
//         };

//         // find hits of form 1:fw, 2:rc
//         merge_lists(first_fw1, last_fw1, first_rc2, last_rc2, back_inserter);
//         // find hits of form 1:rc, 2:fw
//         merge_lists(first_rc1, last_rc1, first_fw2, last_fw2, back_inserter);

//         map_cache_out.map_type = (map_cache_out.accepted_hits.size() > 0) ? MappingType::MAPPED_PAIR
//                                                                           : MappingType::UNMAPPED;
//         // std::cout << "map cache sizes" << map_cache_out.accepted_hits.size() << std::endl;
//         // std::cout << "map cache2 " << map_cache_out.accepted_hits[0] << std::endl;
//     } else if ((num_accepted_left > 0) and !had_matching_kmers_right) {
//         // just return the left mappings
//         std::swap(map_cache_left.accepted_hits, map_cache_out.accepted_hits);
//         map_cache_out.map_type = (map_cache_out.accepted_hits.size() > 0)
//                                      ? MappingType::MAPPED_FIRST_ORPHAN
//                                      : MappingType::UNMAPPED;
//     } else if ((num_accepted_right > 0) and !had_matching_kmers_left) {
//         // just return the right mappings
//         std::swap(map_cache_right.accepted_hits, map_cache_out.accepted_hits);
//         map_cache_out.map_type = (map_cache_out.accepted_hits.size() > 0)
//                                      ? MappingType::MAPPED_SECOND_ORPHAN
//                                      : MappingType::UNMAPPED;
//     } else {
//         // return nothing
//     }

//     if (map_cache_out.accepted_hits.size() > 0) {
//         std::vector<mapping::util::simple_hit> accepted_hits;
//         uint32_t max_num_hits = 0;
//         for (const auto& hit:map_cache_out.accepted_hits)  {
//             max_num_hits = std::max(hit.num_hits, max_num_hits);
//         }
//         for (const auto& hit:map_cache_out.accepted_hits)  {
//             if (hit.num_hits >= max_num_hits) {
//                 accepted_hits.emplace_back(hit);
//             }
//         }
//         map_cache_out.accepted_hits = accepted_hits;    
//     }
    
//     // std::cout << "hits right\n";
//     // print_hits(map_cache_right.accepted_hits);
//     // std::cout << "hits left\n";
//     // print_hits(map_cache_left.accepted_hits);
// }

// }  // namespace util

// }  // namespace mapping
