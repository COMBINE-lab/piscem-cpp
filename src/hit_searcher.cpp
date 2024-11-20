#include "../include/hit_searcher.hpp"
#include "../external/sshash/include/bit_vector_iterator.hpp"
#include "../include/streaming_query.hpp"
#include "../external/sshash/include/util.hpp"
#include <cmath>
#include <limits>
#include <optional>

// using spp:sparse_hash_map;

namespace mindex {
// polute the namespace --- put this in the functions that need it.
namespace kmers = combinelib::kmers;

void hit_searcher::setAltSkip(uint32_t as) { altSkip = as; }

enum class LastSkipType : uint8_t { NO_HIT = 0, SKIP_READ = 1, SKIP_UNI = 2 };

struct FastHitInfo {
    // return true if it is currently valid to do
    // a fast check, and false otherwise.
    inline bool valid() { return fast_check; }

    // set the state of this FastHitInfo instance
    // true if it is valid to check upon next query
    // and false otherwise.
    inline void valid(bool v) { fast_check = v; }

    bool fast_check{false};
    int32_t offset{0};
    uint64_t ref_kmer{0};
    KmerMatchType expected_match_type{KmerMatchType::NO_MATCH};
    bool is_confirmatory{false};
};

// Idea is to move the logic of the search into here.
// We should also consider the optimizations that can
// be done here (like having small checks expected to)
// be on the same contig bypass a hash lookup.
struct SkipContext {
    SkipContext(std::string& read, reference_index* pfi_in, int32_t k_in, uint32_t alt_skip_in)
        : kit1(read)
        , kit_tmp(read)
        , pfi(pfi_in)
        , ref_contig_it(sshash::bit_vector_iterator(pfi_in->contigs(), 0))
        , read_len(static_cast<int32_t>(read.length()))
        , read_target_pos(0)
        , read_current_pos(0)
        , read_prev_pos(0)
        , safe_skip(1)
        , k(k_in)
        , expected_cid(invalid_cid)
        , last_skip_type(LastSkipType::NO_HIT)
        , miss_it(0)
        , global_contig_pos(-1)
        , alt_skip(alt_skip_in)
        , hit_found(false) {}


    inline int32_t get_k() const { return k; }

    inline bool is_exhausted() { return kit1 == kit_end; }

    inline CanonicalKmer& curr_kmer() { return kit1->first; }

    inline bool is_good_query() { return !(kit1->first.is_homopolymer()); }

    // tries to increment the iterator by a single nucleotide, but returns 
    // the number of nucleotides we actually advanced. It could be greater 
    // than 1 because of `N`s.
    inline int increment_read_iter() {
      auto pos_before_advancement = kit1->second;
      ++kit1;
      return (kit1->second - pos_before_advancement);
    }

    inline pufferfish::CanonicalKmerIterator get_iter() const { return kit1; }
    inline void set_iter(pufferfish::CanonicalKmerIterator& kit_in) { kit1 = kit_in; }

    // tries to increment the iterator by `amount`, but returns 
    // the number of nucleotides we actually advanced. It could be greater 
    // than `amount` because of `N`s.
    inline int advance_read_iter(int32_t amount) {
      auto pos_before_advancement = kit1->second;
      kit1 += amount;
      return (kit1->second - pos_before_advancement);
    }
 
    // 
    inline KmerMatchType check_match() {
      return kit1->first.isEquivalent(fast_hit.ref_kmer);
    }

    inline bool query_kmer(piscem::streaming_query& qc) {
        phits = pfi->query(kit1, qc);
        return !phits.empty();
    }

    inline bool query_kmer_complex(piscem::streaming_query& qc,
                           bool& obtained_confirmatory_hit) {
        bool found_match = false;
        if (fast_hit.valid()) {
            auto keq = kit1->first.isEquivalent(fast_hit.ref_kmer);
            // if the k-mer we found at this position on the query matched
            // the reference k-mer on the contig at the position we expect
            // (in any orientation), then we know where this k-mer occurs in
            // the index. Note, if it doesn't match the expected orientation
            // then it might not be a confirmatory hit, but we can still avoid
            // having to query the index for it explicitly.
            if (keq != KmerMatchType::NO_MATCH) {
                found_match = true;
                // if the hit is in the expected orientation, then it's a
                // confirmatory hit.
                // Note: consider if the below should actually be "&=" instead of "=" or not.
                fast_hit.is_confirmatory = (keq == fast_hit.expected_match_type);
                // how the k-mer hits the contig (true if k-mer in fwd orientation,
                // false otherwise)
                bool hit_fw = (keq == KmerMatchType::IDENTITY_MATCH);
                phits.contigOrientation_ = hit_fw;
                phits.globalPos_ += fast_hit.offset;
                phits.contigPos_ += fast_hit.offset;
            }
            // the next query can't be resolved with a fast hit
            fast_hit.valid(false);
        }

        // if we didn't find the query k-mer at the expected reference
        // position (or we were not set up to do a fast hit) then
        // query the index.
        if (!found_match) {
            obtained_confirmatory_hit = false;
            phits = pfi->query(kit1, qc);
            hit_found = (hit_found or !phits.empty());
        } else {
            // we obtained a confirmatory result iff we satisfied
            // the query via a fast hit, and the match confirmed
            // our expectetation.
            obtained_confirmatory_hit = fast_hit.is_confirmatory;
        }

        return !phits.empty();
    }

    // Returns true if the current hit occurred
    // on a unitig other than that which was expected.
    // If there is no expecation about the unitig
    // (i.e. if `clear_expectation()` was called before
    // obtaining this hit) or if the current unitig
    // matches expectation, then it returns false.
    inline bool hit_is_unexpected() {
        return (expected_cid != invalid_cid and phits.contigIdx_ != expected_cid);
    }

    // Returns:
    // * True if there is a target position (end of unitig / read)
    //   jump that is _beyond_ the current hit's position on the read.
    // * False if the current hit is at the position to which
    //   the jump was intended.
    // This value can be used to determine the appropriate course
    // of action during safe skipping.
    inline bool hit_is_before_target_pos() {
        return (miss_it > 0 and read_target_pos > kit1->second);
    }

    // Returns the current position of the skip context
    // on the read.
    inline int32_t read_pos() { return kit1->second; }

    // Returns the most recent projected hits object
    // obtained by this skip context.
    inline projected_hits proj_hits() { return phits; }

    // Clears out any expectation we have about the unitig
    // on which the next hit should occur.
    inline void clear_expectation() { expected_cid = invalid_cid; }

    // Computes the amount by which we will move along
    // the read and query the reference until we hit the
    // original target position.  This skip is used in the
    // case that we query and find that we have arrived
    // at an unexpected unitig (i.e. JE(HD) ).
    inline int32_t compute_safe_skip() {
        // read_target_pos = //hit_is_before_target_pos() ? read_target_pos : read_target_pos-1;

        // There are 2 cases we may encounter here.
        // Either:
        // (1) The current hit was found at the original
        // jump location.  In that case we got
        // an unexpected hit. To avoid iterating to that
        // already queried hit again, the target pos
        // will be the position before that hit
        // (2) The current hit was found at some position
        // before the original jump location.  This can
        // happen if we have an intervening miss (i.e. JE(M))
        // before we find the discordant hit (i.e. JE(HD)).
        // In that case, we leave the target position
        // the same.

        // In either case, we have *already* checked the k-mer
        // at read_target_pos, and it was either a hit (1) or a miss (2).
        // So, we would like to avoid doing that check again.
        // Thus, we traverse up to position (read_target_pos - 1)

        int32_t safe_skip_end = read_target_pos - 1;
        read_prev_pos = kit_tmp->second;
        int32_t dist = (safe_skip_end)-read_prev_pos;

        safe_skip = 1;
        // if the distance is large enough
        if (dist > 2) {
            safe_skip = (safe_skip_end - read_prev_pos) / 2;
            safe_skip = static_cast<int32_t>(safe_skip < 1 ? 1 : safe_skip);
        }
        return dist;
    }

    // Backup kit1 to the position of kit_tmp, the backup
    // position on the read.  Save the current location of
    // kit1 in kit_swap so we can restore it later.
    inline void rewind() {
        // start from the next position on the read and
        // walk until we hit the read end or the position
        // that we wanted to skip to
        kit_swap = kit1;
        kit1 = kit_tmp;
        kit_tmp = kit_swap;
    }

    // Advance kit1 and kit_tmp to (at least) the position following
    // read_target_pos.  This function is intended to be called at
    // the end of a `safe_skip` walk when the query at the original
    // read_target_pos yielded a miss.
    // Side effect : Resets the miss_it counter to 0.
    inline void advance_to_target_pos_successor() {
        int32_t diff = (read_target_pos - kit1->second) + 1;
        kit1 += diff;
        kit_tmp = kit1;
    }

    // Reset the miss counter to 0.
    inline void clear_miss_counter() { miss_it = 0; }

    // Reset kit1 to the place it was at when
    // `rewind()` was called.
    inline void fast_forward() {
        // NOTE: normally we would assume that
        // kit1->second >= kit_tmp->second
        // but, this can be violated if we skipped farther than we asked in the iterator
        // because there were 'N's in the read.
        // in that case case, ignore the fact that fast forward
        // could be a rewind.  We always want to skip back to
        // where we were before.  However, in this case we
        // ASSUME that no hit *after* kit_tmp has been added
        // to the set of collected hits.
        kit1 = kit_tmp;
    }

    /**
     * Advance by the safe skip amount calculated in
     * `compute_safe_skip`.
     * Returns true if we are still before or at the expected
     * end point, and false otherwise.
     */
    inline bool advance_safe() {
        int32_t safe_skip_end = read_target_pos - 1;
        int32_t remaining_len = safe_skip_end - kit1->second;
        safe_skip = (remaining_len <= safe_skip) ? remaining_len : safe_skip;
        safe_skip = (safe_skip < 1) ? 1 : safe_skip;
        kit1 += safe_skip;
        return (kit1->second <= safe_skip_end) and (kit1 != kit_end);
    }

    /**
     * Given that we have obtained a hit in response to a query,
     * and given that that hit was not unexpected (i.e. the jump
     * that lead to it was not of type JE(HD)), then this function
     * will:
     * (1) calculate the next position to be checked.
     * (2) advance kit1 to that position
     * (3) preserve the current hit iterator state in kit_tmp
     * (4) set up the appropriate context for a fast hit check if appropriate
     */
    inline void advance_from_hit(bool switching_unitig = false) {
        int32_t skip = 1;
        // the offset of the hit on the read
        int32_t read_offset = kit1->second;
        // the skip that would take us to the last queryable position of the read
        int32_t read_skip = (read_len - k) - read_offset;

        if (switching_unitig) {
          skip = 1;
          kit1 += skip;
          kit_tmp = kit1;
          clear_miss_counter();
          clear_expectation();
          return;
        }

        // If we got here after a miss, and we have an expectation, and this hit
        // matches it, then we want to avoid replaying this whole scenario again.
        // Because this hit is on the same unitig as the one that eventually led
        // us here, we will compute the same attempted skip.  In this case, we
        // found what we expected, but after 1 or more misses.  In that case,
        // we're satisfied with the rest of this unitig so we simply move on
        // to the position after our original skip position.
        if ( 
            ((miss_it > 0) and (expected_cid != invalid_cid) and
             (expected_cid == phits.contigIdx_))) {
            skip = read_target_pos - read_offset + 1;
            skip = (skip < 1) ? 1 : skip;
            kit1 += skip;
            kit_tmp = kit1;
            clear_miss_counter();
            clear_expectation();
            return;
        }

        // any valid skip should always be at least 1 base
        read_skip = (read_skip < 1) ? 1 : read_skip;

        size_t cStartPos = phits.globalPos_ - phits.contigPos_;
        size_t cEndPos = cStartPos + phits.contigLen_;
        size_t cCurrPos = phits.globalPos_;
        global_contig_pos = static_cast<int64_t>(cCurrPos);

        int32_t ctg_skip = 1;
        // fw ori
        if (phits.contigOrientation_) {
            ctg_skip = static_cast<int64_t>(cEndPos) - (static_cast<int64_t>(cCurrPos + k));
        } else {  // rc ori
            ctg_skip = static_cast<int32_t>(phits.contigPos_);
        }
        // we're already at the end of the contig
        bool at_contig_end = (ctg_skip == 0);
        // if we're at the end of the contig
        // we'll set the contig skip to be one
        // (i.e. look for the next k-mer), but we won't
        // set an expectation on what that contig should be
        if (at_contig_end) {
            ctg_skip = 1;
            clear_expectation();
        } else {  // otherwise, we expect the next contig to be the same
            expected_cid = phits.contigIdx_;
        }

        // remember where we are coming from
        kit_tmp = kit1;

        // remember the reason we are skipping
        last_skip_type = (read_skip <= ctg_skip) ? LastSkipType::SKIP_READ : LastSkipType::SKIP_UNI;

        // skip will be min of read and contig
        skip = (last_skip_type == LastSkipType::SKIP_READ) ? read_skip : ctg_skip;

        // skip must be at least 1
        skip = skip < 1 ? 1 : skip;
        // onward
        int32_t pos_before_skip = kit1->second;
        kit1 += skip;

        read_target_pos = (kit1->second > read_len - k) ? read_len - k : kit1->second;
        clear_miss_counter();

        // was the skip we got the one we expected?
        // that is, was the intervening part of the read
        // free of `N` characters?
        bool expected_skip = (kit1->second - pos_before_skip) == skip;

        // if we didn't get the expected skip, then un-set our expectation
        // about where we will land.
        // NOTE : This should be relatively infrequent, but check the effect
        if (!expected_skip) { clear_expectation(); }

        // if we got the skip we expected, then
        // set ourselves up for a fast check in case we see
        // what we expect to see.
        if (expected_skip and (expected_cid != invalid_cid) and (kit1 != kit_end)) {
            /*
            if (2*cCurrPos > ref_contig_it.size()) {
              std::cout << "cCurrPos = " << 2*cCurrPos << ", ref_contig_len = " <<
            ref_contig_it.size() << "\n"; std::cout << "expected_cid = " << expected_cid << ", skip
            = " << skip << "\n";
            }
            */

            // this will be eligible for a fast hit check
            // and if we get one it will be purely confirmatory
            fast_hit.is_confirmatory = true;
            fast_hit.valid(true);

            if (phits.contigOrientation_) {
                // if match is fw, go to the next k-mer in the contig
                cCurrPos += skip;
                fast_hit.offset = skip;
                fast_hit.expected_match_type = KmerMatchType::IDENTITY_MATCH;
            } else {
                cCurrPos -= skip;
                fast_hit.offset = -skip;
                fast_hit.expected_match_type = KmerMatchType::TWIN_MATCH;
            }

            // set the ref contig iterator position and read off
            // the reference k-mer
            ref_contig_it.at(2 * cCurrPos);
            fast_hit.ref_kmer = ref_contig_it.read(2 * k);
        }
    }

    inline uint32_t alternative_skip_amount() const { return hit_found ? alt_skip : 1; }

    inline void advance_from_miss() {
        int32_t skip = 1;

        // distance from backup position
        int32_t dist = read_target_pos - kit_tmp->second;

        switch (last_skip_type) {
            // we could have not yet seen a hit
            // should move alt_skip at a time
            case LastSkipType::NO_HIT: {
                // int32_t dist_to_end = (read_len - (kit1->second + k));
                kit1 += alternative_skip_amount();
                return;
            } break;

            // we could have seen a hit, and tried to jump to the end of the read / uni
            // and the previous search was either that hit, or a miss
            case LastSkipType::SKIP_READ:
            case LastSkipType::SKIP_UNI: {
                // if distance from backup is 1, that means we already tried
                // last position, so we should just move to the next position
                int32_t pos_before_skip = kit_tmp->second;
                if (dist <= 1) {
                    // if we're past the position we tried to jump to
                    // then we're willing to skip a little faster

                    // debug print
                    // std::cerr << "dist = " << dist << ", read_target_pos = " << read_target_pos
                    // << ", kit_tmp->second" << kit_tmp->second << "\n";

                    // skip goes directly to 2 b/c if dist == 1 skip 1 gets us to a position
                    // we already evaluated. Likewise, if dist == 0 (shouldn't happen)
                    // we need to move to the next un-searched position.
                    skip = (dist == 0) ? 1 : 2;  // (dist < 0) ? 2 : 1;
                    kit_tmp += skip;
                    kit1 = kit_tmp;
                } else {
                    // otherwise move the backup position toward us
                    // and move the current point to the backup
                    skip = (miss_it < 4) ? dist / 2 : 2;
                    skip = (skip < 2) ? 2 : skip;
                    kit_tmp += skip;
                    kit1 = kit_tmp;
                }

                // If we can advance to the next query
                // point without having to consult the index
                // then do it.
                if (  // if we have an expecation
                    (expected_cid != invalid_cid) and
                    // and the prev hit was on the expected unitig
                    (phits.contigIdx_ == expected_cid) and
                    // and we are still before the final target position
                    (kit1->second < read_target_pos)) {
                    // Here, we compute the actual amount we skipped by, since
                    // there may have been intervening 'N's in the read and
                    // we could have skipped more than the expected amount.
                    int actual_skip = (kit1->second - pos_before_skip);

                    // If this was the first miss we encountered
                    // in the most recent round of search, then
                    // reset the fast_hit offset.
                    if (miss_it == 0) { fast_hit.offset = 0; }

                    if (phits.contigOrientation_) {
                        fast_hit.offset += actual_skip;
                        global_contig_pos += actual_skip;
                        fast_hit.expected_match_type = KmerMatchType::IDENTITY_MATCH;
                    } else {
                        fast_hit.offset -= actual_skip;
                        global_contig_pos -= actual_skip;
                        fast_hit.expected_match_type = KmerMatchType::TWIN_MATCH;
                    }
                    fast_hit.is_confirmatory = false;
                    fast_hit.valid(true);
                    ref_contig_it.at(2 * global_contig_pos);
                    fast_hit.ref_kmer = ref_contig_it.read(2 * k);
                }
                // if we pass the read target position, then
                // we no longer have an expectation of what
                // we should see.
                if (kit1->second >= read_target_pos) {
                    //fast_hit.valid(false);
                    clear_expectation();
                }
                // increment the miss iterator
                miss_it += 1;
                return;
            } break;
        }
        // we could have previously not seen a hit either
        // that doesn't matter since it is recursively one of the above cases.
    }

    pufferfish::CanonicalKmerIterator kit1;
    pufferfish::CanonicalKmerIterator kit_tmp;
    pufferfish::CanonicalKmerIterator kit_end;
    pufferfish::CanonicalKmerIterator kit_swap;
    reference_index* pfi = {nullptr};
    sshash::bit_vector_iterator ref_contig_it;
    int32_t read_len;
    int32_t read_target_pos;
    int32_t read_current_pos;
    int32_t read_prev_pos;
    int32_t safe_skip;
    int32_t k;
    FastHitInfo fast_hit;
    uint32_t expected_cid;
    LastSkipType last_skip_type{LastSkipType::NO_HIT};
    int32_t miss_it;
    int64_t global_contig_pos;
    projected_hits phits;
    static constexpr uint32_t invalid_cid{std::numeric_limits<uint32_t>::max()};
    // the amount by which we skip on a missed lookup.
    uint32_t alt_skip;
    // this is false until we find the first hit on the read
    // and is subsequently true. It is used to determine
    // the amount by which we should move forward on a miss
    // (1 if no hit is found yet, else altSkip).
    bool hit_found;
};

struct SkipInfoT {
  int direction{0}; // 1 for forward, -1 for rc
  int32_t offset{0}; // amount to move
}; 
// lightweight struct with just the info necessary 
// for a single skip.

//
// Walk safely (one k-mer at a time until skip_ctx reaches the end_readPos)
//
inline void walk_safely_until(
  SkipContext& skip_ctx,  // the skip context we'll use for searching
  piscem::streaming_query& qc, 
  int end_read_pos, // the position that we will shearch until
  std::vector<std::pair<int, projected_hits>>& raw_hits // the structure where we'll aggregate results
) {
  const int32_t k = skip_ctx.get_k();
  int64_t dist_to_contig_end = 0;
  // while we have not yet reached the last k-mer 
  // of `read`
  while (!skip_ctx.is_exhausted() && skip_ctx.read_pos() <= end_read_pos) {
    // check if the search was a hit
    if (skip_ctx.is_good_query() and skip_ctx.query_kmer(qc)) {
      // in this branch of the if/else, we found a hit 
      // for the current k-mer.

      // record some relevant information;
      // the position on the read of the matching k-mer
      auto read_pos = skip_ctx.read_pos();
      auto initial_search_pos = read_pos;
      // the projected hits object for this hit
      // which includes the contig id and position
      // of the matching k-mer, as well as the orientation
      auto phit_info = skip_ctx.proj_hits();
     
      // compute the relevant information about this hit
      // where the contig starts and ends, and where our 
      // hit is on the contig.
      size_t cStartPos = phit_info.globalPos_ - phit_info.contigPos_;
      size_t cEndPos = cStartPos + phit_info.contigLen_;
      int64_t cCurrPos = static_cast<int64_t>(phit_info.globalPos_);
      
      // determine if we should add this hit. If there are currently 
      // no hits for this read, or if this hit reaches further than 
      // any other (i.e. it's read position is greater than any hit
      // we have seen before), then add it.
      if (raw_hits.empty() or (read_pos > raw_hits.back().first)) {
        auto proj_hits = skip_ctx.proj_hits();
        // this hit was not looking for a match on a known contig
        // rather, it resulted from an "open" search.
        proj_hits.resulted_from_open_search = true;
        raw_hits.push_back({read_pos, proj_hits});
      }

      // determine the distance to the end of the contig 
      // and try to walk the k-mers from the current position
      // to the contig end.
      int32_t direction = 1;
      // fw ori
      if (phit_info.contigOrientation_) {
        dist_to_contig_end = static_cast<int64_t>(cEndPos) - (static_cast<int64_t>(cCurrPos + k));
      } else {  // rc ori
        dist_to_contig_end = static_cast<int64_t>(phit_info.contigPos_);
        direction = -1;
      }

      // otherwise, take the careful path
      bool matches = true;
      bool ended_on_match = false;
      // we'll use this to hold the last valid k-mer match that was 
      // observered
      std::pair<int, projected_hits> last_valid_hit = raw_hits.back();
      // now, we will walk starting from the hit we just saw above 
      // until we reach the end of the read, the end of the contig, or a mismatch.
      while (!skip_ctx.is_exhausted() and matches and dist_to_contig_end > 0) {
        // increment the read iterator by the minimal amount possible.
        int inc_amt = skip_ctx.increment_read_iter();
        dist_to_contig_end -= inc_amt;
        // check to make sure this last move didn't overshoot
        // the end of the contig.
        if (dist_to_contig_end >= 0) {
          int32_t inc_offset = (direction * inc_amt);
          cCurrPos += inc_offset;
          assert(("cCurrPos >= 0", cCurrPos >= 0));

          // set the ref contig iterator position and read off
          // the reference k-mer
          skip_ctx.ref_contig_it.at(2 * cCurrPos);
          skip_ctx.fast_hit.ref_kmer = skip_ctx.ref_contig_it.read(2 * k);
          auto match_type = skip_ctx.check_match();
          matches = (match_type != KmerMatchType::NO_MATCH);

          // if this was a match then store the hit in case we don't 
          // see another before a missing kmer, mismatch, or the end 
          // of the contig.
          if (matches) {
            bool hit_fw = (match_type == KmerMatchType::IDENTITY_MATCH);
            // update the last_valid_hit variable with the 
            // projected hit corresponding to this matching k-mer
            auto& phit = last_valid_hit.second;
            phit.resulted_from_open_search = false;
            phit.contigOrientation_ = hit_fw;
            phit.globalPos_ += inc_offset;
            phit.contigPos_ += inc_offset;
            // set the read position for this hit
            last_valid_hit.first = skip_ctx.read_pos();

            // if this was a bookending k-mer for a contig then record 
            // that we eneded on a match.
            ended_on_match = (dist_to_contig_end == 0);
          } else {
            // the last k-mer we looked for was 
            // not a match.
            break;
          }

        } else {
          // we overshot the end of the contig
          matches = false;
        }
      }
 
      // determine if we had a valid match after the first one 
      // that started this interval. If we did, then add it
      if (last_valid_hit.first > raw_hits.back().first) {
        raw_hits.push_back( last_valid_hit );
        assert(("read_pos > last_valid_hit.first", skip_ctx.read_pos() > last_valid_hit.first));
      } 

     
      // if we ended the above while loop on a match to the 
      // end of a unitig, then we need to increment the skip 
      // context here to get to the next available k-mer.
      // Likewise, if we didn't advance the iterator at all
      // since the initial open search, then we need to do it
      // here.
      if (ended_on_match or (skip_ctx.read_pos() == initial_search_pos)) {
        skip_ctx.increment_read_iter();
      } 
    } else { // otherwise this open k-mer search resulted in a miss. 
      skip_ctx.increment_read_iter();
    }
  }
}


// move forward on the current unitig by the proposed amount
// and check if we have a match. returns true if we have a match
// and false if we don't. The template parameter controls,
// if there is a match, whether we add it to raw hits or 
// roll back the state of skip_ctx.
template <bool add_hit_if_successful> 
inline bool check_direct_match(
			SkipContext& skip_ctx,
			int32_t k,
			int direction,
			int32_t dist,
			int64_t& curr_pos,
			std::vector<std::pair<int, projected_hits>>& raw_hits) 
{

  int32_t inc_offset = (direction * dist);
  curr_pos += inc_offset;
  skip_ctx.ref_contig_it.at(2 * curr_pos);
  skip_ctx.fast_hit.ref_kmer = skip_ctx.ref_contig_it.read(2 * k);

  auto direct_phit = raw_hits.back().second;
  auto prev_hit_fw = direct_phit.hit_fw_on_contig();

  if constexpr (add_hit_if_successful) {
    direct_phit.resulted_from_open_search = false;
    direct_phit.globalPos_ += inc_offset;
    direct_phit.contigPos_ += inc_offset;
  }

  // check if what we find at the given position is a match 
  // or not.
  auto match_type = skip_ctx.check_match();
  bool matches = (match_type != KmerMatchType::NO_MATCH);
  int read_pos = skip_ctx.read_pos();

  // in this branch, we moved forward on the contig and found 
  // a match.
  bool hit_fw = (match_type == KmerMatchType::IDENTITY_MATCH);
  if (matches and (hit_fw == prev_hit_fw)) {
    if constexpr (add_hit_if_successful) {
      direct_phit.contigOrientation_ = hit_fw;
      raw_hits.push_back({read_pos, direct_phit});
      skip_ctx.increment_read_iter();
    } 
    // otherwise, we found the match, but we're not adding it.
    // let the caller know we were successful.
    return true;
  } else {
    // not a match
    return false;
  }
}

// This method performs k-mer / hit collection
// using a custom implementation of the corresponding

struct EveryKmer {
  int64_t cStartPos;
  int64_t cEndPos;
  int64_t cCurrPos;
  int32_t direction;
  int64_t dist_to_contig_end;
  int32_t k;
  bool new_state;

  EveryKmer(int32_t k): 
    cStartPos(std::numeric_limits<int64_t>::max()),
    cEndPos(std::numeric_limits<int64_t>::max()),
    cCurrPos(std::numeric_limits<int64_t>::max()),
    direction(std::numeric_limits<int32_t>::max()),
    dist_to_contig_end(std::numeric_limits<int64_t>::max()),
    k(k), new_state(true)
  {

  }
  EveryKmer(projected_hits ph, size_t k) :
    cStartPos(static_cast<int64_t>(ph.globalPos_ - ph.contigPos_)),
    cEndPos(static_cast<int64_t>(ph.globalPos_ - ph.contigPos_ + ph.contigLen_)),
    cCurrPos(static_cast<int64_t>(ph.globalPos_))
    {
      if (ph.contigOrientation_) {
        direction = 1;
        dist_to_contig_end = static_cast<int64_t>(cEndPos) - (static_cast<int64_t>(cCurrPos + k));
      }
      else {
        direction = -1;
        dist_to_contig_end = static_cast<int64_t>(ph.contigPos_);
      }
      k = k;
      new_state = (dist_to_contig_end <= 0) ? true : false;
    }
  
  inline void set_state(projected_hits &ph) {
    cStartPos = static_cast<int64_t>(ph.globalPos_ - ph.contigPos_);
    cEndPos = static_cast<int64_t>(cStartPos + ph.contigLen_);
    cCurrPos = static_cast<int64_t>(ph.globalPos_);
    if (ph.contigOrientation_) {
      direction = 1;
      dist_to_contig_end = static_cast<int64_t>(cEndPos) - (static_cast<int64_t>(cCurrPos + k));
    } else {
      direction = -1;
      dist_to_contig_end = static_cast<int64_t>(ph.contigPos_);
    }
    new_state = (dist_to_contig_end <= 0) ? true : false;
  }

  inline void query_kmer(pufferfish::CanonicalKmerIterator& kit,
  mindex::reference_index *pfi, std::vector<std::pair<int, projected_hits>> &raw_hits,
  sshash::bit_vector_iterator &ref_contig_it,
  piscem::streaming_query& qc) {
    (void)ref_contig_it;
    qc.reset_state();
    auto ph = pfi->query(kit, qc);
    
    if (!ph.empty()) {
        // std::cerr << "found hit between " << kit->second << " and " << ph.globalPos_ << "\n";
        raw_hits.push_back(std::make_pair(kit->second, ph));
        set_state(ph);
        // std::cerr << "direction " << direction << ", dist_to_end" << dist_to_contig_end << ", cCurrPos: " << cCurrPos << "\n";
        ref_contig_it.at(2*ph.globalPos_);
        uint64_t ref_kmer = ref_contig_it.read(2*k);
        (void)ref_kmer;
        // std::cerr << "\t(MATCH IS) ref_kmer = " << sshash::util::uint_kmer_to_string(ref_kmer, k) << ", kit = " << kit->first.to_str() << "\n";
    }
    else {
      new_state = true;
    }
  }

  inline bool check_match(sshash::bit_vector_iterator &ref_contig_it,
                  pufferfish::CanonicalKmerIterator& kit) {
    
    int64_t cpos = cCurrPos + direction;
    
    ref_contig_it.at(2*cpos);
    auto ref_kmer = ref_contig_it.read(2 * k);
	
    /*
    std::cerr << "\t(k = " << k << ") checking hit between " << kit->second << " and " << cpos << "\n";
    auto k2 = kit->first;
    k2.fromNum(ref_kmer);
    std::cerr << "\t\t ref_kmer = " << k2.to_str() << ", kit = " << kit->first.to_str() << "\n";
    */
    
    auto match_type = kit->first.isEquivalent(ref_kmer);
    bool matches = (match_type != KmerMatchType::NO_MATCH);
    
    return matches;
  }

  inline void add_next(std::vector<std::pair<int, projected_hits>> &raw_hits,
    pufferfish::CanonicalKmerIterator& kit) {
    auto ph = raw_hits.back().second;
    ph.globalPos_ += direction;
    ph.contigPos_ += direction;
    raw_hits.push_back(std::make_pair(kit->second, ph));
    set_state(ph);
  }
};


// This method performs k-mer / hit collection 
// using a custom implementation of the corresponding 
// part of the pseudoalignment algorithm as described in (1).
// Specifically, it attempts to find a small set of
// k-mers along the fragment (`read`) that are shared with
// a set of unitigs in the compacted colored de Bruijn graph,
// using the structure of the graph to avoid queries that
// are unlikely to change the resulting set of unitigs
// that are discovered.  This function only fills out the
// set of hits, and does not itself implement any
// consensus mechanism.  This hit collection mechanism
// prioritizes speed compared to e.g. the uniMEM collection
// strategy implemented in the `operator()` method of this
// class.  Currently, this strategy is only used in the
// `--sketch` mode of alevin.  One may refer to the
// [release notes](https://github.com/COMBINE-lab/salmon/releases/tag/v1.4.0)
// of the relevant version of salmon and links therein
// for a more detailed discussion of the downstream effects of
// different strategies.
//
// The function fills out the appropriate raw hits
// member of this MemCollector instance.
//
// If the `isLeft` flag is set to true, then the left
// raw hits are filled and can be retreived with
// `get_left_hits()`.  Otherwise, the right raw hits are
// filled and can be retrieved with `get_right_hits()`.
//
// This function returns `true` if at least one hit was
// found shared between the read and reference and
// `false` otherwise.
//
// [1] Bray NL, Pimentel H, Melsted P, Pachter L.
// Near-optimal probabilistic RNA-seq quantification.
// Nat Biotechnol. 2016;34(5):525-527.
bool hit_searcher::get_raw_hits_sketch(std::string& read,
                                       piscem::streaming_query& qc, 
                                       mindex::SkippingStrategy strat,
                                       bool isLeft,
                                       bool verbose) {
  (void)verbose;
  bool strict_mode = (strat == mindex::SkippingStrategy::STRICT);

  projected_hits phits;
  auto& raw_hits = isLeft ? left_rawHits : right_rawHits;

  CanonicalKmer::k(k);
  int32_t k = static_cast<int32_t>(CanonicalKmer::k());
  SkipContext skip_ctx(read, pfi_, k, altSkip);
  skip_ctx.fast_hit.valid(false);

  // for this new read, restart the streaming query
  qc.reset_state();

  if (strict_mode) {
    int read_end_pos = read.length() - k;
    walk_safely_until(skip_ctx, qc, read_end_pos, raw_hits);
  } else {
    constexpr bool verbose = false;
    int64_t dist_to_contig_end = 0;
    // while we have not yet reached the last k-mer 
    // of `read`
    while (!skip_ctx.is_exhausted()) {
      // check if the search was a hit
      if (skip_ctx.is_good_query() and skip_ctx.query_kmer(qc)) {
        // in this branch of the if/else, we found a hit 
        // for the current k-mer.

        // record some relevant information;
        // the position on the read of the matching k-mer
        auto read_pos = skip_ctx.read_pos();
        
        // the projected hits object for this hit
        // which includes the contig id and position
        // of the matching k-mer, as well as the orientation
        auto phit = skip_ctx.proj_hits();

        // compute the relevant information about this hit
        // where the contig starts and ends, and where our 
        // hit is on the contig.
        size_t cStartPos = phit.globalPos_ - phit.contigPos_;
        size_t cEndPos = cStartPos + phit.contigLen_;
        int64_t cCurrPos = static_cast<int64_t>(phit.globalPos_);

        // determine if we should add this hit. 
        // If there are currently no hits for this read, or if this 
        // hit reaches further than any other 
        // (i.e. it's read position is greater than any hit we have seen before), 
        // then add it.
        if (raw_hits.empty() or (read_pos > raw_hits.back().first)) {
          auto proj_hits = skip_ctx.proj_hits();
          // this hit was not looking for a match on a known contig
          // rather, it resulted from an "open" search.
          proj_hits.resulted_from_open_search = true;
          raw_hits.push_back({read_pos, proj_hits});
        }

        // determine the distance to the end of the contig 
        // and try to walk the k-mers from the current position
        // to the contig end.
        int32_t direction = 1;
        // fw ori
        if (phit.contigOrientation_) {
          dist_to_contig_end = static_cast<int64_t>(cEndPos) - (static_cast<int64_t>(cCurrPos + k));
        } else {  // rc ori
          dist_to_contig_end = static_cast<int64_t>(phit.contigPos_);
          direction = -1;
        }

        // If we can successfully jump to the end of this unitig then just do it.
        int64_t dist_to_read_end = (read.size() - k) - read_pos;
        int32_t skip_dist = static_cast<uint32_t>(std::min(dist_to_read_end, dist_to_contig_end));
        // if we are already at the end of the read or the unitig, then there is 
        // nothing to do here.
        if (skip_dist > 1) {
          // before we attempt the skip, backup our iterator and our 
          // current contig position.
          auto backup_kit = skip_ctx.get_iter();
          auto backup_cpos = cCurrPos;

          auto neighbor_dist = skip_ctx.increment_read_iter();
          // if we jumped past the end by trying to walk
          // one nucleotide. We're done
          if (skip_ctx.is_exhausted()) { 
            continue;
          }

          if (neighbor_dist < skip_dist) {
            bool found_match = check_direct_match<false>(
              skip_ctx, k, direction, neighbor_dist, 
              cCurrPos, raw_hits);
            if (!found_match) {
              // go to the top of the loop and do a regular search
              continue;
            } else { 
              // found a match, so backup the skip_ctx and 
              // current unitig position to prepare for the 
              // actual skip.
              skip_ctx.set_iter(backup_kit);
              cCurrPos = backup_cpos;
            }
          }

          // try to move forward the expected amount.
          auto actual_dist = skip_ctx.advance_read_iter(skip_dist);
          if constexpr (verbose) {
            std::cerr << "skip_dist: " << skip_dist << ", actual_skip: " << actual_dist << "\n";
          }
          // if we jumped past the end
          if (skip_ctx.is_exhausted()) { 
            // just give up.
            //break; 
            // try the slow way
            skip_ctx.set_iter(backup_kit);
            skip_ctx.increment_read_iter();
            walk_safely_until(skip_ctx, qc, read_pos + skip_dist, raw_hits);
            continue;
          }

          // if this led to a valid iterator, save it for 
          // later.
          auto alt_kit = skip_ctx.get_iter();
          // save the read position for later as well.
          int32_t next_read_pos = skip_ctx.read_pos();

          // if the skip was the expected amount, then just jump ahead 
          // and check it directly.  If it matches, take it and continue
          // otherwise fallback to a standard search.
          if (actual_dist == skip_dist) {
            bool found_match = check_direct_match<true>(
              skip_ctx, k, direction, skip_dist, 
              cCurrPos, raw_hits);
            if (found_match) { continue; }
          }

          // if we got here, either the iterator advanced too far 
          // (so we can't do a quick_check) or the quick check failed.
          // Either way, we have to do a regular search.
          bool alt_found = skip_ctx.query_kmer(qc);

          if constexpr (verbose) {
            std::cerr << "query at read_pos: " << skip_ctx.read_pos() << " " << (alt_found ? "was" : "was not") << " found in the index\n";
          }

          projected_hits alt_phit; 
          // Otherwise, accept the hit only 
          // if the hit occurs on the same contig (in the same orientation?)
          if (alt_found) {
            auto check_phit = skip_ctx.proj_hits();
            bool accept_hit = (check_phit.contig_id() == phit.contig_id())
              and (check_phit.hit_fw_on_contig() == phit.hit_fw_on_contig())
              and ((direction > 0) ? (check_phit.contig_pos() > phit.contig_pos()) 
              : (check_phit.contig_pos() < phit.contig_pos()));
            alt_phit = check_phit;

            // if the hit matches our expectation, accept it
            // here and return to the top of the search loop
            if (accept_hit) { 
              read_pos = skip_ctx.read_pos();
              phit = check_phit; 
              phit.resulted_from_open_search = false;
              raw_hits.push_back({read_pos, phit}); 
              skip_ctx.increment_read_iter();
              continue;
            }
          }  
          // if the hit doesn't match our expectation, or was 
          // empty, then we end up here.

          // we got here and we found a hit for our jump position that 
          // did *not* land on the expected contig (or was not present at all).
          // In this case, check the center k-mer.
          bool mid_acceptable = false;
          if (skip_dist > 4) {
            skip_dist = skip_dist / 2;
            // backup to the original iterator
            skip_ctx.set_iter(backup_kit);
            // move forward to the target.
            skip_ctx.advance_read_iter(skip_dist);

            // if we didn't jump past the end
            if (!skip_ctx.is_exhausted()) { 
              bool mid_found = skip_ctx.query_kmer(qc);
              if (mid_found) {
                auto mid_phit = skip_ctx.proj_hits();
                if (mid_phit.contig_id() == phit.contig_id()) {
                  // matched our first contig
                  // so in this case add the mid_phit and *if it exists* 
                  // the alt_phit
                  mid_phit.resulted_from_open_search = false;
                  raw_hits.push_back({skip_ctx.read_pos(), mid_phit});
                  if (alt_found) {
                    alt_phit.resulted_from_open_search = true;
                    raw_hits.push_back({alt_kit->second, alt_phit});
                  }
                  mid_acceptable = true;
                } else if (alt_found and mid_phit.contig_id() == alt_phit.contig_id()) {
                  // matched our second contig
                  // in this case, no point in adding the mid hit (it's redundant), 
                  // just add the alt_phit.
                  alt_phit.resulted_from_open_search = true;
                  raw_hits.push_back({alt_kit->second, alt_phit});
                  mid_acceptable = true;
                } else {
                  mid_acceptable = false;
                }
              }
            }

          }

          if (mid_acceptable){
            // we will have already added the actual hit 
            // above. Here we set the skip context iterator
            // to our prospective jump position and continue
            // from there.
            
            // NOTE: Consider the implication of the ++ here,
            // this means we will skip right past the inital 
            // jump point (which we did search for), even 
            // if the match we had was to the center k-mer.
            ++alt_kit;

            skip_ctx.set_iter(alt_kit);
            continue;
          } else {
            // if we didn't find acceptable middle 
            // hit, or the skip_distance was <= 4.
            skip_ctx.set_iter(backup_kit);
            // twice because we already did the neighbor 
            // check above, and we got here so it must 
            // have passed.
            skip_ctx.increment_read_iter();
            skip_ctx.increment_read_iter();
            walk_safely_until(skip_ctx, qc, next_read_pos, raw_hits);
            continue;
          }
        } 
        skip_ctx.increment_read_iter();

      } else { // otherwise this open k-mer search resulted in a miss. 
        skip_ctx.increment_read_iter();
      }
    }
  }  
  return !raw_hits.empty();
}




// This method performs k-mer / hit collection
// using a custom implementation of the corresponding
// part of the pseudoalignment algorithm as described in (1).
// Specifically, it attempts to find a small set of
// k-mers along the fragment (`read`) that are shared with
// a set of unitigs in the compacted colored de Bruijn graph,
// using the structure of the graph to avoid queries that
// are unlikely to change the resulting set of unitigs
// that are discovered.  This function only fills out the
// set of hits, and does not itself implement any
// consensus mechanism.  This hit collection mechanism
// prioritizes speed compared to e.g. the uniMEM collection
// strategy implemented in the `operator()` method of this
// class.  Currently, this strategy is only used in the
// `--sketch` mode of alevin.  One may refer to the
// [release notes](https://github.com/COMBINE-lab/salmon/releases/tag/v1.4.0)
// of the relevant version of salmon and links therein
// for a more detailed discussion of the downstream effects of
// different strategies.
//
// The function fills out the appropriate raw hits
// member of this MemCollector instance.
//
// If the `isLeft` flag is set to true, then the left
// raw hits are filled and can be retreived with
// `get_left_hits()`.  Otherwise, the right raw hits are
// filled and can be retrieved with `get_right_hits()`.
//
// This function returns `true` if at least one hit was
// found shared between the read and reference and
// `false` otherwise.
//
// [1] Bray NL, Pimentel H, Melsted P, Pachter L.
// Near-optimal probabilistic RNA-seq quantification.
// Nat Biotechnol. 2016;34(5):525-527.
//
bool hit_searcher::get_raw_hits_sketch_orig(std::string& read,
                                       piscem::streaming_query& qc,
                                       mindex::SkippingStrategy strat,
                                       bool isLeft,
                                       bool verbose) {
    (void)verbose;
    (void)strat;
    projected_hits phits;
    auto& raw_hits = isLeft ? left_rawHits : right_rawHits;

    CanonicalKmer::k(k);
    int32_t k = static_cast<int32_t>(CanonicalKmer::k());
    SkipContext skip_ctx(read, pfi_, k, altSkip);

    // for this new read, restart the streaming query
    qc.reset_state();

    // while it is possible to search further
    while (!skip_ctx.is_exhausted()) {
        // if we had a hit
        bool confirmatory_fast_hit = false;
        if (skip_ctx.query_kmer_complex(qc, confirmatory_fast_hit)) {
            // record this hit
            int32_t read_pos = skip_ctx.read_pos();
            auto proj_hits = skip_ctx.proj_hits();

            // if the hit was not inline with what we were
            // expecting.
            if (skip_ctx.hit_is_unexpected()) {
                // there are two scenarios now. Either:
                // (1) the hit we found was at the inital place
                // we tried to jump to, but just on the wrong unitig.
                // OR
                // (2) we encountered >= 1 miss in between and so
                // the hit we found that led us here is before the
                // original jump point.
                //
                // In case (1), we want to add the current hit
                // *after* doing safe advances to the end position, and
                // then do a normal `advance_from_hit()`; in case (2)
                // we want to add the current hit *before* doing advances
                // to the end position, and then do a normal `advance_from_miss()`.

                bool hit_at_end = !skip_ctx.hit_is_before_target_pos();

                // we are in case (2)
                if (!hit_at_end) { raw_hits.push_back(std::make_pair(read_pos, proj_hits)); }

                int32_t dist_to_target = skip_ctx.compute_safe_skip();

                // special case is the prev hit was the
                // one right before this, in this case
                // no need to do this work, we already
                // have both hits.  In that case, just
                // proceed as usual.
                if (dist_to_target > 0) {
                    // start from the next position on the read
                    skip_ctx.rewind();

                    // walk until we hit the read end or the position
                    // that we wanted to skip to, collecting the
                    // hits we find along the way.
                    while (skip_ctx.advance_safe()) {
                        bool confirmatory_fast_hit = false;
                        if (skip_ctx.query_kmer_complex(qc, confirmatory_fast_hit)) {
                            raw_hits.push_back(
                                std::make_pair(skip_ctx.read_pos(), skip_ctx.proj_hits()));
                        }
                    }

                    // now we examined the positions in between.
                    // if we were in case (1), jump back to the valid
                    // hit that we first encountered.  If we were in
                    // case (2), advance to at least the point right
                    // after the initial jump.
                    if (hit_at_end) {
                        skip_ctx.fast_forward();
                    } else {
                        skip_ctx.advance_to_target_pos_successor();
                        skip_ctx.clear_miss_counter();
                    }
                }

                // if we got here, then either we are done, or we
                // reached our target position so we don't
                // have an expectation of what should come next.
                skip_ctx.clear_expectation();

                // If we were in case (1)
                if (hit_at_end) {
                    // set the phits for the skip_ctx so the
                    // next jump can be properly computed
                    skip_ctx.phits = proj_hits;
                    // add the hit here
                    raw_hits.push_back(std::make_pair(read_pos, proj_hits));
                    // advance as if this was a normal hit
                    skip_ctx.advance_from_hit(true);
                } else {
                    // We were in case 2
                    // now advance as if this was a normal miss.
                    skip_ctx.advance_from_miss();
                }
            } else {
                // We got a hit and either we had no expectation
                // or the hit was in accordance with our expectation.

                // push this hit and advance
                if (!confirmatory_fast_hit) {
                    raw_hits.push_back(std::make_pair(read_pos, proj_hits));
                }
                skip_ctx.advance_from_hit();
            }
        } else {
            // if we got here, we looked for a match and didn't find one
            skip_ctx.advance_from_miss();
        }
    }

    return raw_hits.size() != 0;
  }
  

bool hit_searcher::get_raw_hits_sketch_everykmer(std::string &read,
                  piscem::streaming_query& qc,
                  bool isLeft,
                  bool verbose) {
    clear();
    (void) verbose;
    auto& raw_hits = isLeft ? left_rawHits : right_rawHits;
    pufferfish::CanonicalKmerIterator kit(read), kit_end;
    CanonicalKmer::k(k);
    int32_t k = static_cast<int32_t>(CanonicalKmer::k());

    // qc.reset_state();
    EveryKmer evs(k);
    
    // auto ref_contig_it = sshash::bit_vector_iterator(pfi_->contigs(), 0);

    // Look at every kmer: if new state, do index query
    // Else: move forward by 1 kmer on reference and see if it matches the contig kmer.
    //       If it matches the kmer on the contig then update proj hit and add it
    //       Otherwise do the index query
    // New state is set to true for the following - 
    //  1) The first kmer of read
    //  2) The immedidate contig reference kmer does not match the read kmer following the previous hit
    //  3) The distance to contig end is 0
    // At each kmer check we are resetting the phit variables, where the new_state is also defined
    // size_t size_hits = 0;
    while(kit != kit_end) {
      auto ph = pfi_->query(kit, qc);
      if (!ph.empty()) {
        // std::cerr << "found hit between " << kit->second << " and " << ph.globalPos_ << "\n";
        raw_hits.push_back(std::make_pair(kit->second, ph));
      }
  //     if (evs.new_state)  {
  //       evs.query_kmer(kit, pfi_, raw_hits, ref_contig_it, qc);
	// // new_state_cnt++;
  //     } else {
  //         bool matches = evs.check_match(ref_contig_it, kit);
  //         if (matches) {
  //           evs.add_next(raw_hits, kit);
	//     // matches_cnt++;
  //         } else {
  //           evs.query_kmer(kit, pfi_, raw_hits, ref_contig_it, qc); 
	//     // non_matches_cnt++;
  //         }
  //     }
      //auto ph = raw_hits.back().second;
      //auto &refs = ph.refRange;
      // std::cout << "refs len " << refs.size() << std::endl;
      // for (auto v : refs) {
      //       const auto &ref_pos_ori = ph.decode_hit(v);
      //       uint32_t tid = sshash::util::transcript_id(v);
      //       int32_t pos = static_cast<int32_t>(ref_pos_ori.pos);
      //       std::cout << "tid " << tid << " pos " << pos << std::endl;
      // }
  
      ++kit;
    }
    return raw_hits.size() != 0;
}

void hit_searcher::clear() {
    left_rawHits.clear();
    right_rawHits.clear();
}

}  // namespace mindex
