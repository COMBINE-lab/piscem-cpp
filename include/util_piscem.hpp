#include "bitsery/bitsery.h"
#include "bitsery/brief_syntax.h"
#include "../external/sshash/include/hash_util.hpp"
#include "../external/sshash/include/util.hpp"
#include "../include/boost/unordered/concurrent_flat_map.hpp"

#pragma once

#include <vector>
#include <cassert>
#include <fstream>
#include <cmath>  // for std::ceil on linux

namespace piscem {
    class unitig_end_cache_t {
    public:
      explicit unitig_end_cache_t(size_t max_size) : 
        m_max_size(max_size), m_unitig_end_map(max_size) {}

      unitig_end_cache_t() = delete;
      unitig_end_cache_t(const unitig_end_cache_t& other) = delete;
      unitig_end_cache_t(unitig_end_cache_t&& other) = delete;
      unitig_end_cache_t& operator=(const unitig_end_cache_t& other) = delete;
      unitig_end_cache_t& operator=(unitig_end_cache_t&& other) = delete;

      inline boost::concurrent_flat_map<uint64_t, sshash::lookup_result>* get_map() {
        return &m_unitig_end_map;
      }
      inline size_t get_capacity() const { return m_max_size; }
      inline size_t get_map_size() const { return m_unitig_end_map.size(); }
    private:
      size_t m_max_size{0};
      boost::concurrent_flat_map<uint64_t, sshash::lookup_result> m_unitig_end_map;
    };
}

namespace sshash {

    namespace constants {
        constexpr uint64_t hashcode_bits = 64;
    }

    struct RobinHoodHash {
  // from https://github.com/martinus/robin-hood-hashing/blob/master/src/include/robin_hood.h
        inline size_t operator()(uint64_t x) const noexcept {
            // tried lots of different hashes, let's stick with murmurhash3. It's simple, fast, well tested,
            // and doesn't need any special 128bit operations.
            x ^= x >> 33U;
            x *= UINT64_C(0xff51afd7ed558ccd);
            x ^= x >> 33U;

            // not doing the final step here, because this will be done by keyToIdx anyways
            // x *= UINT64_C(0xc4ceb9fe1a85ec53);
            // x ^= x >> 33U;
            return static_cast<size_t>(x);
        }
    };
    struct labeled_poison_occ_t {
    // we need to know what unitig this occurs on
    // what position on that unitig
    // and what the k-mer is (maybe also orientation?)
    uint64_t canonical_kmer{std::numeric_limits<uint64_t>::max()};
    uint32_t unitig_id{std::numeric_limits<uint32_t>::max()};
    uint32_t unitig_pos{std::numeric_limits<uint32_t>::max()};
    };

    inline bool operator==(const labeled_poison_occ_t& a, const labeled_poison_occ_t& b) {
        return (a.canonical_kmer == b.canonical_kmer) and 
            (a.unitig_id == b.unitig_id) and
            (a.unitig_pos == b.unitig_pos);
    }

    inline std::ostream& operator<<(std::ostream& os, const labeled_poison_occ_t& a) {
        os << "{ kmer (u64): " << a.canonical_kmer << ", unitig_id (u32): " << a.unitig_id << ", unitig_pos (u32): " << a.unitig_pos << "\n";
        return os;
    }

    // an occurrence of a poison k-mer
    struct poison_occ_t {
    // we need to know what unitig this occurs on
    // what position on that unitig
    // and what the k-mer is (maybe also orientation?)
    uint32_t unitig_id{std::numeric_limits<uint32_t>::max()};
    uint32_t unitig_pos{std::numeric_limits<uint32_t>::max()};

    poison_occ_t() :
        unitig_id(std::numeric_limits<uint32_t>::max()),
        unitig_pos(std::numeric_limits<uint32_t>::max()) {}

    poison_occ_t(const labeled_poison_occ_t& o) : 
        unitig_id(o.unitig_id),
        unitig_pos(o.unitig_pos) {}

    /*poison_occ_t& operator=(const labeled_poison_occ_t& o) :
        unitig_id(o.unitig_id),
        unitig_pos(o.unitig_pos) { return *this; }*/

    poison_occ_t(const poison_occ_t& o) = default;
    poison_occ_t(poison_occ_t&& o) = default;
    poison_occ_t& operator=(const poison_occ_t& o) = default;
    poison_occ_t& operator=(poison_occ_t&& o) = default;

    // for bitsery to allow seralizing this type
    template<typename S>
        void serialize(S& s) { s(unitig_id, unitig_pos); }
    };

    inline bool operator==(const poison_occ_t& a, const poison_occ_t& b) {
        return (a.unitig_id == b.unitig_id) and
            (a.unitig_pos == b.unitig_pos);
    }

    inline std::ostream& operator<<(std::ostream& os, const poison_occ_t& a) {
        os << "{ unitig_id (u32): " << a.unitig_id << ", unitig_pos (u32): " << a.unitig_pos << "\n";
        return os;
    }



namespace util {

    class PiscemIndexUtils {
        public:

        inline static uint64_t ref_shift(uint64_t ref_shift_in) {
            std::swap(_ref_shift, ref_shift_in);
            return ref_shift_in;
        }

        inline static uint64_t ref_shift() { return _ref_shift; }

        inline static uint64_t pos_mask(uint64_t pos_mask_in) {
            std::swap(_pos_mask, pos_mask_in);
            return pos_mask_in;
        }

        inline static uint64_t pos_mask() { return _pos_mask; }
        
        private:
        // the amount we have to shift right
        // to get just the reference on which
        // the hit occurs
        inline static uint64_t _ref_shift;

        // once we shift 1 bit (to get rid of position)
        // the maks we have to apply to remove the 
        // upper reference bits
        inline static uint64_t _pos_mask;
    };

    struct contig_span {
        pthash::compact_vector::iterator start;
        pthash::compact_vector::iterator stop;
        size_t len=0;

        inline pthash::compact_vector::iterator begin() { return start; }
        inline pthash::compact_vector::iterator end() { return stop; }
        inline bool empty() const { return len == 0; }
        inline size_t size() const { return len; }
    };
    struct ec_span {
        pthash::compact_vector::iterator start;
        pthash::compact_vector::iterator stop;
        size_t len=0;

        inline pthash::compact_vector::iterator begin() { return start; }
        inline pthash::compact_vector::iterator end() { return stop; }
        inline bool empty() const { return len == 0; }
        inline size_t size() const { return len; }
    };

    constexpr uint64_t pos_masks[] = {
        0x0, 0x1, 0x3, 0x7, 0xf, 0x1f, 0x3f, 0x7f, 0xff, 0x1ff, 0x3ff,
        0x7ff, 0xfff, 0x1fff, 0x3fff, 0x7fff, 0xffff, 0x1ffff, 0x3ffff,
        0x7ffff, 0xfffff, 0x1fffff, 0x3fffff, 0x7fffff, 0xffffff,
        0x1ffffff, 0x3ffffff, 0x7ffffff, 0xfffffff, 0x1fffffff, 0x3fffffff,
        0x7fffffff, 0xffffffff, 0x1ffffffff, 0x3ffffffff, 0x7ffffffff,
        0xfffffffff, 0x1fffffffff, 0x3fffffffff, 0x7fffffffff, 0xffffffffff,
        0x1ffffffffff, 0x3ffffffffff, 0x7ffffffffff, 0xfffffffffff,
        0x1fffffffffff, 0x3fffffffffff, 0x7fffffffffff, 0xffffffffffff,
        0x1ffffffffffff, 0x3ffffffffffff, 0x7ffffffffffff, 0xfffffffffffff,
        0x1fffffffffffff, 0x3fffffffffffff, 0x7fffffffffffff, 0xffffffffffffff,
        0x1ffffffffffffff, 0x3ffffffffffffff, 0x7ffffffffffffff,
        0xfffffffffffffff, 0x1fffffffffffffff, 0x3fffffffffffffff,
        0x7fffffffffffffff};

    inline uint32_t transcript_id(uint64_t e) { 
        return static_cast<uint32_t>((e >> PiscemIndexUtils::ref_shift() ));
    }
    inline uint32_t pos(uint64_t e) { 
        return static_cast<uint32_t>((e >> 1) & PiscemIndexUtils::pos_mask() ); 
    }
    inline bool orientation(uint64_t e) { return (e & 0x1); }

    inline uint64_t encode_contig_entry(uint64_t refctr, uint64_t current_offset, bool is_fw) {
        // e starts out with the reference index
        uint64_t e = refctr;
        // we shift this left by _ref_shift (which is pos_bits + 1)
        e <<= PiscemIndexUtils::ref_shift();
        // shift the current offset left by 1 and add it to the representation
        e |= (current_offset << 1);
        // set the orientation bit
        e |= is_fw ? 1 : 0;
        return e;
    }

    struct Position {
        uint32_t transcript_id_;  // reference id
        uint32_t pos_;            // position in reference

        // bool orien;
        // Position() {
        //    transcript_id_ = std::numeric_limits<decltype(transcript_id_)>::max();
        //    pos_ = std::numeric_limits<decltype(pos_)>::max();
        //}
        /* can't have this and be POD
        Position(uint32_t tid, uint32_t tpos, bool torien) {
            transcript_id_ = tid;
            pos_ = tpos;
            setOrientation(torien);
            // orien = torien;
        }
        */

        // The most significant bit carry
        // the orientation information
        void setOrientation(bool orientation) {
            if (orientation) {
                pos_ |= 1 << 31;
            } else {
                pos_ &= 0x7FFFFFFF;
            }
        }

        inline uint32_t transcript_id() { return transcript_id_; }
        inline uint32_t pos() { return (pos_ & 0x7FFFFFFF); }
        inline bool orientation() { return (pos_ & 0x80000000); }

        template <typename Visitor>
        void visit(Visitor& visitor) {
            visitor.visit(transcript_id_);
            visitor.visit(pos_);
        }

        void update(uint32_t tid, uint32_t tpos, bool torien) {
            transcript_id_ = tid;
            pos_ = tpos;
            setOrientation(torien);
        }
    };

    static inline void check_hash_collision_probability(uint64_t size) {
        /*
            From: https://preshing.com/20110504/hash-collision-probabilities/
            Given a universe of size U (total number of possible hash values),
            which is U = 2^b for b-bit hash codes,
            the collision probability for n keys is (approximately):
                1 - e^{-n(n-1)/(2U)}.
            For example, for U=2^32 (32-bit hash codes), this probability
            gets to 50% already for n = 77,163 keys.
            We can approximate 1-e^{-X} with X when X is sufficiently small.
            Then our collision probability is
                n(n-1)/(2U) ~ n^2/(2U).
            So it can derived that ~1.97B keys and 64-bit hash codes,
            we have a probability of collision that is ~0.1 (10%), which may not be
            so small for certain applications.
            For n = 2^30, the probability of collision is ~0.031 (3.1%).
        */
        if (constants::hashcode_bits == 64 and size > (1ULL << 30)) {
            throw std::runtime_error(
                "Using 64-bit hash codes with more than 2^30 keys can be dangerous due to "
                "collisions: use 128-bit hash codes instead.");
        }
    }

}
}
