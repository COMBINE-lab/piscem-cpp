#pragma once

#include <vector>
#include <cassert>
#include <fstream>
#include <cmath>  // for std::ceil on linux

#include "../external/pthash/include/pthash.hpp"

namespace sshash {

namespace constants {
constexpr uint64_t max_k = 31;  // max *odd* size that can be packed into 64 bits
constexpr uint64_t invalid_uint64 = uint64_t(-1);
constexpr uint32_t invalid_uint32 = uint32_t(-1);
constexpr uint64_t seed = 1;
constexpr uint64_t hashcode_bits = 64;
constexpr double c = 3.0;  // for PTHash
constexpr uint64_t min_l = 6;
constexpr uint64_t max_l = 12;
static const std::string default_tmp_dirname(".");
constexpr bool forward_orientation = 0;
constexpr bool backward_orientation = 1;
}  // namespace constants

typedef pthash::murmurhash2_64 base_hasher_type;

typedef pthash::single_phf<base_hasher_type,               // base hasher
                           pthash::dictionary_dictionary,  // encoder type
                           true                            // minimal output
                           >
    pthash_mphf_type;


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

// an occurrence of a poison k-mer
struct poison_occ_t {
  // we need to know what unitig this occurs on
  // what position on that unitig
  // and what the k-mer is (maybe also orientation?)
  uint64_t canonical_kmer{std::numeric_limits<uint64_t>::max()};
  uint32_t unitig_id{std::numeric_limits<uint32_t>::max()};
  uint32_t unitig_pos{std::numeric_limits<uint32_t>::max()};
};

inline bool operator==(const poison_occ_t& a, const poison_occ_t& b) {
  return (a.canonical_kmer == b.canonical_kmer) and 
         (a.unitig_id == b.unitig_id) and
         (a.unitig_pos == b.unitig_pos);
}

inline std::ostream& operator<<(std::ostream& os, const poison_occ_t& a) {
  os << "{ kmer (u64): " << a.canonical_kmer << ", unitig_id (u32): " << a.unitig_id << ", unitig_pos (u32): " << a.unitig_pos << "\n";
  return os;
}


struct streaming_query_report {
    streaming_query_report()
        : num_kmers(0), num_positive_kmers(0), num_searches(0), num_extensions(0) {}
    uint64_t num_kmers;
    uint64_t num_positive_kmers;
    uint64_t num_searches;
    uint64_t num_extensions;
};

struct lookup_result {
    lookup_result()
        : kmer_id(constants::invalid_uint64)
        , kmer_id_in_contig(constants::invalid_uint32)
        , kmer_orientation(constants::forward_orientation)
        , contig_id(constants::invalid_uint32)
        , contig_size(constants::invalid_uint32) {}
    uint64_t kmer_id;            // "absolute" kmer-id
    uint32_t kmer_id_in_contig;  // "relative" kmer-id: 0 <= kmer_id_in_contig < contig_size
    uint32_t kmer_orientation;
    uint32_t contig_id;
    uint32_t contig_size;
};

inline std::ostream& operator<<(std::ostream& os, const lookup_result& r) {
  os << " { kmer_id: " << r.kmer_id  << ", "
     << "kmer_id_in_contig: " << r.kmer_id_in_contig 
     << ", kmer_orientation: " << r.kmer_orientation
     << ", contig_id: " << r.contig_id 
     << ", contig_size: " << r.contig_size << "}\n";
  return os;
}

[[maybe_unused]] static bool equal_lookup_result(lookup_result expected, lookup_result got) {
    if (expected.kmer_id != got.kmer_id) {
        std::cout << "expected kmer_id " << expected.kmer_id << " but got " << got.kmer_id
                  << std::endl;
        return false;
    }
    if (expected.kmer_id_in_contig != got.kmer_id_in_contig) {
        std::cout << "expected kmer_id_in_contig " << expected.kmer_id_in_contig << " but got "
                  << got.kmer_id_in_contig << std::endl;
        return false;
    }
    if (got.kmer_id != constants::invalid_uint64 and
        expected.kmer_orientation != got.kmer_orientation) {
        std::cout << "expected kmer_orientation " << expected.kmer_orientation << " but got "
                  << got.kmer_orientation << std::endl;
        return false;
    }
    if (expected.contig_id != got.contig_id) {
        std::cout << "expected contig_id " << expected.contig_id << " but got " << got.contig_id
                  << std::endl;
        return false;
    }
    if (expected.contig_size != got.contig_size) {
        std::cout << "expected contig_size " << expected.contig_size << " but got "
                  << got.contig_size << std::endl;
        return false;
    }
    return true;
}

struct build_configuration {
    build_configuration()
        : k(31)
        , m(17)
        , seed(constants::seed)

        , l(constants::min_l)
        , c(constants::c)

        , canonical_parsing(false)
        , weighted(false)
        , verbose(true)
        , num_threads(1)
        , tmp_dirname(constants::default_tmp_dirname) {}

    uint64_t k;  // kmer size
    uint64_t m;  // minimizer size
    uint64_t seed;

    uint64_t l;  // drive dictionary trade-off
    double c;    // drive PTHash trade-off

    bool canonical_parsing;
    bool weighted;
    bool verbose;
    
    uint64_t num_threads; // number of threads to use during construction
    std::string tmp_dirname;

    void print() const {
        std::cout << "k = " << k << ", m = " << m << ", seed = " << seed << ", l = " << l
                  << ", c = " << c
                  << ", canonical_parsing = " << (canonical_parsing ? "true" : "false")
                  << ", weighted = " << (weighted ? "true" : "false") << std::endl;
    }
};

struct buckets_statistics {
    static const uint64_t max_bucket_size = 4 * 1024;
    static const uint64_t max_string_size = 256;

    buckets_statistics(uint64_t num_buckets, uint64_t num_kmers, uint64_t num_super_kmers = 0)
        : m_num_buckets(num_buckets)
        , m_num_kmers(num_kmers)
        // , m_num_super_kmers(num_super_kmers)
        , m_max_num_kmers_in_super_kmer(0)
        , m_max_num_super_kmers_in_bucket(0) {
        (void)num_super_kmers;
        m_bucket_sizes.resize(max_bucket_size + 1, 0);
        m_total_kmers.resize(max_bucket_size + 1, 0);
        m_string_sizes.resize(max_string_size + 1, 0);
    }

    void add_num_super_kmers_in_bucket(uint64_t num_super_kmers_in_bucket) {
        if (num_super_kmers_in_bucket < max_bucket_size + 1) {
            m_bucket_sizes[num_super_kmers_in_bucket] += 1;
        }
    }

    void add_num_kmers_in_super_kmer(uint64_t num_super_kmers_in_bucket,
                                     uint64_t num_kmers_in_super_kmer) {
        if (num_super_kmers_in_bucket < max_bucket_size + 1) {
            m_total_kmers[num_super_kmers_in_bucket] += num_kmers_in_super_kmer;
        }
        if (num_kmers_in_super_kmer > m_max_num_kmers_in_super_kmer) {
            m_max_num_kmers_in_super_kmer = num_kmers_in_super_kmer;
        }
        if (num_super_kmers_in_bucket > m_max_num_super_kmers_in_bucket) {
            m_max_num_super_kmers_in_bucket = num_super_kmers_in_bucket;
        }
        if (num_kmers_in_super_kmer < max_string_size + 1)
            m_string_sizes[num_kmers_in_super_kmer] += 1;
    }

    uint64_t num_kmers() const { return m_num_kmers; }
    uint64_t num_buckets() const { return m_num_buckets; }
    uint64_t max_num_super_kmers_in_bucket() const { return m_max_num_super_kmers_in_bucket; }

    void print() const {
        // full statistics
        // std::cout << " === bucket statistics === \n";
        // for (uint64_t bucket_size = 1, prev_bucket_size = 0, prev_kmers_in_buckets = 0,
        //               kmers_in_buckets = 0;
        //      bucket_size != max_bucket_size + 1; ++bucket_size) {
        //     if (m_bucket_sizes[bucket_size] > 0) {
        //         std::cout << "buckets with " << bucket_size
        //                   << " super_kmers=" << m_bucket_sizes[bucket_size] << "("
        //                   << (m_bucket_sizes[bucket_size] * 100.0) / m_num_buckets
        //                   << "%)|total_kmers=" << m_total_kmers[bucket_size] << "("
        //                   << (m_total_kmers[bucket_size] * 100.0) / m_num_kmers << "%)"
        //                   << "|avg_num_kmers_per_bucket="
        //                   << static_cast<double>(m_total_kmers[bucket_size]) /
        //                          m_bucket_sizes[bucket_size]
        //                   << "|avg_num_kmers_per_string="
        //                   << static_cast<double>(m_total_kmers[bucket_size]) /
        //                          (m_bucket_sizes[bucket_size] * bucket_size)
        //                   << std::endl;
        //         kmers_in_buckets += m_total_kmers[bucket_size];
        //     }
        //     if (bucket_size == 4 or bucket_size == 8 or bucket_size == 16 or bucket_size == 32 or
        //         bucket_size == 64 or bucket_size == 128 or bucket_size == 256 or
        //         bucket_size == 512 or bucket_size == 1024 or bucket_size == max_bucket_size) {
        //         assert(kmers_in_buckets >= prev_kmers_in_buckets);

        //         std::cout << " *** kmers in buckets of size > " << prev_bucket_size
        //                   << " and <= " << bucket_size << ": "
        //                   << kmers_in_buckets - prev_kmers_in_buckets << "("
        //                   << (100.0 * (kmers_in_buckets - prev_kmers_in_buckets)) / m_num_kmers
        //                   << "%)" << std::endl;
        //         std::cout << " *** kmers in buckets of size <= " << bucket_size << ": "
        //                   << kmers_in_buckets << "(" << (100.0 * kmers_in_buckets) / m_num_kmers
        //                   << "%)" << std::endl;

        //         prev_bucket_size = bucket_size;
        //         prev_kmers_in_buckets = kmers_in_buckets;
        //     }
        // }

        std::cout << " === bucket statistics (less) === \n";
        for (uint64_t bucket_size = 1; bucket_size != 16 + 1; ++bucket_size) {
            if (m_bucket_sizes[bucket_size] > 0) {
                std::cout << "buckets with " << bucket_size << " super_kmers = "
                          << (m_bucket_sizes[bucket_size] * 100.0) / m_num_buckets << "%"
                          << std::endl;
            }
        }
        std::cout << "max_num_super_kmers_in_bucket " << m_max_num_super_kmers_in_bucket
                  << std::endl;

        // std::cout << " === super_kmer statistics === \n";
        // uint64_t total_super_kmers = 0;
        // uint64_t total_kmers = 0;
        // for (uint64_t string_size = 1; string_size != max_string_size + 1; ++string_size) {
        //     if (m_string_sizes[string_size] > 0) {
        //         std::cout << "super_kmers with " << string_size
        //                   << " kmer=" << m_string_sizes[string_size] << "("
        //                   << (m_string_sizes[string_size] * 100.0) / m_num_super_kmers
        //                   << "%)|total_kmers=" << (string_size * m_string_sizes[string_size]) <<
        //                   "("
        //                   << (string_size * m_string_sizes[string_size] * 100.0) / m_num_kmers
        //                   << "%)" << std::endl;
        //         total_super_kmers += m_string_sizes[string_size];
        //         total_kmers += string_size * m_string_sizes[string_size];
        //     }
        // }
        // std::cout << "total_super_kmers " << total_super_kmers << "/" << m_num_super_kmers << "("
        //           << (total_super_kmers * 100.0) / m_num_super_kmers << "%)" << std::endl;
        // std::cout << "total_kmers " << total_kmers << "/" << m_num_kmers << " ("
        //           << (total_kmers * 100.0) / m_num_kmers << "%)" << std::endl;
        // std::cout << "max_num_kmers_in_super_kmer " << m_max_num_kmers_in_super_kmer <<
        // std::endl;
    }

private:
    uint64_t m_num_buckets;
    uint64_t m_num_kmers;
    // uint64_t m_num_super_kmers;
    uint64_t m_max_num_kmers_in_super_kmer;
    uint64_t m_max_num_super_kmers_in_bucket;
    std::vector<uint64_t> m_bucket_sizes;
    std::vector<uint64_t> m_total_kmers;
    std::vector<uint64_t> m_string_sizes;
};

namespace util {

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

// For the time being, assume < 4B contigs
// and that each contig is < 4B bases
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

/* return the position of the most significant bit */
static inline uint32_t msb(uint32_t x) {
    assert(x > 0);
    return 31 - __builtin_clz(x);
}

static inline uint32_t ceil_log2_uint32(uint32_t x) { return (x > 1) ? msb(x - 1) + 1 : 0; }

[[maybe_unused]] static bool ends_with(std::string const& str, std::string const& pattern) {
    if (pattern.size() > str.size()) return false;
    return std::equal(pattern.begin(), pattern.end(), str.end() - pattern.size());
}

// for a sorted list of size n whose universe is u
[[maybe_unused]] static uint64_t elias_fano_bitsize(uint64_t n, uint64_t u) {
    // return n * ((u > n ? (std::ceil(std::log2(static_cast<double>(u) / n))) : 0) + 2);
    uint64_t l = uint64_t((n && u / n) ? pthash::util::msb(u / n) : 0);
    uint64_t high_bits = n + (u >> l) + 1;
    uint64_t low_bits = n * l;
    return high_bits + low_bits;
}

/*
char decimal  binary
 A     65     01000-00-1 -> 00
 C     67     01000-01-1 -> 01
 G     71     01000-11-1 -> 11
 T     84     01010-10-0 -> 10
*/
// static uint64_t char_to_uint64(char c) { return (c >> 1) & 3; }

// static char uint64_to_char(uint64_t x) {
//     assert(x <= 3);
//     static char nucleotides[4] = {'A', 'C', 'T', 'G'};
//     return nucleotides[x];
// }

/*
    Traditional mapping.
*/
static inline uint64_t char_to_uint64(char c) {
    switch (c) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
    }
    assert(false);
    return -1;
}

static inline char uint64_to_char(uint64_t x) {
    switch (x) {
        case 0:
            return 'A';
        case 1:
            return 'C';
        case 2:
            return 'G';
        case 3:
            return 'T';
    }
    assert(false);
    return 0;
}

/****************************************************************************
    The following two functions preserves the lexicographic order of k-mers,
    that is: if g and t are two k-mers and g < t lexicographically,
    then also id(g) < id(t).
*/
[[maybe_unused]] static uint64_t string_to_uint64(char const* str, uint64_t k) {
    assert(k <= 32);
    uint64_t x = 0;
    for (uint64_t i = 0; i != k; ++i) {
        x <<= 2;
        x += char_to_uint64(str[i]);
    }
    return x;
}
[[maybe_unused]] static void uint64_to_string(uint64_t x, char* str, uint64_t k) {
    assert(k <= 32);
    for (int i = k - 1; i >= 0; --i) {
        str[i] = uint64_to_char(x & 3);
        x >>= 2;
    }
}
/****************************************************************************/

[[maybe_unused]] static std::string uint64_to_string(uint64_t x, uint64_t k) {
    assert(k <= 32);
    std::string str;
    str.resize(k);
    uint64_to_string(x, str.data(), k);
    return str;
}

[[maybe_unused]] static uint64_t string_to_uint64_no_reverse(char const* str, uint64_t k) {
    assert(k <= 32);
    uint64_t x = 0;
    for (uint64_t i = 0; i != k; ++i) x += char_to_uint64(str[i]) << (2 * i);
    return x;
}

static void uint64_to_string_no_reverse(uint64_t x, char* str, uint64_t k) {
    assert(k <= 32);
    for (uint64_t i = 0; i != k; ++i) {
        str[i] = uint64_to_char(x & 3);
        x >>= 2;
    }
}

[[maybe_unused]] static std::string uint64_to_string_no_reverse(uint64_t x, uint64_t k) {
    assert(k <= 32);
    std::string str;
    str.resize(k);
    uint64_to_string_no_reverse(x, str.data(), k);
    return str;
}

/*
    taken from Blight:
    it works with the map
    A -> 00; C -> 01; G -> 11; T -> 10
    Example:
    reverse_complement("ACTCACG") = CGTGAGT
    in binary:
    reverse_complement("00011001000111") = 01111011001110
*/
/*
[[maybe_unused]] static uint64_t compute_reverse_complement(uint64_t x, uint64_t size) {
    assert(size <= 32);
    // Complement, swap byte order
    uint64_t res = __builtin_bswap64(x ^ 0xaaaaaaaaaaaaaaaa);
    // Swap nuc order in bytes
    const uint64_t c1 = 0x0f0f0f0f0f0f0f0f;
    const uint64_t c2 = 0x3333333333333333;
    res = ((res & c1) << 4) | ((res & (c1 << 4)) >> 4);  // swap 2-nuc order in bytes
    res = ((res & c2) << 2) | ((res & (c2 << 2)) >> 2);  // swap nuc order in 2-nuc
    // Realign to the right
    res >>= 64 - 2 * size;
    return res;
}
*/

[[maybe_unused]] static uint64_t compute_reverse_complement(uint64_t x, uint64_t k) {
    uint64_t res = ~x;
    res = (res >> 2 & 0x3333333333333333) | (res & 0x3333333333333333) << 2;
    res = (res >> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) << 4;
    res = (res >> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) << 8;
    res = (res >> 16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16;
    res = (res >> 32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32;

    return res >> (2 * (32 - k));
}


// forward character map. A -> A, C -> C, G -> G, T -> T. rest maps to zero.
static const char canonicalize_basepair_forward_map[256] = {
    0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,
    0, 0, 0, 0, 0, 0, 0, 65, 0, 67, 0, 0, 0, 71, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 84, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};

// reverse character map. A -> T, C -> G, G -> C, T -> A. rest maps to zero.
static const char canonicalize_basepair_reverse_map[256] = {
    0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,
    0, 0, 0, 0, 0, 0, 0, 84, 0, 71, 0, 0, 0, 67, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 65, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};

[[maybe_unused]] static void compute_reverse_complement(char const* input, char* output,
                                                        uint64_t size) {
    for (uint64_t i = 0; i != size; ++i) {
        int c = input[i];
        output[size - i - 1] = canonicalize_basepair_reverse_map[c];
    }
}

static inline bool is_valid(int c) { return canonicalize_basepair_forward_map[c]; }

[[maybe_unused]] static bool is_valid(char const* str, uint64_t size) {
    for (uint64_t i = 0; i != size; ++i) {
        int c = str[i];
        if (canonicalize_basepair_forward_map[c] == 0) return false;
        // if (c != 'A' and c != 'C' and c != 'G' and c != 'T') return false;
    }
    return true;
}

// struct byte_range {
//     char const* begin;
//     char const* end;
// };

struct murmurhash2_64 {
    // generic range of bytes
    // static inline uint64_t hash(byte_range range, uint64_t seed) {
    //     return pthash::MurmurHash2_64(range.begin, range.end - range.begin, seed);
    // }

    // // specialization for std::string
    // static inline uint64_t hash(std::string const& val, uint64_t seed) {
    //     return MurmurHash2_64(val.data(), val.size(), seed);
    // }

    // specialization for uint64_t
    static inline uint64_t hash(uint64_t val, uint64_t seed) {
        return pthash::MurmurHash2_64(reinterpret_cast<char const*>(&val), sizeof(val), seed);
    }
};

template <typename Hasher = murmurhash2_64>
static uint64_t compute_minimizer(uint64_t kmer, uint64_t k, uint64_t m, uint64_t seed) {
    assert(m < 32);
    assert(m <= k);
    uint64_t min_hash = uint64_t(-1);
    uint64_t minimizer = uint64_t(-1);
    uint64_t mask = (uint64_t(1) << (2 * m)) - 1;
    for (uint64_t i = 0; i != k - m + 1; ++i) {
        uint64_t sub_kmer = kmer & mask;
        uint64_t hash = Hasher::hash(sub_kmer, seed);
        if (hash < min_hash) {
            min_hash = hash;
            minimizer = sub_kmer;
        }
        kmer >>= 2;
    }
    return minimizer;
}

/* not used: just for debug */
template <typename Hasher = murmurhash2_64>
static std::pair<uint64_t, uint64_t> compute_minimizer_pos(uint64_t kmer, uint64_t k, uint64_t m,
                                                           uint64_t seed) {
    assert(m < 32);
    assert(m <= k);
    uint64_t min_hash = uint64_t(-1);
    uint64_t minimizer = uint64_t(-1);
    uint64_t mask = (uint64_t(1) << (2 * m)) - 1;
    uint64_t pos = 0;
    for (uint64_t i = 0; i != k - m + 1; ++i) {
        uint64_t sub_kmer = kmer & mask;
        uint64_t hash = Hasher::hash(sub_kmer, seed);
        if (hash < min_hash) {
            min_hash = hash;
            minimizer = sub_kmer;
            pos = i;
        }
        kmer >>= 2;
    }
    return {minimizer, pos};
}

}  // namespace util

// taken from tlx
static inline std::istream& appendline(std::istream& is, std::string& str, char delim = '\n') {
    size_t size = str.size();
    size_t capacity = str.capacity();
    std::streamsize rest = capacity - size;

    if (rest == 0) {
        // if rest is zero, already expand string
        capacity = std::max(static_cast<size_t>(8), capacity * 2);
        rest = capacity - size;
    }

    // give getline access to all of capacity
    str.resize(capacity);

    // get until delim or rest is filled
    is.getline(const_cast<char*>(str.data()) + size, rest, delim);

    // gcount includes the delimiter
    size_t new_size = size + is.gcount();

    // is failbit set?
    if (!is) {
        // if string ran out of space, expand, and retry
        if (is.gcount() + 1 == rest) {
            is.clear();
            str.resize(new_size);
            str.reserve(capacity * 2);
            return appendline(is, str, delim);
        }
        // else fall through and deliver error
    } else if (!is.eof()) {
        // subtract delimiter
        --new_size;
    }

    // resize string to fit its contents
    str.resize(new_size);
    return is;
}

struct buffered_lines_iterator {
    static const uint64_t BUFFER_SIZE = 1024;

    buffered_lines_iterator(std::istream& is, uint64_t buffer_size = BUFFER_SIZE)
        : m_is(is), m_buffer_size(buffer_size), m_read_chars(0) {}

    bool fill_buffer(std::string& buffer,
                     bool force = false /* force reading of m_buffer_size characters */
    ) {
        bool empty_line_was_read = false;
        uint64_t size = buffer.size();
        uint64_t target_size = size + m_buffer_size;
        if (force) target_size += m_buffer_size;

        buffer.resize(target_size);

        char* ptr = buffer.data() + size;
        while (size != target_size) {
            // read until '\n' or rest is filled
            uint64_t rest = target_size - size;
            m_is.getline(ptr, rest, '\n');
            uint64_t read_chars = m_is.gcount();
            m_read_chars += read_chars;

            if (!m_is) {
                if (read_chars + 1 == rest) {  // '\n' not found
                    m_is.clear();
                    size += read_chars;
                    break;
                }
            } else if (!eof()) {
                assert(read_chars > 0);
                --read_chars;  // discard the delimiter
            }

            if (read_chars == 0) {  // empty line was read
                empty_line_was_read = true;
                break;
            }

            size += read_chars;
            ptr += read_chars;
        }

        buffer.resize(size);
        return empty_line_was_read;
    }

    bool eof() const { return m_is.eof(); }

    uint64_t read_chars() const { return m_read_chars; }

private:
    std::istream& m_is;
    uint64_t m_buffer_size;
    uint64_t m_read_chars;
};

}  // namespace sshash
