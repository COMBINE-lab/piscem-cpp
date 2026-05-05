#ifndef LEAN_STREAMING_QUERY_HPP
#define LEAN_STREAMING_QUERY_HPP

#include "essentials.hpp"
#include "../external/sshash/include/streaming_query.hpp"
#include "../external/sshash/include/util.hpp"
#include "basic_contig_table.hpp"
#include "util_piscem.hpp"
#include <limits>
#include <cstring>

namespace piscem {

enum class KmerMatchResult : uint8_t {
    NO_MATCH = 0,
    IDENTITY_MATCH = 1,
    TWIN_MATCH = 2
};

class lean_read_iterator {
    using dict_t = piscem::piscem_dictionary;
    using kmer_t = dict_t::kmer_type;

    // sshash streaming query engine
    sshash::streaming_query<dict_t, false> m_engine;

    // contig table for locate
    sshash::basic_contig_table const* m_bct;
    uint64_t m_prev_contig_id;
    sshash::util::contig_span m_ctg_span;

    // rolling k-mer state (maintained independently of sshash engine)
    uint64_t m_fw;
    uint64_t m_rc;
    uint64_t m_k;
    uint64_t m_fw_shift;   // 2 * (k - 1), precomputed
    uint64_t m_rc_mask;    // (1 << 2k) - 1, precomputed

    // read state
    char const* m_seq;
    int32_t m_seq_len;
    int32_t m_pos;          // current k-mer start position on the read
    int32_t m_last_invalid; // position of the most recent non-ACGT character
    bool m_valid;           // current k-mer has no Ns (k chars since last invalid)
    bool m_exhausted;

    // engine sync state
    bool m_engine_synced;   // true if engine has been called for current position

    static constexpr uint64_t invalid_contig_id = std::numeric_limits<uint64_t>::max();

    static inline uint64_t char_to_2bit(char c) {
        return kmer_t::char_to_uint(c);
    }

    static inline uint64_t complement_2bit(uint64_t b) {
#ifdef SSHASH_USE_TRADITIONAL_NUCLEOTIDE_ENCODING
        // A=00, C=01, G=10, T=11 => complement = 3 - x = ~x & 3
        return (~b) & 0x3;
#else
        // A=00, C=01, T=10, G=11 => complement: XOR with 0b10
        return b ^ 0x2;
#endif
    }

    static inline bool is_valid_char(char c) {
        return kmer_t::is_valid(c);
    }

    inline void roll_forward(char new_char) {
        uint64_t b = char_to_2bit(new_char);
        m_fw = (m_fw >> 2) | (b << m_fw_shift);
        m_rc = ((m_rc << 2) | complement_2bit(b)) & m_rc_mask;
    }

    inline void build_kmer_at(char const* kmer_start) {
        m_fw = 0;
        for (uint64_t i = 0; i < m_k; ++i) {
            m_fw |= char_to_2bit(kmer_start[i]) << (2 * i);
        }
        m_rc = 0;
        for (uint64_t i = 0; i < m_k; ++i) {
            uint64_t b = char_to_2bit(kmer_start[m_k - 1 - i]);
            m_rc |= complement_2bit(b) << (2 * i);
        }
    }

    inline void refresh_contig_span(uint64_t string_id) {
        if (string_id != m_prev_contig_id) {
            auto start_pos = m_bct->m_ctg_offsets.access(string_id);
            auto end_pos = m_bct->m_ctg_offsets.access(string_id + 1);
            size_t len = end_pos - start_pos;
            m_ctg_span = {m_bct->m_ctg_entries.get_iterator_at(start_pos),
                          m_bct->m_ctg_entries.get_iterator_at(start_pos + len), len};
            m_prev_contig_id = string_id;
        }
    }

public:
    lean_read_iterator(dict_t const* d, sshash::basic_contig_table const& bct)
        : m_engine(d)
        , m_bct(&bct)
        , m_prev_contig_id(invalid_contig_id)
        , m_ctg_span()
        , m_fw(0), m_rc(0)
        , m_k(d->k())
        , m_fw_shift(2 * (d->k() - 1))
        , m_rc_mask((uint64_t(1) << (2 * d->k())) - 1)
        , m_seq(nullptr), m_seq_len(0)
        , m_pos(0), m_last_invalid(-1)
        , m_valid(false), m_exhausted(true)
        , m_engine_synced(false) {}

    lean_read_iterator(const lean_read_iterator&) = delete;
    lean_read_iterator(lean_read_iterator&&) = default;

    // Start iterating over a new read sequence.
    // Finds the first valid k-mer (no Ns), building the rolling kmer words.
    // Resets the sshash engine state for a new read.
    void start(char const* seq, int32_t len) {
        m_seq = seq;
        m_seq_len = len;
        m_pos = -1;
        m_last_invalid = -1;
        m_exhausted = false;
        m_engine_synced = false;
        m_engine.reset();
        m_prev_contig_id = invalid_contig_id;

        // scan forward to find the first valid k-mer
        find_next_valid(-1);
    }

    // Is the iterator past the end of the read?
    inline bool is_exhausted() const { return m_exhausted; }

    // Does the current k-mer contain only valid bases?
    inline bool kmer_is_valid() const { return m_valid; }

    // Current position on the read (0-based).
    inline int32_t pos() const { return m_pos; }

    // Advance to the next k-mer position. Handles N-skipping:
    // if the next character is invalid, scans forward to the next valid k-mer.
    // Returns the new position (or marks exhausted).
    inline int32_t operator++() {
        if (m_exhausted) return m_pos;
        int32_t next_j = m_pos + static_cast<int32_t>(m_k);
        if (next_j >= m_seq_len) {
            m_exhausted = true;
            return m_pos;
        }

        char c = m_seq[next_j];
        if (is_valid_char(c)) {
            roll_forward(c);
            m_pos++;
            m_engine_synced = false;
        } else {
            m_last_invalid = next_j;
            m_engine.reset();
            m_engine_synced = false;
            find_next_valid(next_j);
        }
        return m_pos;
    }

    // Advance by n positions. For small n, rolls incrementally.
    // For large n (> k), jumps directly and rebuilds from scratch (O(k) vs O(n)).
    // Returns the actual number of positions advanced.
    inline int32_t advance(int32_t n) {
        if (m_exhausted || n <= 0) return 0;
        int32_t start_pos = m_pos;

        if (n > static_cast<int32_t>(m_k)) {
            jump_to(m_pos + n);
        } else {
            for (int32_t i = 0; i < n && !m_exhausted; ++i) {
                operator++();
            }
        }
        return m_pos - start_pos;
    }

    // Jump to an arbitrary position on the read.
    // Rebuilds the k-mer from scratch and resets the engine.
    // If the target position contains an N in the k-mer window,
    // scans forward to the next valid k-mer.
    void jump_to(int32_t target_pos) {
        if (target_pos + static_cast<int32_t>(m_k) > m_seq_len) {
            m_exhausted = true;
            return;
        }
        m_engine.reset();
        m_engine_synced = false;

        // check for Ns in the k-mer window at target_pos
        m_last_invalid = -1;
        for (int32_t j = target_pos; j < target_pos + static_cast<int32_t>(m_k); ++j) {
            if (!is_valid_char(m_seq[j])) {
                m_last_invalid = j;
            }
        }

        if (m_last_invalid == -1) {
            // all chars valid — build the k-mer directly
            m_pos = target_pos;
            m_valid = true;
            build_kmer_at(m_seq + target_pos);
        } else {
            // has an N — scan forward from the last invalid position
            find_next_valid(m_last_invalid);
        }
    }

    // Compare the current rolling k-mer against a reference k-mer (uint64).
    // Does NOT require the sshash engine to be synced.
    inline KmerMatchResult is_equivalent(uint64_t ref_kmer) const {
        if (ref_kmer == m_fw) return KmerMatchResult::IDENTITY_MATCH;
        if (ref_kmer == m_rc) return KmerMatchResult::TWIN_MATCH;
        return KmerMatchResult::NO_MATCH;
    }

    // Access the raw forward and reverse-complement k-mer words.
    inline uint64_t fw_word() const { return m_fw; }
    inline uint64_t rc_word() const { return m_rc; }

    // Perform a full lookup through the sshash engine at the current position.
    // This syncs the engine state and returns the full lookup_result.
    // If the current k-mer is not valid, returns an empty result.
    inline sshash::lookup_result lookup() {
        if (!m_valid) {
            m_engine.reset();
            m_engine_synced = true;
            return sshash::lookup_result();
        }

        if (!m_engine_synced) {
            // The engine may be behind or desynchronized due to skips.
            // We must reset and feed the current k-mer as a fresh start.
            m_engine.reset();
            auto res = m_engine.lookup(m_seq + m_pos);
            m_engine_synced = true;

            if (res.kmer_id != sshash::constants::invalid_uint64) {
                refresh_contig_span(res.string_id);
            }
            return res;
        }

        // Engine is already synced — call lookup for the next streaming position.
        auto res = m_engine.lookup(m_seq + m_pos);

        if (res.kmer_id != sshash::constants::invalid_uint64) {
            refresh_contig_span(res.string_id);
        }
        return res;
    }

    // Perform a streaming lookup: advance one position and lookup.
    // This is the hot-path for sequential scanning.
    // The engine stays synced after this call, so extension logic applies.
    inline sshash::lookup_result next_lookup() {
        operator++();
        if (m_exhausted) return sshash::lookup_result();
        // After operator++, the engine is desynced. But we want streaming behavior
        // where sshash's extension logic fires. So we call lookup on the char pointer
        // — sshash will attempt extension from its internal state.
        auto res = m_engine.lookup(m_seq + m_pos);
        m_engine_synced = true;

        if (res.kmer_id != sshash::constants::invalid_uint64) {
            refresh_contig_span(res.string_id);
        }
        return res;
    }

    // Lookup at the current position in streaming mode.
    // Assumes we arrived here via operator++ and the engine is one step behind.
    inline sshash::lookup_result streaming_lookup() {
        if (m_exhausted || !m_valid) return sshash::lookup_result();
        auto res = m_engine.lookup(m_seq + m_pos);
        m_engine_synced = true;

        if (res.kmer_id != sshash::constants::invalid_uint64) {
            refresh_contig_span(res.string_id);
        }
        return res;
    }

    inline sshash::util::contig_span contig_span() const { return m_ctg_span; }
    inline bool is_present() const {
        return m_engine_synced &&
               m_engine.result().kmer_id != sshash::constants::invalid_uint64;
    }

    uint64_t num_searches() const { return m_engine.num_searches(); }
    uint64_t num_extensions() const { return m_engine.num_extensions(); }
    uint64_t k() const { return m_k; }

private:
    // Scan forward from position last_bad to find the next valid k-mer.
    // last_bad is the index of the most recent invalid character.
    void find_next_valid(int32_t last_bad) {
        int32_t start = last_bad + 1;
        // We need k consecutive valid chars starting from 'start'.
        // Scan characters one by one.
        int32_t j = start;
        int32_t valid_run = 0;

        // If we had a partial valid run before last_bad, we start fresh.
        m_fw = 0;
        m_rc = 0;

        while (j < m_seq_len) {
            char c = m_seq[j];
            if (is_valid_char(c)) {
                // shift the forward kmer
                uint64_t b = char_to_2bit(c);
                m_fw = (m_fw >> 2) | (b << m_fw_shift);
                m_rc = ((m_rc << 2) | complement_2bit(b)) & m_rc_mask;
                valid_run++;

                if (valid_run >= static_cast<int32_t>(m_k)) {
                    m_pos = j - static_cast<int32_t>(m_k) + 1;
                    m_valid = true;
                    m_engine_synced = false;
                    return;
                }
            } else {
                m_last_invalid = j;
                valid_run = 0;
                m_fw = 0;
                m_rc = 0;
            }
            ++j;
        }

        // no valid k-mer found
        m_exhausted = true;
        m_valid = false;
    }
};

}  // namespace piscem

#endif
