#pragma once
#include "../dictionary.hpp"
#include "../minimizer_enumerator.hpp"
#include "../util.hpp"
#include "../Kmer.hpp"

namespace sshash {

struct contig_info_query_canonical_parsing {
    contig_info_query_canonical_parsing(dictionary const* dict)

        : num_searches(0)
        , num_extensions(0)

        , m_dict(dict)

        , m_minimizer_enum(dict->m_k, dict->m_m, dict->m_seed)
        , m_minimizer_enum_rc(dict->m_k, dict->m_m, dict->m_seed)
        , m_minimizer_not_found(false)
        , m_start(true)
        , m_curr_minimizer(constants::invalid)
        , m_prev_minimizer(constants::invalid)
        , m_kmer(constants::invalid)

        , m_shift(2 * (dict->m_k - 1))
        , m_k(dict->m_k)
        , m_m(dict->m_m)
        , m_seed(dict->m_seed)

        , m_string_iterator(dict->m_buckets.strings, 0)
        , m_begin(0)
        , m_end(0)
        , m_pos_in_window(0)
        , m_window_size(0)
        , m_prev_query_offset(0)
        , m_reverse(false)
        , m_prev_contig_id(0)
        , m_prev_contig_length(0)
        , m_prev_contig_offset(0)
        , m_prev_global_pos(0) {
        assert(m_dict->m_canonical_parsing);
    }

    inline void start() { m_start = true; }

    struct query_result {
        bool is_valid;
        bool is_member;
        uint64_t contig_id;
        uint64_t contig_length;
        uint64_t contig_offset;
        uint64_t global_pos;
        bool is_forward;
    };

    // Given a k-mer `k-mer`, return the query_result associated with it.
    // NOTE: think about a more general "caching" mechanism — what if the expected
    query_result get_contig_pos(const uint64_t km, const uint64_t km_rc,
                                const uint64_t query_offset) {
        m_kmer = km;
        m_kmer_rc = km_rc;

        // if the current query offset position is
        // the next position after the stored query
        // offset position, then we can apply the
        // relevant optimizations.  Otherwise, we
        // should consider this as basically a "new"
        // query
        if (!m_start) { m_start = (m_prev_query_offset + 1) != query_offset; }
        m_prev_query_offset = query_offset;

        return _get_contig_pos();
    }

    // Given a k-mer `k-mer`, return the query_result associated with it.
    // NOTE: think about a more general "caching" mechanism — what if the expected
    query_result get_contig_pos(const char* kmer) {
        constexpr uint64_t invalid = std::numeric_limits<uint64_t>::max();
        query_result qr = {false, false, invalid, invalid, invalid, invalid, true};

        /* validation */
        bool is_valid = m_start ? util::is_valid(kmer, m_k) : util::is_valid(kmer[m_k - 1]);
        if (!is_valid) {
            m_start = true;
            return qr;
        }
        /*************/

        /* compute kmer and minimizer */
        if (!m_start) {
            m_kmer >>= 2;
            m_kmer += (util::char_to_uint64(kmer[m_k - 1])) << m_shift;
            assert(m_kmer == util::string_to_uint64_no_reverse(kmer, m_k));
        } else {
            m_kmer = util::string_to_uint64_no_reverse(kmer, m_k);
        }
        m_kmer_rc = util::compute_reverse_complement(m_kmer, m_k);
        return _get_contig_pos();
    }

    // assumes that m_kmer and m_kmer_rc have been set, and gets the result
    query_result _get_contig_pos() {
        constexpr uint64_t invalid = std::numeric_limits<uint64_t>::max();
        query_result qr = {false, false, invalid, invalid, invalid, invalid, true};

        m_curr_minimizer = m_minimizer_enum.next(m_kmer, m_start);
        assert(m_curr_minimizer == util::compute_minimizer(m_kmer, m_k, m_m, m_seed));
        constexpr bool reverse = true;
        uint64_t minimizer_rc = m_minimizer_enum_rc.next<reverse>(m_kmer_rc, m_start);
        assert(minimizer_rc == util::compute_minimizer(m_kmer_rc, m_k, m_m, m_seed));
        m_curr_minimizer = std::min<uint64_t>(m_curr_minimizer, minimizer_rc);
        /******************************/

        /* no optimizations */
        /*
        bool answer = false;
        m_minimizer_not_found = false;
        locate_bucket();
        int ret = is_member(qr);
        if (ret == return_value::MINIMIZER_NOT_FOUND) {
            m_minimizer_not_found = true;
            answer = false;
        } else {
            answer = (ret == return_value::KMER_FOUND);
        }
        */
        /* compute answer */
        bool answer = false;
        if (same_minimizer()) {
            if (m_minimizer_not_found) {
                answer = false;
            } else if (extends()) {
                extend(qr);
                answer = true;
            } else {
                int ret = is_member(qr);
                answer = (ret == return_value::KMER_FOUND);
            }
        } else {
            m_minimizer_not_found = false;
            locate_bucket();

            // Try to extend matching even when we change minimizer.
            if (extends()) {
                extend(qr);
                answer = true;
            } else {
                int ret = is_member(qr);
                if (ret == return_value::MINIMIZER_NOT_FOUND) {
                    m_minimizer_not_found = true;
                    answer = false;
                } else {
                    answer = (ret == return_value::KMER_FOUND);
                }
            }
        }
        /******************/

        /* update state */
        m_prev_minimizer = m_curr_minimizer;
        m_start = false;
        /****************/

        // assert(m_dict->is_member(kmer) == answer);

        qr.is_member = answer;
        // record the query information we will
        // need if the next query is answered via
        // extension.
        if (answer) {
            m_prev_contig_id = qr.contig_id;
            m_prev_contig_length = qr.contig_length;
            m_prev_contig_offset = qr.contig_offset;
            m_prev_global_pos = qr.global_pos;
        }
        return qr;  //{true, answer};
    }

    /* counts */
    uint64_t num_searches, num_extensions;

private:
    dictionary const* m_dict;

    /* (kmer,minimizer) state */
    minimizer_enumerator<> m_minimizer_enum;
    minimizer_enumerator<> m_minimizer_enum_rc;
    bool m_minimizer_not_found;
    bool m_start;
    uint64_t m_curr_minimizer, m_prev_minimizer;
    uint64_t m_kmer, m_kmer_rc;

    /* constants */
    uint64_t m_shift, m_k, m_m, m_seed;

    /* string state */
    bit_vector_iterator m_string_iterator;
    uint64_t m_begin, m_end;
    uint64_t m_pos_in_window, m_window_size;
    uint64_t m_prev_query_offset;
    bool m_reverse;
    uint64_t m_prev_contig_id, m_prev_contig_length, m_prev_contig_offset, m_prev_global_pos;

    enum return_value { MINIMIZER_NOT_FOUND = 0, KMER_FOUND = 1, KMER_NOT_FOUND = 2 };

    inline bool same_minimizer() const { return m_curr_minimizer == m_prev_minimizer; }

    void locate_bucket() {
        uint64_t bucket_id = (m_dict->m_minimizers).lookup(m_curr_minimizer);
        std::tie(m_begin, m_end) = (m_dict->m_buckets).locate_bucket(bucket_id);
    }

    int is_member(query_result& qr) {
        bool check_minimizer = !same_minimizer();
        if (!m_dict->m_skew_index.empty()) {
            uint64_t num_strings_in_bucket = m_end - m_begin;
            uint64_t log2_num_strings_in_bucket = util::ceil_log2_uint32(num_strings_in_bucket);

            // If the current minimizer occurs frequently enough to be held in the
            // skew index.
            if (log2_num_strings_in_bucket > (m_dict->m_skew_index).min_log2) {
                // try to lookup in the forward orientation first
                uint64_t p = m_dict->m_skew_index.lookup(m_kmer, log2_num_strings_in_bucket);
                if (p < num_strings_in_bucket) {
                    int ret = is_member(m_begin + p, m_begin + p + 1, check_minimizer, qr);
                    if (ret != return_value::KMER_NOT_FOUND) return ret;
                    check_minimizer = false;
                }
                // then in the reverse complement orientation
                uint64_t p_rc = m_dict->m_skew_index.lookup(m_kmer_rc, log2_num_strings_in_bucket);
                if (p_rc < num_strings_in_bucket) {
                    int ret = is_member(m_begin + p_rc, m_begin + p_rc + 1, check_minimizer, qr);
                    if (ret != return_value::KMER_NOT_FOUND) return ret;
                }
                return return_value::KMER_NOT_FOUND;
            }
        }
        return is_member(m_begin, m_end, check_minimizer, qr);
    }

    int is_member(uint64_t begin, uint64_t end, bool check_minimizer, query_result& qr) {
        // check every string (super k-mer occurrence) to which this
        // k-mer might belong
        for (uint64_t string_id = begin; string_id != end; ++string_id) {
            uint64_t offset = (m_dict->m_buckets).offsets.access(string_id);
            uint64_t pos_in_string = 2 * offset;
            m_reverse = false;
            m_string_iterator.at(pos_in_string);
            auto [kmer_id, contig_id, offset_end, prev_end] =
                (m_dict->m_buckets).offset_to_contig_info(offset, m_k);
            // auto [kmer_id, offset_end] = (m_dict->m_buckets).offset_to_id(offset, m_k);
            (void)prev_end;
            (void)kmer_id;
            m_pos_in_window = 0;
            m_window_size = std::min<uint64_t>(m_k - m_m + 1, offset_end - offset - m_k + 1);

            while (m_pos_in_window != m_window_size) {
                uint64_t val = m_string_iterator.read(2 * m_k);

                if (check_minimizer and string_id == begin and m_pos_in_window == 0) {
                    uint64_t val_rc = util::compute_reverse_complement(val, m_k);
                    uint64_t minimizer =
                        std::min<uint64_t>(util::compute_minimizer(val, m_k, m_m, m_seed),
                                           util::compute_minimizer(val_rc, m_k, m_m, m_seed));
                    if (minimizer != m_curr_minimizer) return return_value::MINIMIZER_NOT_FOUND;
                }

                m_string_iterator.eat(2);
                m_pos_in_window += 1;
                pos_in_string += 2;
                assert(m_pos_in_window <= m_window_size);

                if (m_kmer == val) {
                    num_searches += 1;
                    qr.is_forward = true;
                    uint64_t nuc_pos = pos_in_string / 2;
                    qr.contig_offset = (nuc_pos - (prev_end + 1));
                    qr.contig_id = contig_id;
                    qr.contig_length = (offset_end - prev_end);
                    qr.global_pos = nuc_pos;
                    return return_value::KMER_FOUND;
                }

                if (m_kmer_rc == val) {
                    m_reverse = true;
                    pos_in_string -= 2;
                    num_searches += 1;
                    m_string_iterator.at(pos_in_string + 2 * (m_k - 1));
                    qr.is_forward = false;
                    uint64_t nuc_pos = pos_in_string / 2;
                    qr.contig_offset = (nuc_pos - prev_end);
                    qr.contig_id = contig_id;
                    qr.contig_length = (offset_end - prev_end);
                    qr.global_pos = nuc_pos;
                    return return_value::KMER_FOUND;
                }
            }
        }

        return return_value::KMER_NOT_FOUND;
    }

    // If `extends()` returned true, then we do the actual
    // extension here.  
    // Move the relevant offsets and iterators depending on 
    // the direction of the hit, and set the relevant state for the 
    // query result.
    inline void extend(query_result& qr) {
        if (m_reverse) {
            m_string_iterator.eat_reverse(2);
            m_pos_in_window -= 1;
            m_prev_contig_offset -= 1;
            m_prev_global_pos -= 1;
            qr.is_forward = false;
            assert(m_pos_in_window >= 1);
        } else {
            m_string_iterator.eat(2);
            m_pos_in_window += 1;
            m_prev_contig_offset += 1;
            m_prev_global_pos += 1;
            qr.is_forward = true;
            assert(m_pos_in_window <= m_window_size);
        }
        qr.contig_id = m_prev_contig_id;
        qr.contig_length = m_prev_contig_length;
        qr.contig_offset = m_prev_contig_offset;
        qr.global_pos = m_prev_global_pos;
    }

    inline bool extends() {
        if (m_reverse) {
            if (m_pos_in_window == 1) return false;
            if (m_kmer_rc == m_string_iterator.read_reverse(2 * m_k)) {
                ++num_extensions;
                return true;
            }
            return false;
        }
        if (m_pos_in_window == m_window_size) return false;
        if (m_kmer == m_string_iterator.read(2 * m_k)) {
            ++num_extensions;
            return true;
        }
        return false;
    }
};

}  // namespace sshash
