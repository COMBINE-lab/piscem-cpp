#pragma once

#include "../external/pthash/include/encoders/bit_vector.hpp"

namespace sshash {

struct bit_vector_iterator {
    bit_vector_iterator() : m_bv(nullptr) {}

    bit_vector_iterator(pthash::bit_vector const& bv, uint64_t pos) : m_bv(&bv) { at(pos); }

    void at(uint64_t pos) {
        m_pos = pos;
        m_buf = 0;
        m_avail = 0;
    }

    inline uint64_t read(uint64_t l) {
        assert(l <= 64);
        if (m_avail < l) fill_buf();
        uint64_t val = 0;
        if (l != 64) {
            val = m_buf & ((uint64_t(1) << l) - 1);
        } else {
            val = m_buf;
        }
        return val;
    }

    inline uint64_t read_reverse(uint64_t l) {
        assert(l <= 64);
        if (m_avail < l) fill_buf_reverse();
        uint64_t val = 0;
        if (l != 64) {
            val = m_buf >> (64 - l);
        } else {
            val = m_buf;
        }
        return val;
    }

    inline void eat(uint64_t l) {
        assert(l <= 64);
        if (m_avail < l) fill_buf();
        if (l != 64) m_buf >>= l;
        m_avail -= l;
        m_pos += l;
    }

    inline void eat_reverse(uint64_t l) {
        assert(l <= 64);
        if (m_avail < l) fill_buf_reverse();
        if (l != 64) m_buf <<= l;
        m_avail -= l;
        m_pos -= l;
    }

    inline uint64_t read_and_advance_by_two(uint64_t l) {
        assert(l <= 64);
        if (m_avail < l) fill_buf();
        uint64_t val = 0;
        if (l != 64) {
            val = m_buf & ((uint64_t(1) << l) - 1);
            m_buf >>= 2;
        } else {
            val = m_buf;
        }
        m_avail -= 2;
        m_pos += 2;
        return val;
    }

    inline uint64_t read_and_advance_by_two_reverse(uint64_t l) {
        assert(l <= 64);
        if (m_avail < l) fill_buf_reverse();
        uint64_t val = 0;
        if (l != 64) {
            val = m_buf >> (64 - l);
            m_buf <<= 2;
        } else {
            val = m_buf;
        }
        m_avail -= 2;
        m_pos -= 2;
        return val;
    }

    inline uint64_t get_next_two_bits() {
        if (m_avail < 2) fill_buf();
        uint64_t val = m_buf & 3;
        m_buf >>= 2;
        m_avail -= 2;
        m_pos += 2;
        return val;
    }

    inline uint64_t take(uint64_t l) {
        assert(l <= 64);
        if (m_avail < l) fill_buf();
        uint64_t val;
        if (l != 64) {
            val = m_buf & ((uint64_t(1) << l) - 1);
            m_buf >>= l;
        } else {
            val = m_buf;
        }
        m_avail -= l;
        m_pos += l;
        return val;
    }

    inline uint64_t position() const { return m_pos; }

private:
    inline void fill_buf() {
        m_buf = m_bv->get_word64(m_pos);
        m_avail = 64;
    }

    inline void fill_buf_reverse() {
        if (m_pos < 64) {
            m_buf = m_bv->get_word64(0);
            m_avail = m_pos;
            m_buf <<= (64 - m_pos);
            return;
        }
        m_buf = m_bv->get_word64(m_pos - 64);
        m_avail = 64;
    }

    pthash::bit_vector const* m_bv;
    uint64_t m_pos;
    uint64_t m_buf;
    uint64_t m_avail;
};

}  // namespace sshash