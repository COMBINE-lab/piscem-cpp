
#ifndef LIBRADICL_READ_RECORD_HPP
#define LIBRADICL_READ_RECORD_HPP



#include "Type.hpp"
#include "Byte_Array.hpp"
#include "Alignment_Record.hpp"
#include "Buffer.hpp"

#include <cstdint>
#include <cassert>


namespace RAD
{

class Read
{

protected:

    Type::u32 aln_count;    // Number of alignment records to follow for this read.

    Byte_Array byte_arr;    // Container for tag-values and alignment records.


public:

    template <typename T_tag_>
    void add_tag(const T_tag_& val);

    void add_aln_rec(const Aln_Record& aln_rec);

    const auto& byte_array() const { return byte_arr; }

    std::size_t size_in_bytes() const { return byte_arr.size(); }

    void dump(Buffer& buf) const;
};


class Single_End_Read : public Read
{

private:

    Type::u16 read_len; // Currently don’t support reads > 65,536bp.


public:

    void set(uint32_t aln_count, uint16_t read_len);
};


class Paired_End_Read : public Read
{

private:

    Type::u16 read_l_end_len;
    Type::u16 read_r_end_len;


public:

    void set(uint32_t aln_count, uint16_t read_l_end_len, uint16_t read_r_end_len);

};


template <typename T_tag_>
inline void Read::add_tag(const T_tag_& val)
{
    static_assert(is_RAD_type<T_tag_>());

    byte_arr.add(val);
}


inline void Read::add_aln_rec(const Aln_Record& aln_rec)
{
    aln_rec.dump(byte_arr);
}


inline void Read::dump(Buffer& buf) const
{
    buf.add(byte_arr);
}


inline void Single_End_Read::set(const uint32_t aln_count, const uint16_t read_len)
{
    this->aln_count = aln_count;
    this->read_len = read_len;

    byte_arr.clear();
    byte_arr.add(this->aln_count);
    byte_arr.add(this->read_len);
}


inline void Paired_End_Read::set(const uint32_t aln_count, const uint16_t read_l_end_len, const uint16_t read_r_end_len)
{
    this->aln_count = aln_count;
    this->read_l_end_len = read_l_end_len;
    this->read_r_end_len = read_r_end_len;

    byte_arr.clear();
    byte_arr.add(this->aln_count);
    byte_arr.add(this->read_l_end_len);
    byte_arr.add(this->read_r_end_len);
}

}



#endif
