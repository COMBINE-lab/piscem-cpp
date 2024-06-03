
#ifndef LIBRADICL_ALIGNMENT_RECORD_HPP
#define LIBRADICL_ALIGNMENT_RECORD_HPP



#include "Type.hpp"
#include "Tags.hpp"
#include "Byte_Array.hpp"

#include <cstdint>


namespace RAD
{

class Aln_Record
{

private:

    Type::u32 ref_id;   // The id of the reference where this alignment occurs.
    Type::u8 aln_type;  // Encodes the type of alignment to follow. Valid alignment types are: 0 — read fw, 1 — read rc.
    Tag_List tag;   // Array of specified alignment-level tags for this alignment.


public:

    void set(uint32_t ref_id, uint8_t aln_type);

    template <typename T_tag_>
    void add_tag(const T_tag_& val);

    void dump(Byte_Array& byte_arr) const;
};


inline void Aln_Record::set(const uint32_t ref_id, const uint8_t aln_type)
{
    this->ref_id = ref_id;
    this->aln_type = aln_type;

    tag.clear();
}


template <typename T_tag_>
inline void Aln_Record::add_tag(const T_tag_& val)
{
    static_assert(is_RAD_type<T_tag_>());

    tag.add(val);
}


inline void Aln_Record::dump(Byte_Array& byte_arr) const
{
    byte_arr.add(ref_id);
    byte_arr.add(aln_type);
    byte_arr.add(tag);
}

}



#endif
