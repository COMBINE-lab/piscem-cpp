
#ifndef LIBRADICL_BYTE_ARRAY_HPP
#define LIBRADICL_BYTE_ARRAY_HPP



#include "Type.hpp"
#include "Tags.hpp"

#include <cstdint>
#include <cstddef>
#include <cassert>


namespace RAD
{

class Byte_Array
{

private:

    static constexpr std::size_t cap_default = 4 * 1024;    // 4KB.
    std::size_t cap;

    uint8_t* arr;
    std::size_t sz;


    template <typename T_> void add_POD(T_ val);

    void expand(std::size_t req_cap);


public:

    explicit Byte_Array(std::size_t cap = cap_default);

    Byte_Array(const Byte_Array&) = delete;
    Byte_Array(Byte_Array&&) = delete;
    Byte_Array& operator=(const Byte_Array&) = delete;
    Byte_Array& operator=(Byte_Array&&) = delete;

    ~Byte_Array();

    auto empty() const { return sz == 0; }

    auto size() const { return sz; }

    const auto* data() const { return arr; }

    void clear() { sz = 0; }

    template <typename T_> void add(const T_& val);
};


template <typename T_>
inline void Byte_Array::add(const T_& val)
{
    static_assert(is_RAD_type<T_>());

    add_POD(val.val());
}


template <>
inline void Byte_Array::add<Type::null>(const Type::null&)
{}


template <>
inline void Byte_Array::add<Type::str>(const Type::str& s)
{
    add(Type::u16(s.val().length()));
    for(const auto v : s.val())
        add_POD(v);
}


template <>
inline void Byte_Array::add<Tag_List>(const Tag_List& tag)
{
    tag.write(*this);
}


template <typename T_>
inline void Byte_Array::add_POD(const T_ val)
{
    static_assert(std::is_pod<T_>());

    if(sz + sizeof(val) > cap)
        expand(sz + sizeof(val));

    assert(sz + sizeof(val) <= cap);
    std::memcpy(reinterpret_cast<char*>(arr + sz), reinterpret_cast<const char*>(&val), sizeof(val));
    sz += sizeof(val);
}


inline void Byte_Array::expand(const std::size_t req_cap)
{
    while(cap < req_cap)
        cap *= 2;

    arr = static_cast<uint8_t*>(std::realloc(arr, req_cap));
}

}



#endif
