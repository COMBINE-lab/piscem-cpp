
#ifndef LIBRADICL_TYPE_HPP
#define LIBRADICL_TYPE_HPP



#include "boost/variant/variant.hpp"

#include <string>
#include <cstdint>
#include <type_traits>


namespace RAD
{

namespace Type
{

// Allowed types in a RAD file.
template <typename T_>
class RAD_Type
{
private:

    T_ val_;

public:

    RAD_Type() {}

    RAD_Type(const T_ val): val_(val) {}

    auto val() const { return val_; }

    auto operator()() const { return val_; }

    static constexpr uint8_t type_id()
    {
        using std::is_same;

        if constexpr(is_same<T_, void>())           return 255;
        if constexpr(is_same<T_, bool>())           return 0;
        if constexpr(is_same<T_, uint8_t>())        return 1;
        if constexpr(is_same<T_, uint16_t>())       return 2;
        if constexpr(is_same<T_, uint32_t>())       return 3;
        if constexpr(is_same<T_, uint64_t>())       return 4;
        if constexpr(is_same<T_, float>())          return 5;
        if constexpr(is_same<T_, double>())         return 6;
        // TODO: no support for arbitrary arrays as of now.
        if constexpr(is_same<T_, std::string>())    return 8;

        return 255;
    }
};

template <>
class RAD_Type<void>
{
public:

    static constexpr uint8_t type_id() { return 255; }
};


typedef RAD_Type<void> null;
typedef RAD_Type<bool> b;
typedef RAD_Type<uint8_t> u8;
typedef RAD_Type<uint16_t> u16;
typedef RAD_Type<uint32_t> u32;
typedef RAD_Type<uint64_t> u64;
typedef RAD_Type<float> f32;
typedef RAD_Type<double> f64;
// TODO: no support for arbitrary arrays as of now.
typedef RAD_Type<std::string> str;


typedef boost::variant<Type::null, Type::b, Type::u8, Type::u16, Type::u32, Type::u64, Type::f32, Type::f64, Type::str> variant_t;

}


template <typename T_>
inline static constexpr auto is_RAD_type()
{
    using std::is_same;
    using namespace Type;
    return  is_same<T_, null>() || is_same<T_, b>() ||
            is_same<T_, u8>() || is_same<T_, u16>() || is_same<T_, u32>() || is_same<T_, u64>() ||
            is_same<T_, f32>() || is_same<T_, f64>() || is_same<T_, str>();
}

}



#endif
