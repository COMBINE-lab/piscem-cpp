
#include "../../include/libradicl/Byte_Array.hpp"

#include <cstdlib>
#include <cassert>


namespace RAD
{

Byte_Array::Byte_Array(const std::size_t cap):
      cap(cap)
    , arr(static_cast<uint8_t*>(std::malloc(cap)))
    , sz(0)
{
    assert(cap > 0);
}


Byte_Array::~Byte_Array()
{
    std::free(arr);
}

}
