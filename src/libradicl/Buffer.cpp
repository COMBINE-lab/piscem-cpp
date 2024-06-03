
#include "../../include/libradicl/Buffer.hpp"

#include <cstdlib>


namespace RAD
{

Buffer::Buffer(const std::size_t cap, std::ofstream& os):
      cap(cap)
    , buf(static_cast<uint8_t*>(std::malloc(cap)))
    , sz(0)
    , os(os)
{}


Buffer::~Buffer()
{
    std::free(buf);
}

}
