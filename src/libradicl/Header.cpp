
#include "../../include/libradicl/Header.hpp"

#include <cassert>


namespace RAD
{

Header::Header(const bool is_paired, const uint64_t ref_count, const std::vector<std::string>& refs, const uint64_t num_chunks):
      is_paired(is_paired)
    , ref_count(ref_count)
    , num_chunks(num_chunks)
{
    assert(refs.size() == ref_count);

    ref_names.reserve(ref_count);
    for(const auto& v : refs)
        ref_names.emplace_back(v);
}


void Header::write(Buffer& buf) const
{
    buf.add(is_paired);
    buf.add(ref_count);
    for(const auto& v : ref_names)
        buf.add(v);
    buf.add(num_chunks);
}

}
