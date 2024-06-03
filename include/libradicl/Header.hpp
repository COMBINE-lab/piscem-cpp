
#ifndef LIBRADICL_HEADER_HPP
#define LIBRADICL_HEADER_HPP



#include "Type.hpp"
#include "Buffer.hpp"

#include <cstdint>
#include <vector>
#include <string>


namespace RAD
{

class Header
{
private:

    const Type::b is_paired;    // Whether the mappings to follow for paired-end or single-end fragments.
    const Type::u64 ref_count;  // Number of reference sequences being aligned against.
    std::vector<Type::str> ref_names;   // Strings encoding the names of the target references.
    Type::u64 num_chunks;   // Number of chunks in the file following the header. If this value is 0, then the number of chunks is not known.

public:

    explicit Header(bool is_paired, uint64_t ref_count, const std::vector<std::string>& refs, uint64_t num_chunks = 0);

    void write(Buffer& buf) const;
};

}



#endif
