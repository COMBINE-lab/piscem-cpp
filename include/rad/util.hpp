#ifndef __RAD_UTIL_HPP__
#define __RAD_UTIL_HPP__

#include <fstream>

#include "../reference_index.hpp"
#include "rad_header.hpp"
#include "rad_writer.hpp"

namespace rad {
namespace util {

size_t write_rad_header(mindex::reference_index& ri, std::ofstream& rad_file) {
    constexpr size_t bc_length = 16;
    constexpr size_t umi_length = 12;

    rad_writer bw;
    //  RADHeader
    rad_header rh;
    // right now, all formats are effectively single-end
    rh.is_paired(false);
    for (size_t i = 0; i < ri.num_refs(); ++i) { rh.add_refname(ri.ref_name(i)); }
    rh.dump_to_bin(bw);

    // where we will write the number of chunks when we know
    // how many there are
    size_t chunk_offset = bw.num_bytes() - sizeof(uint64_t);

    // ### start of tags

    // Tags we will have
    // write the tag meta-information section

    // File-level tag description
    uint16_t file_level_tags{2};
    bw << file_level_tags;

    // cblen
    uint8_t type_id{2};
    bw << std::string("cblen");
    bw << type_id;

    bw << std::string("ulen");
    bw << type_id;

    // read-level tag description
    uint16_t read_level_tags{2};
    bw << read_level_tags;

    // barcode
    bw << std::string("b");
    if (bc_length > 32) {
        type_id = 8;
    } else if (bc_length > 16) {
        // 17 - 32 bases
        type_id = 4;
    } else {
        // <= 16 bases
        type_id = 3;
    }
    bw << type_id;

    // umi
    bw << std::string("u");
    if (umi_length > 32) {
        type_id = 8;
    } else if (umi_length > 16) {
        // 17 - 32 bases
        type_id = 4;
    } else {
        // <= 16 bases
        type_id = 3;
    }
    bw << type_id;

    // alignment-level tag description
    uint16_t aln_level_tags{1};
    bw << aln_level_tags;
    // we maintain orientation
    // bw << std::string("orientation");
    // type_id = 1;
    // bw << type_id;

    // and reference id
    bw << std::string("compressed_ori_refid");
    type_id = 3;
    bw << type_id;

    // ### end of tag definitions

    // the actual file-level tags
    bw << static_cast<uint16_t>(16);
    bw << static_cast<uint16_t>(12);

    rad_file << bw;
    bw.clear();
    return chunk_offset;
}

}  // namespace util
}  // namespace rad

#endif  //__RAD_UTIL_HPP__
