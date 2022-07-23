#ifndef __RAD_UTIL_HPP__
#define __RAD_UTIL_HPP__

#include <fstream>

#include "../include/parallel_hashmap/phmap.h"
#include "../include/mapping/utils.hpp"
#include "../Kmer.hpp"
#include "../reference_index.hpp"
#include "rad_header.hpp"
#include "rad_writer.hpp"

namespace rad {
namespace util {

using umi_kmer_t = combinelib::kmers::Kmer<31, 2>;
using bc_kmer_t = combinelib::kmers::Kmer<31, 3>;

size_t write_rad_header(mindex::reference_index& ri, size_t bc_length, size_t umi_length,
                        std::ofstream& rad_file) {
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
    bw << static_cast<uint16_t>(bc_length);
    bw << static_cast<uint16_t>(umi_length);

    rad_file << bw;
    bw.clear();
    return chunk_offset;
}

inline void write_to_rad_stream(bc_kmer_t& bck, umi_kmer_t& umi,
                                mapping::util::MappingType map_type,
                                std::vector<mapping::util::simple_hit>& accepted_hits,
                                phmap::flat_hash_map<uint64_t, uint32_t>& unmapped_bc_map,
                                uint32_t& num_reads_in_chunk, rad_writer& bw) {
    const uint32_t barcode_len = bc_kmer_t::k();
    const uint32_t umi_len = umi_kmer_t::k();

    if (map_type == mapping::util::MappingType::SINGLE_MAPPED) {
        // number of mappings
        bw << static_cast<uint32_t>(accepted_hits.size());

        // bc
        // if we can fit the barcode into an integer
        if (barcode_len <= 32) {
            if (barcode_len <= 16) {  // can use 32-bit int
                uint32_t shortbck = static_cast<uint32_t>(0x00000000FFFFFFFF & bck.word(0));
                bw << shortbck;
            } else {  // must use 64-bit int
                bw << bck.word(0);
            }
        } else {  // must use a string for the barcode
            std::cerr << "should not happen\n";
        }

        // umi
        if (umi_len <= 16) {  // if we can use 32-bit int
            uint64_t umiint = umi.word(0);
            uint32_t shortumi = static_cast<uint32_t>(0x00000000FFFFFFFF & umiint);
            bw << shortumi;
        } else if (umi_len <= 32) {  // if we can use 64-bit int
            uint64_t umiint = umi.word(0);
            bw << umiint;
        } else {  // must use string
            std::cerr << "should not happen\n";
        }

        for (auto& aln : accepted_hits) {
            uint32_t fw_mask = aln.is_fw ? 0x80000000 : 0x00000000;
            bw << (aln.tid | fw_mask);
        }
        ++num_reads_in_chunk;
    } else {
        unmapped_bc_map[bck.word(0)] += 1;
    }
}

}  // namespace util
}  // namespace rad

#endif  //__RAD_UTIL_HPP__
