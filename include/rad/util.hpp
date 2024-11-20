#ifndef __RAD_UTIL_HPP__
#define __RAD_UTIL_HPP__

#include <fstream>

#include "../parallel_hashmap/phmap.h"
#include "../mapping/utils.hpp"
#include "../mapping/utils_bin.hpp"
#include "../Kmer.hpp"
#include "../reference_index.hpp"
#include "rad_header.hpp"
#include "rad_writer.hpp"

#include "../../external/libradicl/include/RAD_Writer.hpp"
#include "../../external/libradicl/include/Alignment_Record.hpp"
#include "../../external/libradicl/include/Read_Record.hpp"
#include "../../external/libradicl/include/Byte_Array.hpp"
#include "../../external/libradicl/include/Tags.hpp"
#include "../../external/libradicl/include/Tags.hpp"

namespace rad {
namespace util {

using umi_kmer_t = combinelib::kmers::Kmer<31, 2>;
using bc_kmer_t = combinelib::kmers::Kmer<31, 3>;

inline size_t write_rad_header(mindex::reference_index &ri, size_t bc_length,
                               size_t umi_length, std::ofstream &rad_file) {
  rad_writer bw;
  //  RADHeader
  rad_header rh;
  // right now, all formats are effectively single-end
  rh.is_paired(false);
  for (size_t i = 0; i < ri.num_refs(); ++i) {
    rh.add_refname(ri.ref_name(i));
  }
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

inline size_t write_rad_header_bulk(mindex::reference_index &ri, bool is_paired,
                                    std::ofstream &rad_file) {
  rad_writer bw;
  //  RADHeader
  rad_header rh;
  rh.is_paired(is_paired);
  for (size_t i = 0; i < ri.num_refs(); ++i) {
    rh.add_refname(ri.ref_name(i));
  }
  rh.dump_to_bin(bw);

  // where we will write the number of chunks when we know
  // how many there are
  size_t chunk_offset = bw.num_bytes() - sizeof(uint64_t);

  // ### start of tags

  // Tags we will have
  // write the tag meta-information section

  // File-level tag description
  // none right now
  uint16_t file_level_tags{1};
  bw << file_level_tags;

  bw << std::string("ref_lengths");
  uint8_t type_id{7}; // type is array
  bw << type_id;
  // length type is u32
  type_id = 3;
  bw << type_id;
  // element type is u32
  bw << type_id;

  // read-level tag description
  // will hold the type of mappings for the read
  uint16_t read_level_tags{1};
  bw << read_level_tags;

  type_id = 1;
  // fragment mapping type
  // unmapped, mapped single, orphan_left, orphan_right, paired_end
  bw << std::string("frag_map_type");
  bw << type_id;

  // alignment-level tag description
  uint16_t aln_level_tags{3};
  bw << aln_level_tags;

  // read_pair_ori, refid
  bw << std::string("compressed_ori_ref");
  type_id = 3;
  bw << type_id;

  // fragment start position
  bw << std::string("pos");
  type_id = 3;
  bw << type_id;

  // fragment length
  bw << std::string("frag_len");
  type_id = 2;
  bw << type_id;

  // ### end of tag definitions

  // the actual file-level tag
  // we've already recorded the description
  // so here, we give the length and the the
  // elements.
  uint32_t num_refs = static_cast<uint32_t>(ri.num_refs());
  // array length type u32, number of elements is num refs, element type u32
  bw << num_refs;
  for (size_t i = 0; i < num_refs; ++i) {
    uint32_t rl = static_cast<uint32_t>(ri.ref_len(i));
    bw << rl;
  }

  rad_file << bw;
  bw.clear();
  return chunk_offset;
}

inline void write_rad_header_atac(mindex::reference_index& ri, std::vector<std::string>& refs, RAD::Tag_Defn& tag_defn) {
    
    for (size_t i = 0; i < ri.num_refs(); ++i) { refs.emplace_back(ri.ref_name(i)); }

    tag_defn.add_file_tag<RAD::Type::u16>("cblen");
    tag_defn.add_file_tag<RAD::Type::v_u64>("ref_lengths");

    tag_defn.add_read_tag<RAD::Type::u32>("barcode");

    tag_defn.add_aln_tag<RAD::Type::u32>("ref");
    tag_defn.add_aln_tag<RAD::Type::u8>("type");
    tag_defn.add_aln_tag<RAD::Type::u32>("start_pos");
    tag_defn.add_aln_tag<RAD::Type::u16>("frag_len");
}

inline void write_to_rad_stream(bc_kmer_t &bck, umi_kmer_t &umi,
                    mapping::util::MappingType map_type,
                    std::vector<mapping::util::simple_hit> &accepted_hits,
                    phmap::flat_hash_map<uint64_t, uint32_t> &unmapped_bc_map,
                    uint32_t &num_reads_in_chunk, rad_writer &bw) {
  if (map_type == mapping::util::MappingType::UNMAPPED) {
    unmapped_bc_map[bck.word(0)] += 1;
    return;
  }

  const uint32_t barcode_len = bc_kmer_t::k();
  const uint32_t umi_len = umi_kmer_t::k();

  // number of mappings
  bw << static_cast<uint32_t>(accepted_hits.size());

  // bc
  // if we can fit the barcode into an integer
  if (barcode_len <= 32) {
    if (barcode_len <= 16) { // can use 32-bit int
      uint32_t shortbck =
        static_cast<uint32_t>(0x00000000FFFFFFFF & bck.word(0));
      bw << shortbck;
    } else { // must use 64-bit int
      bw << bck.word(0);
    }
  } else { // must use a string for the barcode
    std::cerr << "should not happen\n";
  }

  // umi
  if (umi_len <= 16) { // if we can use 32-bit int
    uint64_t umiint = umi.word(0);
    uint32_t shortumi = static_cast<uint32_t>(0x00000000FFFFFFFF & umiint);
    bw << shortumi;
  } else if (umi_len <= 32) { // if we can use 64-bit int
    uint64_t umiint = umi.word(0);
    bw << umiint;
  } else { // must use string
    std::cerr << "should not happen\n";
  }

  for (auto &aln : accepted_hits) {
    uint32_t fw_mask = 0x00000000;
    switch (map_type) {
    case mapping::util::MappingType::SINGLE_MAPPED:
      // single end case is the easiest, orientation is
      // whatever the read is.
      fw_mask = aln.is_fw ? 0x80000000 : 0x00000000;
      break;
    case mapping::util::MappingType::MAPPED_FIRST_ORPHAN:
      // if we mapped the first read in the pair
      // take it's orientation
      fw_mask = aln.is_fw ? 0x80000000 : 0x00000000;
      break;
    case mapping::util::MappingType::MAPPED_SECOND_ORPHAN:
      // if we mapped the second read in the pair
      // take the *opposite* of the mapped read's
      // orientation (representing what we'd expect
      // for read1 in a concordant pair)
      fw_mask = aln.is_fw ? 0x00000000 : 0x80000000;
      break;
    case mapping::util::MappingType::MAPPED_PAIR:
      // if the read is mapped in a proper pair, we
      // simply take the orientation of read 1 as
      // the orientation of the fragment
      fw_mask = aln.is_fw ? 0x80000000 : 0x00000000;
      break;
    case mapping::util::MappingType::UNMAPPED:
      // don't do anything here
      // NOTE: should not happen!
      break;
    }
    bw << (aln.tid | fw_mask);
  }
  ++num_reads_in_chunk;
}

inline void write_to_rad_stream_bulk(mapping::util::MappingType map_type,
                         std::vector<mapping::util::simple_hit> &accepted_hits,
                         uint32_t &num_reads_in_chunk, rad_writer &bw) {
  if (map_type == mapping::util::MappingType::UNMAPPED) {
    // do nothing here
    return;
  }

  // otherwise, we always write the number of mappings
  bw << static_cast<uint32_t>(accepted_hits.size());

  // then the read-level tag (fragment mapping type)
  bw << static_cast<uint8_t>(map_type);

  // for each fragment
  for (auto &aln : accepted_hits) {
    // top 2 bits are fw,rc ori
    uint32_t fw_mask = aln.is_fw ? 0x80000000 : 0x00000000;
    uint32_t mate_fw_mask = aln.mate_is_fw ? 0x40000000 : 0x00000000;
    // bottom 30 bits are target id
    bw << ((0x3FFFFFFF & aln.tid) | fw_mask | mate_fw_mask);

    int32_t leftmost_pos = 0;
    // placeholder value for no fragment length
    uint16_t frag_len = std::numeric_limits<uint16_t>::max();

    switch (map_type) {
    case mapping::util::MappingType::SINGLE_MAPPED:
      // then the posittion must be that of the only
      // mapped read.
      leftmost_pos = std::max(0, aln.pos);
      break;
    case mapping::util::MappingType::MAPPED_FIRST_ORPHAN:
      leftmost_pos = std::max(0, aln.pos);
      break;
    case mapping::util::MappingType::MAPPED_SECOND_ORPHAN:
      // it's not mate pos b/c in this case we
      // simply returned the right accepted hits
      // as the accepted hits
      leftmost_pos = std::max(0, aln.pos);
      break;
    case mapping::util::MappingType::MAPPED_PAIR:
      // if we actually have a paird fragment get the
      // leftmost position
      leftmost_pos = std::min(aln.pos, aln.mate_pos);
      frag_len = aln.frag_len();
      // if the leftmost position is < 0, then adjust
      // the overhang by setting the start position to 0
      // and subtracting the overhang from the fragment
      // length.
      if (leftmost_pos < 0) {
        frag_len = aln.frag_len() + leftmost_pos;
        leftmost_pos = 0;
      }
      break;
    case mapping::util::MappingType::UNMAPPED:
      // don't do anything here
      break;
    }

    bw << static_cast<uint32_t>(leftmost_pos);
    bw << frag_len;
  }
  ++num_reads_in_chunk;
}



inline void write_to_rad_stream_atac(bc_kmer_t& bck, mapping::util::MappingType map_type,
                                     std::vector<mapping::util::simple_hit>& accepted_hits,
                                     phmap::flat_hash_map<uint64_t, uint32_t>& unmapped_bc_map,
                                     uint32_t& num_reads_in_chunk, std::string& strbuff, 
                                     std::string& barcode, mindex::reference_index& ri, 
                                     RAD::RAD_Writer& rw, RAD::Token& token, bool tn5_shift) {
                                        
    if (map_type == mapping::util::MappingType::UNMAPPED) {
        unmapped_bc_map[bck.word(0)] += 1;
        // do nothing here
        return;
    }
    RAD::Read read_rec;
    
    read_rec.set(accepted_hits.size());
    
    const uint32_t barcode_len = bc_kmer_t::k();
    if (barcode_len <= 32) {
        if (barcode_len <= 16) {  // can use 32-bit int
            uint32_t shortbck = static_cast<uint32_t>(0x00000000FFFFFFFF & bck.word(0));
            read_rec.add_tag(RAD::Type::u32(shortbck));
        } else {  // must use 64-bit int
            read_rec.add_tag(RAD::Type::u64(bck.word(0)));
        }
    } else {
            std::cerr << "should not happen\n";
    }
    for (auto& aln : accepted_hits) {
        RAD::Aln_Record aln_rec;
        uint8_t type{0};
        // top 2 bits are fw,rc ori
        // uint32_t fw_mask = aln.is_fw ? 0x80000000 : 0x00000000;
        // uint32_t mate_fw_mask = aln.mate_is_fw ? 0x40000000 : 0x00000000;
        // bottom 30 bits are target id
        // strbuff += std::to_string((0x3FFFFFFF & aln.tid) | fw_mask | mate_fw_mask);
        strbuff += ri.ref_name(aln.tid);
        strbuff += "\t";
        int32_t leftmost_pos = 0;
        // placeholder value for no fragment length
        uint16_t frag_len = std::numeric_limits<uint16_t>::max();

        switch (map_type) {
            case mapping::util::MappingType::SINGLE_MAPPED:
                // then the posittion must be that of the only
                // mapped read.
                leftmost_pos = std::max(0, aln.pos);
                type = 1;
                break;
            case mapping::util::MappingType::MAPPED_FIRST_ORPHAN:
                leftmost_pos = std::max(0, aln.pos);
                type = 2;
                break;
            case mapping::util::MappingType::MAPPED_SECOND_ORPHAN:
                // it's not mate pos b/c in this case we
                // simply returned the right accepted hits
                // as the accepted hits
                leftmost_pos = std::max(0, aln.pos);
                type = 3;
                break;
            case mapping::util::MappingType::MAPPED_PAIR:
                // if we actually have a paird fragment get the
                // leftmost position
                leftmost_pos = std::min(aln.pos, aln.mate_pos);
                frag_len = aln.frag_len();
                type = 4;
                // if the leftmost position is < 0, then adjust
                // the overhang by setting the start position to 0
                // and subtracting the overhang from the fragment
                // length.
                if (leftmost_pos < 0) {
                    frag_len = aln.frag_len() + leftmost_pos;
                    leftmost_pos = 0;
                }
                break;
            case mapping::util::MappingType::UNMAPPED:
                type = 8;
                // don't do anything here
                break;
        }
        if (tn5_shift) {
          leftmost_pos += 4;
          frag_len -= 9;
        }
        aln_rec.clear();
        aln_rec.add_tag(RAD::Type::u32(aln.tid));
        aln_rec.add_tag(RAD::Type::u8(type));
        aln_rec.add_tag(RAD::Type::u32(leftmost_pos));
        aln_rec.add_tag(RAD::Type::u16(frag_len));
        read_rec.add_aln_rec(aln_rec);
        
        strbuff += std::to_string(leftmost_pos);
        strbuff += "\t";
        strbuff += std::to_string(leftmost_pos + frag_len);
        strbuff += "\t";
        strbuff += barcode;
        strbuff += "\t";
        strbuff += std::to_string(accepted_hits.size());
        strbuff += "\n";
    }
    
    rw.add(read_rec, token);
    ++num_reads_in_chunk;
}


}  // namespace util
}  // namespace rad


#endif  //__RAD_UTIL_HPP__
