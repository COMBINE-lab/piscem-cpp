#include "../include/libradicl/RAD_Writer.hpp"
#include "../include/libradicl/Alignment_Record.hpp"
#include "../include/libradicl/Read_Record.hpp"
#include "../include/libradicl/Byte_Array.hpp"
#include "../include/libradicl/Tags.hpp"
#include "../include/libradicl/Header.hpp"
#include "../src/libradicl/Header.cpp"
#include "../src/libradicl/Tags.cpp"
#include "../src/libradicl/RAD_Writer.cpp"
#include "../src/libradicl/Byte_Array.cpp"
#include "../src/libradicl/Buffer.cpp"
#include "../include/reference_index.hpp"

int main() {
  std::string index_basename = "/fs/cbcb-lab/rob/students/noor/Atacseq/piscem_analysis/hg38_ind_k23/hg38_ind_k23";
  mindex::reference_index ri(index_basename);
  std::vector<std::string> refs;
  for (size_t i = 0; i < ri.num_refs(); ++i) { refs.emplace_back(ri.ref_name(i)); }

  const RAD::Header header(true, refs.size(), refs);
  RAD::Tag_Defn tag_defn;

  tag_defn.add_file_tag<RAD::Type::u32>("ref_lengths");
  tag_defn.add_file_tag<RAD::Type::str>("cblen");

  tag_defn.add_read_tag<RAD::Type::u16>("barcode");

  tag_defn.add_aln_tag<RAD::Type::u32>("start_pos");
  tag_defn.add_aln_tag<RAD::Type::u16>("frag_len");
  tag_defn.add_aln_tag<RAD::Type::u32>("ref_name");

  RAD::RAD_Writer rw(header, tag_defn, "out.rad");

  rw.close();

  return 0;
}