#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"
#include "../include/sc/util.hpp"




TEST_CASE("10x_v2_is_right") {

  std::string r1{"ACGGATGGTGAGGTTGCCGTATGAGAGG"};
  std::string r2{"ACGAGACATTTGCACATGACTGAATGGGGGGTTGTCCACACCTAGAGTTATAGGCATATGC"};

  chromium_v2 geom;

  CHECK(*geom.extract_bc(r1, r2) == r1.substr(0, 16));
  CHECK(*geom.extract_umi(r1, r2) == r1.substr(16, 10));
  auto mappable_seqs = geom.get_mappable_read_sequences(r1, r2);
  CHECK(mappable_seqs.seq1 == nullptr);
  CHECK(*mappable_seqs.seq2 == r2);

  CHECK(mappable_seqs.get_reads_to_map() == ReadsToMap::SECOND);
  CHECK(geom.is_bio_paired_end() == false);
}

TEST_CASE("10x_v4_3p_is_right") {

  std::string r1{"ACGGATGGTGAGGTTGCCGTATGAGAGG"};
  std::string r2{"ACGAGACATTTGCACATGACTGAATGGGGGGTTGTCCACACCTAGAGTTATAGGCATATGC"};

  chromium_v4_3p geom;

  CHECK(*geom.extract_bc(r1, r2) == r1.substr(0, 16));
  CHECK(*geom.extract_umi(r1, r2) == r1.substr(16, 12));
  auto mappable_seqs = geom.get_mappable_read_sequences(r1, r2);
  CHECK(mappable_seqs.seq1 == nullptr);
  CHECK(*mappable_seqs.seq2 == r2);

  CHECK(mappable_seqs.get_reads_to_map() == ReadsToMap::SECOND);
  CHECK(geom.is_bio_paired_end() == false);
}

TEST_CASE("10x_v4_5p_is_right") {

  std::string r1{"ACGGATGGTGAGGTTGCCGTATGAGAGGACTACTGGTACAGGGTGTGTCATGGTGACGTGA"};
  std::string r2{"ACGAGACATTTGCACATGACTGAATGGGGGGTTGTCCACACCTAGAGTTATAGGCATATGC"};

  chromium_v4_5p geom;
  constexpr size_t tso_len = 13;

  CHECK(*geom.extract_bc(r1, r2) == r1.substr(0, 16));
  CHECK(*geom.extract_umi(r1, r2) == r1.substr(16, 12));
  auto mappable_seqs = geom.get_mappable_read_sequences(r1, r2);
  CHECK(*mappable_seqs.seq1 == r1.substr(28 + tso_len));
  CHECK(*mappable_seqs.seq2 == r2);

  CHECK(mappable_seqs.get_reads_to_map() == ReadsToMap::BOTH);
  CHECK(geom.is_bio_paired_end() == true);
}

TEST_CASE("10x_v3_is_right") {

  std::string r1{"ACGGATGGTGAGGTTGCCGTATGAGAGG"};
  std::string r2{"ACGAGACATTTGCACATGACTGAATGGGGGGTTGTCCACACCTAGAGTTATAGGCATATGC"};

  chromium_v3 geom;

  CHECK(*geom.extract_bc(r1, r2) == r1.substr(0, 16));
  CHECK(*geom.extract_umi(r1, r2) == r1.substr(16, 12));
  auto mappable_seqs = geom.get_mappable_read_sequences(r1, r2);
  CHECK(mappable_seqs.seq1 == nullptr);
  CHECK(*mappable_seqs.seq2 == r2);

  CHECK(mappable_seqs.get_reads_to_map() == ReadsToMap::SECOND);
  CHECK(geom.is_bio_paired_end() == false);
}

TEST_CASE("10x_v3_5p_is_right") {

  std::string r1{"ACGGATGGTGAGGTTGCCGTATGAGAGGACTACTGGTACAGGGTGTGTCATGGTGACGTGA"};
  std::string r2{"ACGAGACATTTGCACATGACTGAATGGGGGGTTGTCCACACCTAGAGTTATAGGCATATGC"};

  chromium_v3_5p geom;
  constexpr size_t tso_len = 13;

  CHECK(*geom.extract_bc(r1, r2) == r1.substr(0, 16));
  CHECK(*geom.extract_umi(r1, r2) == r1.substr(16, 12));
  auto mappable_seqs = geom.get_mappable_read_sequences(r1, r2);
  CHECK(*mappable_seqs.seq1 == r1.substr(28 + tso_len));
  CHECK(*mappable_seqs.seq2 == r2);

  CHECK(mappable_seqs.get_reads_to_map() == ReadsToMap::BOTH);
  CHECK(geom.is_bio_paired_end() == true);
}

TEST_CASE("custom_single_read_is_correct") {
  std::string geom_str = "1{b[18]u[14]x:}2{r:}";
  auto custom_geo_ptr = single_cell::util::parse_custom_geometry(geom_str);
  std::string r1{"ACGGATGGTGAGGTTGCCGTATGAGAGGACTACTGGTACAGGGTGTGTCATTCCGTGACTTTTTCACAG"};
  std::string r2{"ACGAGACATTTGCACATGACTGAATGGGGGGTTGTCCACACCTAGAGTTATAGGCATATGC"};

  CHECK(*custom_geo_ptr->extract_bc(r1, r2) == r1.substr(0, 18));
  CHECK(*custom_geo_ptr->extract_umi(r1, r2) == r1.substr(18, 14));
  auto mappable_seqs = custom_geo_ptr->get_mappable_read_sequences(r1, r2);
  CHECK(mappable_seqs.seq1 == nullptr);
  CHECK(*mappable_seqs.seq2 == r2);

  CHECK(mappable_seqs.get_reads_to_map() == ReadsToMap::SECOND);
  CHECK(custom_geo_ptr->is_bio_paired_end() == false);
}

TEST_CASE("custom_pe_read_is_correct") {
  std::string geom_str = "1{b[18]u[14]x[13]r:}2{r:}";
  auto custom_geo_ptr = single_cell::util::parse_custom_geometry(geom_str);
  std::string r1{"ACGGATGGTGAGGTTGCCGTATGAGAGGACTACTGGTACAGGGTGTGTCATTCCGTGACTTTTTCACAG"};
  std::string r2{"ACGAGACATTTGCACATGACTGAATGGGGGGTTGTCCACACCTAGAGTTATAGGCATATGC"};

  constexpr size_t tso_len = 13;

  CHECK(*custom_geo_ptr->extract_bc(r1, r2) == r1.substr(0, 18));
  CHECK(*custom_geo_ptr->extract_umi(r1, r2) == r1.substr(18, 14));
  auto mappable_seqs = custom_geo_ptr->get_mappable_read_sequences(r1, r2);
  CHECK(*mappable_seqs.seq1 == r1.substr(32 + tso_len));
  CHECK(*mappable_seqs.seq2 == r2);

  CHECK(mappable_seqs.get_reads_to_map() == ReadsToMap::BOTH);
  CHECK(custom_geo_ptr->is_bio_paired_end() == true);
}


TEST_CASE("custom_just_r1_is_correct") {
  std::string geom_str = "1{b[18]u[7]r:}2{u[7]x:}";
  auto custom_geo_ptr = single_cell::util::parse_custom_geometry(geom_str);
  std::string r1{"ACGGATGGTGAGGTTGCCGTATGAGAGGACTACTGGTACAGGGTGTGTCATTCCGTGACTTTTTCACAG"};
  std::string r2{"ACGAGACATTTGCACATGACTGAATGGGGGGTTGTCCACACCTAGAGTTATAGGCATATGC"};

  CHECK(custom_geo_ptr->validate());
  CHECK(*custom_geo_ptr->extract_bc(r1, r2) == r1.substr(0, 18));
  CHECK(*custom_geo_ptr->extract_umi(r1, r2) == r1.substr(18, 7) + r2.substr(0, 7));
  auto mappable_seqs = custom_geo_ptr->get_mappable_read_sequences(r1, r2);
  CHECK(*mappable_seqs.seq1 == r1.substr(25));
  CHECK(mappable_seqs.seq2 == nullptr);

  CHECK(mappable_seqs.get_reads_to_map() == ReadsToMap::FIRST);
  CHECK(custom_geo_ptr->is_bio_paired_end() == false);
}

TEST_CASE("invalid_geo_is_null") {
  std::string geom_str = "1{b[18]u[7]r:}2{u[7]:}";
  auto custom_geo_ptr = single_cell::util::parse_custom_geometry(geom_str);
  CHECK(custom_geo_ptr == nullptr);
}
