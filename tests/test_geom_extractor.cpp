#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"
#include "../include/sc/util.hpp"


TEST_CASE("10x_v3_is_right") {

  std::string r1{"ACGGATGGTGAGGTTGCCGTATGAGAGG"};
  std::string r2{"ACGAGACATTTGCACATGACTGAATGGGGGGTTGTCCACACCTAGAGTTATAGGCATATGC"};

  chromium_v3 geom;

  CHECK(*geom.extract_bc(r1, r2) == r1.substr(0, 16));
  CHECK(*geom.extract_umi(r1, r2) == r1.substr(16, 12));
  auto mappable_seqs = geom.get_mappable_read_sequences(r1, r2);
  CHECK(mappable_seqs.seq1 == nullptr);
  CHECK(*mappable_seqs.seq2 == r2);
}

TEST_CASE("10x_v2") {

  std::string r1{"ACGGATGGTGAGGTTGCCGTATGAGAGG"};
  std::string r2{"ACGAGACATTTGCACATGACTGAATGGGGGGTTGTCCACACCTAGAGTTATAGGCATATGC"};

  chromium_v2 geom;

  CHECK(*geom.extract_bc(r1, r2) == r1.substr(0, 16));
  CHECK(*geom.extract_umi(r1, r2) == r1.substr(16, 10));
  auto mappable_seqs = geom.get_mappable_read_sequences(r1, r2);
  CHECK(mappable_seqs.seq1 == nullptr);
  CHECK(*mappable_seqs.seq2 == r2);
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
}


