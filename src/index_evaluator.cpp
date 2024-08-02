#include "../external/sshash/external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "../include/reference_index.hpp"
#include "../include/util_piscem.hpp"
#include "../external/sshash/include/util.hpp"
#include "../include/mapping/utils.hpp"
#include "../include/parallel_hashmap/phmap.h"
#include "spdlog_piscem/spdlog.h"

#include <iostream>
#include <vector>
#include <numeric>
#include <cstdio>

std::vector<uint64_t> histogram(mindex::reference_index& ri){
	const auto& ct = ri.get_contig_table();
	
	phmap::flat_hash_map<uint64_t, uint64_t> freq_map;

	auto cit = ct.m_ctg_offsets.at(0);// sshash::ef_sequence<false>::iterator(&ct.m_ctg_offsets, 0);
	uint64_t max_freq = 0;
	uint64_t prev = 0;
	while (cit.has_next()) {
		auto n = cit.next();
		auto v = (n-prev);
		auto& f = freq_map[v];
		f += 1;
		if (v > max_freq) { max_freq = v; }
		prev = n;
	}

	std::cerr << "max_freq = " << max_freq << "\n";
	std::vector<uint64_t> freqs(max_freq+1, 0);
	for (size_t i = 0; i <= max_freq; ++i){
		freqs[i] = freq_map[i];
	}
	return freqs;
}

int main(int argc, char* argv[]) {

	std::ios_base::sync_with_stdio(false);
  spdlog_piscem::set_level(spdlog_piscem::level::warn);
  cmd_line_parser::parser parser(argc, argv);

  /* mandatory arguments */
  parser.add("input_filename",
             "input index prefix.", "-i", true);
  parser.add("verbose",
             "output extra information from the index.", "-v", false, true);
 
  if (!parser.parse()) return 1;

	auto index_prefix = parser.get<std::string>("input_filename");

	mindex::reference_index ri(index_prefix);

  auto num_refs = ri.num_refs(); 
  if (parser.get<bool>("verbose")) {
    std::cout << "{\n";
    std::cout << "\"num_refs\": " << num_refs << ",\n";
    std::cout << "\"ref_entries\": {\n";
    for (size_t i = 0; i < num_refs; ++i) {
      auto& ni = ri.ref_name(i);
      std::cout << "\"" << ni << "\": " << ri.ref_len(i);
      if (i < num_refs - 1) { 
        std::cout << ",\n";
      }
    }
    std::cout << "},\n";
    std::cout << "histogram : [";
    auto freqs = histogram(ri);
    for (size_t i = 0; i < freqs.size(); ++i) {
      std::cout << freqs[i];
      if (i != freqs.size()-1) { std::cout << ", "; }
    }
    std::cout << "],\n";
    std::cout << "}\n"; 
  }

}
