#include "../external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "../include/reference_index.hpp"
#include "../include/util.hpp"
#include "../include/mapping/utils.hpp"
#include "../include/parallel_hashmap/phmap.h"

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
    cmd_line_parser::parser parser(argc, argv);

    /* mandatory arguments */
    parser.add("input_filename",
               "input index prefix.");
    if (!parser.parse()) return 1;

	auto index_prefix = parser.get<std::string>("input_filename");

	mindex::reference_index ri(index_prefix);

	auto freqs = histogram(ri);
	for (size_t i = 0; i < freqs.size(); ++i) {
		std::cout << freqs[i];
		if (i != freqs.size()-1) { std::cout << "\t"; }
	}
	std::cout << "\n";
}
