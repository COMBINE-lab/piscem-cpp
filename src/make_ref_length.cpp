#include "../include/reference_index.hpp"
#include "../include/cli11/CLI11.hpp"
#include <atomic>
#include <chrono>
#include <iostream>
#include <vector>

void write_ref_len(mindex::reference_index& ri,
	std::string& output_path) {
	std::string temp = "";
	ghc::filesystem::path o_path(output_path);
	std::ofstream ref_file = o_path / "ref.tsv";
	// ghc::filesystem::path 
	for(auto i = 0; i < ri.num_refs(); i++) {
        temp += ri.ref_name(i) +  "\t" + std::to_string(ri.ref_len(i)) + "\n";
    }
	ref_file << temp;
	ref_file.close();
}

int main(int argc, char** argv) {
	std::string index_basename;
	std::string output_dir;
	CLI::App app{"writing reference length"};
    	app.add_option("-i,--index", index_basename, "input index prefix")->required();
		app.add_option("-o,--output", output_dir, "output dir")->required();
	CLI11_PARSE(app, argc, argv);
	mindex::reference_index ri(index_basename, false);
	write_ref_len(ri, output_dir);
	return 0;
}

