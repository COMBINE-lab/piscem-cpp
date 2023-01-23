#ifndef __RAD_HEADER__
#define __RAD_HEADER__

#include "rad_writer.hpp"
#include <cinttypes>
#include <cstring>

class rad_header {
private:
    bool is_paired_{false};
    uint64_t ref_count_{0};
    std::vector<std::string> ref_names;
    uint64_t num_chunks_{0};

public:
    // bool from_file();

    // adds n to the list of reference names and
    // returns the total number of reference names
    inline uint64_t add_refname(const std::string& n) {
        ref_names.emplace_back(n);
        ++ref_count_;
        return ref_count_;
    }

    inline bool is_paired() const { return is_paired_; }
    inline void is_paired(bool ip) { is_paired_ = ip; }

    inline bool dump_to_bin(rad_writer& bw) {
        uint8_t p = is_paired_ ? 1 : 0;
        bw << p;
        bw << ref_count_;
        for (auto& n : ref_names) { bw << n; }
        bw << num_chunks_;
        return true;
    }
};

#endif  // __RAD_HEADER__
