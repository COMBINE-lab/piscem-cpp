#ifndef __SC_PROTOCOL__
#define __SC_PROTOCOL__

// hardcode 10x v3 for now
class sc_protocol {
public:
    // We'd really like an std::optional<string&> here, but C++17
    // said no to that.
    std::string* extract_bc(std::string& r1, std::string& r2) {
        return (r1.length() >= bc_len) ? (bc.assign(r1, 0, bc_len), &bc) : nullptr;
    }

    // We'd really like an std::optional<string&> here, but C++17
    // said no to that.
    std::string* extract_umi(std::string& r1, std::string& r2) {
        return (r1.length() >= (bc_len + umi_len)) ? (umi.assign(r1, bc_len, umi_len), &umi)
                                                   : nullptr;
    }

    std::string* extract_mappable_read(std::string& r1, std::string& r2) { return &r2; }

private:
    std::string umi;
    std::string bc;
    const size_t bc_len = 16;
    const size_t umi_len = 12;
};


#endif // __SC_PROTOCOL__
