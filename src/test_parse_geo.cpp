#include <string>
#include <utility>
#include <iostream>
#include "../include/itlib/small_vector.hpp"
#include "../include/peglib.h"
#include "../include/cli11/CLI11.hpp"

enum class geo_tag_type { BC, UMI, READ, FIXED, DISCARD };
struct geo_part {
    geo_tag_type ttype{geo_tag_type::DISCARD};
    int64_t len{-1};
    friend std::ostream& operator<<(std::ostream& os, const geo_part& gp);
};

std::ostream& operator<<(std::ostream& os, const geo_part& gp) {
    switch (gp.ttype) {
        case geo_tag_type::BC: {
            os << "BC [" << gp.len << "]";
            break;
        }
        case geo_tag_type::UMI: {
            os << "UMI [" << gp.len << "]";
            break;
        }
        case geo_tag_type::READ: {
            os << "R [" << gp.len << "]";
            break;
        }
        case geo_tag_type::FIXED: {
            os << "F [" << gp.len << "]";
            break;
        }
        case geo_tag_type::DISCARD: {
            os << "X [" << gp.len << "]";
            break;
        }
    }
    return os;
}

struct protocol_state {
    geo_tag_type curr_geo_type{geo_tag_type::DISCARD};
    size_t current_read{0};
    std::pair<int64_t, int64_t> length_range;
    std::vector<int64_t> lengths;
    std::vector<geo_part> geo_parts;
    std::vector<geo_part> geo_parts_r1;
    std::vector<geo_part> geo_parts_r2;
};

struct str_slice {
    int32_t offset{-1};
    int32_t len{-1};
};

class protocol {
public:
    protocol(protocol_state& ps) {
        // convert the list of lengths to a
        // list of offsets and lengths
        size_t current_offset{0};
        for (auto& gp : ps.geo_parts_r1) {
            if (gp.ttype == geo_tag_type::BC) {
                bc_slices_r1.push_back({current_offset, gp.len});
            } else if (gp.ttype == geo_tag_type::UMI) {
                umi_slices_r1.push_back({current_offset, gp.len});
            } else if (gp.ttype == geo_tag_type::READ) {
                read_slices_r1.push_back({current_offset, gp.len});
            }
            current_offset += gp.len;
        }
        current_offset = 0;
        for (auto& gp : ps.geo_parts_r2) {
            if (gp.ttype == geo_tag_type::BC) {
                bc_slices_r2.push_back({current_offset, gp.len});
            } else if (gp.ttype == geo_tag_type::UMI) {
                umi_slices_r2.push_back({current_offset, gp.len});
            } else if (gp.ttype == geo_tag_type::READ) {
                read_slices_r2.push_back({current_offset, gp.len});
            }
            current_offset += gp.len;
        }

        // bc and umi lengths must be
        // bounded in parsing, so we can
        // accumulate them this way
        int64_t bc_len{0};
        for (auto& p : bc_slices_r1) { bc_len += p.len; }
        for (auto& p : bc_slices_r2) { bc_len += p.len; }

        int64_t umi_len{0};
        for (auto& p : umi_slices_r1) { umi_len += p.len; }
        for (auto& p : umi_slices_r2) { umi_len += p.len; }

        // read lengths can be unbounded
        int64_t read_len{0};
        for (auto& p : read_slices_r1) {
            if (p.len <= 0) {
                read_len = -1;
                break;
            }
            read_len += p.len;
        }
        if (read_len > 0) {
            for (auto& p : read_slices_r2) {
                if (p.len <= 0) {
                    read_len = -1;
                    break;
                }
                read_len += p.len;
            }
        }

        std::cerr << "bc_len = " << bc_len << ", umi_len = " << umi_len
                  << ", read_len = " << read_len << "\n";
    }

    bool validate() {
        bool valid = has_barcode;
        valid = (valid and has_umi);
        valid = (valid and has_biological_read);
    }

    // We'd really like an std::optional<string&> here, but C++17
    // said no to that.
    std::string* extract_bc(std::string& r1, std::string& r2) {
        bc_buffer.clear();
        const auto r1_len = r1.size();
        const auto r2_len = r2.size();
        for (auto& bp : bc_slices_r1) {
            if (bp.offset + bp.len > r1_len) { return nullptr; }
            bc_buffer.append(r1, bp.offset, bp.len);
        }
        for (auto& bp : bc_slices_r2) {
            if (bp.offset + bp.len > r2_len) { return nullptr; }
            bc_buffer.append(r2, bp.offset, bp.len);
        }
        return &bc_buffer;
    }

    // We'd really like an std::optional<string&> here, but C++17
    // said no to that.
    std::string* extract_umi(std::string& r1, std::string& r2) {
        umi_buffer.clear();
        const auto r1_len = r1.size();
        const auto r2_len = r2.size();
        for (auto& up : umi_slices_r1) {
            if (up.offset + up.len > r1_len) { return nullptr; }
            umi_buffer.append(r1, up.offset, up.len);
        }
        for (auto& up : umi_slices_r2) {
            if (up.offset + up.len > r2_len) { return nullptr; }
            umi_buffer.append(r2, up.offset, up.len);
        }
        return &umi_buffer;
    }

    std::string* extract_mappable_read(std::string& r1, std::string& r2) { 
      bool uses_r2 = !read_slices_r2.empty();
      bool uses_r1 = !read_slices_r1.empty();
      if ( uses_r2 and read_slices_r2[0].len == -1) {
        return &r2;
      } else if ( uses_r1 and read_slices_r1[0].len == -1) {
        return &r1;
      } else if ( uses_r2 ) {
        read_buffer.clear();
        size_t r2_len = r2.size();
        for (auto& up : read_slices_r2) {
            if (up.offset + up.len > r2_len) { return nullptr; }
            size_t len = (up.len == -1) ? std::string::npos : up.len;
            umi_buffer.append(r2, up.offset, len);
        }
      } else if ( uses_r1 ) {
        read_buffer.clear();
        size_t r1_len = r1.size();
        for (auto& up : read_slices_r1) {
            if (up.offset + up.len > r1_len) { return nullptr; }
            size_t len = (up.len == -1) ? std::string::npos : up.len;
            umi_buffer.append(r1, up.offset, len);
        }
      }
      return &read_buffer;
    }

    void print() {
        std::cerr << "barcode:\n";
        if (!bc_slices_r1.empty()) { std::cerr << "\t1 "; }
        for (auto& s : bc_slices_r1) { std::cerr << s.offset << ":" << s.len << " "; }
        std::cerr << "\n";
        if (!bc_slices_r2.empty()) { std::cerr << "\t2 "; }
        for (auto& s : bc_slices_r2) { std::cerr << s.offset << ":" << s.len << " "; }

        std::cerr << "umi:\n";
        if (!umi_slices_r1.empty()) { std::cerr << "\t1 "; }
        for (auto& s : umi_slices_r1) { std::cerr << s.offset << ":" << s.len << " "; }
        std::cerr << "\n";
        if (!umi_slices_r2.empty()) { std::cerr << "\t2 "; }
        for (auto& s : umi_slices_r2) { std::cerr << s.offset << ":" << s.len << " "; }

        std::cerr << "read:\n";
        if (!read_slices_r1.empty()) { std::cerr << "\t1 "; }
        for (auto& s : read_slices_r1) { std::cerr << s.offset << ":" << s.len << " "; }
        std::cerr << "\n";
        if (!read_slices_r2.empty()) { std::cerr << "\t2 "; }
        for (auto& s : read_slices_r2) { std::cerr << s.offset << ":" << s.len << " "; }
    }

private:
    bool has_biological_read{false};
    bool has_umi{false};
    bool has_barcode{false};
    std::string bc_buffer;
    std::string umi_buffer;
    std::string read_buffer;

    itlib::small_vector<str_slice, 8> bc_slices_r1;
    itlib::small_vector<str_slice, 8> umi_slices_r1;
    itlib::small_vector<str_slice, 8> read_slices_r1;

    itlib::small_vector<str_slice, 8> bc_slices_r2;
    itlib::small_vector<str_slice, 8> umi_slices_r2;
    itlib::small_vector<str_slice, 8> read_slices_r2;
};

int main(int argc, char* argv[]) {
    std::string geom_string;
    CLI::App app{"Test geometry parsing"};
    app.add_option("-g,--custom-geometry", geom_string, "geometry string")->required();
    CLI11_PARSE(app, argc, argv);

    peg::parser parser(R"(
Specification <- Read1Description Read2Description
Read1Description <- '1{' (UnboundedDescription / (BoundedDescription{1,10} UnboundedDescription{0,1})) '}'
Read2Description <- '2{' (UnboundedDescription / (BoundedDescription{1,10} UnboundedDescription{0,1})) '}'
BoundedDescription <- Barcode'['Lengths']' / UMI'['Lengths']' / Fixed'['Sequence']' / Discard'['Lengths']' / Read'['Lengths']'
UnboundedDescription <-  Discard':' / Read':'
Barcode <-  'b'
UMI <- 'u'
Discard <- 'x'
Fixed <- 'f'
Read <- 'r'
Sequence <- [ATGC]+
Lengths <- (Length '-' Length) / Length
Length <- [1-9][0-9]*
)");

    parser["Read1Description"] = [](const peg::SemanticValues& sv, std::any& dt) {
        std::cerr << "read 1 description :: ";
        std::cerr << "type name : " << dt.type().name() << std::endl;
        std::cerr << "line_info : " << sv.line_info().first << ", " << sv.line_info().second
                  << "\n";
        auto& ps = *std::any_cast<protocol_state*>(dt);
        ps.geo_parts_r1 = ps.geo_parts;
        ps.geo_parts.clear();
    };

    parser["Read2Description"] = [](const peg::SemanticValues& sv, std::any& dt) {
        std::cerr << "read 2 description :: ";
        std::cerr << "type name : " << dt.type().name() << std::endl;
        std::cerr << "line_info : " << sv.line_info().first << ", " << sv.line_info().second
                  << "\n";
        auto& ps = *std::any_cast<protocol_state*>(dt);
        ps.geo_parts_r2 = ps.geo_parts;
        ps.geo_parts.clear();
    };

    parser["UnboundedDescription"] = [](const peg::SemanticValues& sv, std::any& dt) {
        std::cerr << "unbounded description :: ";
        std::cerr << "type name : " << dt.type().name() << std::endl;
        std::cerr << "line_info : " << sv.line_info().first << ", " << sv.line_info().second
                  << "\n";
        auto& ps = *std::any_cast<protocol_state*>(dt);
        switch (sv.choice()) {
            case 0: {  // x
                ps.geo_parts.push_back({geo_tag_type::DISCARD, -1});
            } break;
            case 1: {  // r
                ps.geo_parts.push_back({geo_tag_type::READ, -1});
            } break;
            default:
                break;
        }
    };

    parser["BoundedDescription"] = [](const peg::SemanticValues& sv, std::any& dt) {
        std::cerr << "bounded description :: ";
        std::cerr << "type name : " << dt.type().name() << std::endl;
        std::cerr << "line_info : " << sv.line_info().first << ", " << sv.line_info().second
                  << "\n";
        auto& ps = *std::any_cast<protocol_state*>(dt);
        switch (sv.choice()) {
            case 0: {  // b
                std::cerr << "dtype :: barcode\n";
                ps.geo_parts.push_back({geo_tag_type::BC, 0});
                ps.geo_parts.back().len = std::any_cast<int64_t>(sv[1]);
            } break;
            case 1: {  // u
                std::cerr << "dtype :: umi\n";
                ps.geo_parts.push_back({geo_tag_type::UMI, 0});
                ps.geo_parts.back().len = std::any_cast<int64_t>(sv[1]);
            } break;
            case 2: {  // f
                std::cerr << "dtype :: fixed\n";
                ps.geo_parts.push_back({geo_tag_type::FIXED, 0});
                ps.geo_parts.back().len = std::any_cast<int64_t>(sv[1]);
            } break;
            case 3: {  // x
                std::cerr << "dtype :: discard\n";
                ps.geo_parts.push_back({geo_tag_type::DISCARD, 0});
                ps.geo_parts.back().len = std::any_cast<int64_t>(sv[1]);
            } break;
            case 4: {  // r
                std::cerr << "dtype :: read\n";
                ps.geo_parts.push_back({geo_tag_type::READ, 0});
                ps.geo_parts.back().len = std::any_cast<int64_t>(sv[1]);
            } break;
            default:
                break;
        }
    };

    parser["Lengths"] = [](const peg::SemanticValues& sv, std::any& dt) {
        std::cerr << "LENGTHS\n";
        std::cerr << "type name : " << dt.type().name() << std::endl;
        std::cerr << "line_info : " << sv.line_info().first << ", " << sv.line_info().second
                  << "\n";
        auto& ps = *std::any_cast<protocol_state*>(dt);
        int64_t len{-1};
        switch (sv.choice()) {
            case 0: {
                // uint64_t min_len = sv[0].token_to_number<int64_t>();
                // uint64_t max_len = sv[1].token_to_number<int64_t>();
                std::cerr << "variable length barcodes are not currently supported.\n";
            } break;
            default: {
                len = sv.token_to_number<int64_t>();
            } break;
        }
        return len;
    };
    parser["Length"] = [](const peg::SemanticValues& sv, std::any& dt) {
        std::cerr << "LENGTH\n";
        std::cerr << "type name : " << dt.type().name() << std::endl;
        std::cerr << "line_info : " << sv.line_info().first << ", " << sv.line_info().second
                  << "\n";
        return sv.token_to_number<int64_t>();
    };

    protocol_state ps;
    std::any dt = &ps;
    auto r = parser.parse(geom_string, dt);
    if (parser.parse(geom_string, dt)) {
        std::cout << "success!\n";
        std::cerr << "READ 1 : ";
        for (auto& gp : ps.geo_parts_r1) { std::cerr << gp << ", "; }
        std::cerr << "\n";
        std::cerr << "READ 2 : ";
        for (auto& gp : ps.geo_parts_r2) { std::cerr << gp << ", "; }
        std::cerr << "\n";
        protocol p(ps);
        p.print();
    } else {
        std::cerr << "description invalid!\n";
    }
}

/*
     parser["Length"] = [](const peg::SemanticValues &sv)
    {
      auto val = static_cast<size_t>(std::stoull(sv.token()));
      if(val<=0){
        std::cerr << "Lengths should be > 0. Exiting." << std::endl;
        exit(1);
      };
      return val;
    };


    parser["Read1Description"] = [&](const peg::SemanticValues &sv)
    {
        auto val = std::stoi(sv.token());
        proto.readNumber = val;
    };

    parser["ReadNumber2"] = [&](const peg::SemanticValues &sv)
    {
        auto val = std::stoi(sv.token());
        proto.readNumber = val;
        nPatterns = 1;
    };

    parser["Description"] = [&](const peg::SemanticValues &sv)
    {
        auto val = std::string(sv.token());
        if(val == "r") {
            if (!customGeo.bioReadFound){
                customGeo.bioRead = proto.readNumber;
                customGeo.bioReadFound = true;
            } else {
              // aopt.jointLog->error("Only contigous biological read expected.\nExiting now.");
              std::cerr << "Only contigous biological read expected.\nExiting now." << std::endl;
              exit(1);
            }
            modifyRegex(proto.readNumber, customGeo.reg, nPatterns, customGeo.bioPat);
        }
    };

    parser["Type"] = [&](const peg::SemanticValues &sv)
    {
        auto val = std::string(sv.token())[0];
        proto.type = static_cast<customReadpartType>(val);
    };

    parser["Fixed"] = [&](const peg::SemanticValues &sv)
    {
        auto val = std::string(sv.token())[0];
        proto.type = static_cast<customReadpartType>(val);
    };

    parser["Sequence"] = [&](const peg::SemanticValues &sv)
    {
        auto val = std::string(sv.token());
        modifyRegex(proto.readNumber, val, customGeo.reg, nPatterns);
    };

    parser["Lengths"] = [&](const peg::SemanticValues &sv)
    {
        switch (sv.choice())
        {
        case 0:
        {
            auto val = std::make_pair(
                peg::any_cast<size_t>(sv[0]), peg::any_cast<size_t>(sv[1]));
            if(val.second <= val.first) {
                //  aopt.jointLog->error("In length range [X-Y], Y should be > X.\nExiting now");
                std::cerr << "In length range [X-Y], Y should be > X.\nExiting now" << std::endl;
                exit(1);
            }
            if (proto.type == customReadpartType::bc) {
                customGeo.minBcLen += val.first;
                customGeo.maxBcLen += val.second;
                modifyRegex(proto.readNumber, proto.type, customGeo.reg, customGeo.b, nPatterns,
   val.first, val.second); } else if (proto.type == customReadpartType::umi) { customGeo.minUmiLen
   += val.first; customGeo.maxUmiLen += val.second; modifyRegex(proto.readNumber, proto.type,
   customGeo.reg, customGeo.u, nPatterns, val.first, val.second); } else {
                modifyRegex(proto.readNumber, proto.type, customGeo.reg, customGeo.u, nPatterns,
   val.first, val.second); // case 'x'
            }
        }
        break;
        case 1:
        {
            auto val = peg::any_cast<size_t>(sv[0]);
            if (proto.type == customReadpartType::bc)  {
                customGeo.minBcLen += val; // a fixed length bc increases the min length too, not
   updating can cause logical error when there are 2 pos of bcs customGeo.maxBcLen += val;
                modifyRegex(proto.readNumber, proto.type, customGeo.reg, customGeo.b, nPatterns,
   val); } else if (proto.type == customReadpartType::umi){ customGeo.minUmiLen += val;
                customGeo.maxUmiLen += val;
                modifyRegex(proto.readNumber, proto.type, customGeo.reg, customGeo.u, nPatterns,
   val); } else { modifyRegex(proto.readNumber, proto.type, customGeo.reg, customGeo.u, nPatterns,
   val); // case 'x'
            }
        }
        break;
        }
    };
*/
