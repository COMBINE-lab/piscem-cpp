#ifndef __RAD_WRITER__
#define __RAD_WRITER__

#include <vector>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <cstring>

/**
 * low-level interface for writing information direclty into a binary
 * stream (vector of char).
 */
class rad_writer {
private:
    std::vector<char> _bin_data;
    static constexpr const uint16_t max_str_len{std::numeric_limits<uint16_t>::max()};

public:
    rad_writer(){};
    rad_writer(size_t reserve_size) { _bin_data.reserve(reserve_size); }
    // copy bin data to this record
    rad_writer(const std::vector<char>& bin_data) : _bin_data(bin_data){};
    // or just move if possible
    rad_writer(std::vector<char>&& bin_data) : _bin_data(move(bin_data)){};

    void clear() { _bin_data.clear(); }

    template <typename T>
    bool write_integer_at_offset(size_t offset, const T& v) {
        if (offset + sizeof(v) >= _bin_data.size()) { return false; }
        std::memcpy(_bin_data.data() + offset, &v, sizeof(v));
        return true;
    }

    rad_writer& operator<<(const bool& inval) {
        char* in_char_ptr = const_cast<char*>(reinterpret_cast<const char*>(&inval));
        std::copy(in_char_ptr, in_char_ptr + sizeof(inval), std::back_inserter(_bin_data));
        return *this;
    }
    rad_writer& operator<<(const uint8_t& inval) {
        char* in_char_ptr = const_cast<char*>(reinterpret_cast<const char*>(&inval));
        std::copy(in_char_ptr, in_char_ptr + sizeof(inval), std::back_inserter(_bin_data));
        return *this;
    }
    rad_writer& operator<<(const uint16_t& inval) {
        char* in_char_ptr = const_cast<char*>(reinterpret_cast<const char*>(&inval));
        std::copy(in_char_ptr, in_char_ptr + sizeof(inval), std::back_inserter(_bin_data));
        return *this;
    }
    rad_writer& operator<<(const uint32_t& inval) {
        char* in_char_ptr = const_cast<char*>(reinterpret_cast<const char*>(&inval));
        std::copy(in_char_ptr, in_char_ptr + sizeof(inval), std::back_inserter(_bin_data));
        return *this;
    }
    rad_writer& operator<<(const uint64_t& inval) {
        char* in_char_ptr = const_cast<char*>(reinterpret_cast<const char*>(&inval));
        std::copy(in_char_ptr, in_char_ptr + sizeof(inval), std::back_inserter(_bin_data));
        return *this;
    }
    rad_writer& operator<<(const int32_t& inval) {
        char* in_char_ptr = const_cast<char*>(reinterpret_cast<const char*>(&inval));
        std::copy(in_char_ptr, in_char_ptr + sizeof(inval), std::back_inserter(_bin_data));
        return *this;
    }

    rad_writer& operator<<(const double& inval) {
        char* in_char_ptr = const_cast<char*>(reinterpret_cast<const char*>(&inval));
        std::copy(in_char_ptr, in_char_ptr + sizeof(inval), std::back_inserter(_bin_data));
        return *this;
    }
    rad_writer& operator<<(const std::string& inval) {
        if (inval.size() >= max_str_len) {
            std::cerr << "ERROR!! DOESN'T SUPPORT STRING LENGTH LONGER THAN " << max_str_len
                      << ". String length: " << inval.size() << "\n";
            std::exit(1);
        }
        (*this) << static_cast<uint16_t>(inval.size());
        //_bin_data.push_back(tmp);
        //(*this) << inval.size();
        // std::cout << inval.size() << " " << inval.c_str()  << " " << inval <<
        // "\n";
        char* in_char_ptr = const_cast<char*>(inval.c_str());
        std::copy(in_char_ptr, in_char_ptr + inval.size(), std::back_inserter(_bin_data));
        return *this;
    }

#ifdef STX_NO_STD_STRING_VIEW
/*
rad_writer& operator<<(const stx::string_view& inval) {
  if (inval.size() >= max_str_len) {
    std::cerr << "ERROR!! DOESN'T SUPPORT STRING LENGTH LONGER THAN 255. "
                 "String length: "
              << inval.size() << "\n";
    std::exit(1);
  }
  (*this) << static_cast<uint16_t>(inval.size());
  //(*this) << inval.size();
  // std::cout << inval.size() << " " << inval.data()  << " " << inval <<
  // "\n";
  char* in_char_ptr = const_cast<char*>(inval.data());
  std::copy(in_char_ptr, in_char_ptr + inval.size(),
            std::back_inserter(_bin_data));
  return *this;
}
*/
#endif

    uint64_t num_bytes() { return _bin_data.size(); }

    // support for logging directly from spdlog
    template <typename OStream>
    friend OStream& operator<<(OStream& os, const rad_writer& bin_record) {
        std::ostream_iterator<char> out_iter(os);
        std::copy(bin_record._bin_data.begin(), bin_record._bin_data.end(), out_iter);
        return os;
    }
};

template bool rad_writer::write_integer_at_offset<uint16_t>(size_t offset, const uint16_t& v);

template bool rad_writer::write_integer_at_offset<uint32_t>(size_t offset, const uint32_t& v);

template bool rad_writer::write_integer_at_offset<uint64_t>(size_t offset, const uint64_t& v);

#endif  //__RAD_WRITER__
