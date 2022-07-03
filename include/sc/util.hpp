#ifndef __SC_UTIL_HPP__
#define __SC_UTIL_HPP__

namespace single_cell {

  namespace util {

    enum class BarCodeRecovered : uint8_t { OK, RECOVERED, NOT_RECOVERED };

    BarCodeRecovered recover_barcode(std::string& sequence) {
      size_t pos = sequence.find_first_not_of("ACTGactg");
      if (pos == std::string::npos) { return BarCodeRecovered::OK; }

      // Randomly assigning 'A' to first base with 'N'
      sequence[pos] = 'A';
      size_t invalid_pos = sequence.find_first_not_of("ACTGactg", pos);
      return (invalid_pos == std::string::npos) ? BarCodeRecovered::RECOVERED
        : BarCodeRecovered::NOT_RECOVERED;
    }

  }
}


#endif // __SC_UTIL_HPP__
