#include "../external/sshash/src/dictionary.cpp"
#include "../external/sshash/src/builder/build.cpp"
#include "../external/sshash/src/info.cpp"
#include "../external/sshash/include/kmer.hpp"

namespace sshash {
  template struct dictionary<dna_uint_kmer_t<uint64_t>, decoded_offsets>;
}


