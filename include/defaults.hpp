#ifndef __PISCEM_DEFAULTS__
#define __PISCEM_DEFAULTS__

namespace piscem {

  namespace defaults {
    constexpr size_t nthread = 16;
    const std::string skipping_strat = "permissive";
    constexpr uint32_t max_ec_card = 4096;
    constexpr uint32_t max_hit_occ = 256;
    constexpr uint32_t max_hit_occ_recover = 1024;
    constexpr uint32_t max_read_occ = 2500;
  }

}

#endif // __PISCEM_DEFAULTS__
