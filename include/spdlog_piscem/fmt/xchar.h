//
// Copyright(c) 2016 Gabi Melman.
// Distributed under the MIT License (http://opensource.org/licenses/MIT)
//

#pragma once
//
// include bundled or external copy of fmtlib's xchar support
//

#if !defined(SPDLOG_PISCEM_USE_STD_FORMAT)
#    if !defined(SPDLOG_PISCEM_FMT_EXTERNAL)
#        ifdef SPDLOG_PISCEM_HEADER_ONLY
#            ifndef FMT_HEADER_ONLY
#                define FMT_HEADER_ONLY
#            endif
#        endif
#        include <spdlog_piscem/fmt/bundled/xchar.h>
#    else
#        include <fmt/xchar.h>
#    endif
#endif
