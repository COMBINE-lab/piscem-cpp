// Copyright(c) 2015-present, Gabi Melman & spdlog contributors.
// Distributed under the MIT License (http://opensource.org/licenses/MIT)

#pragma once

#ifndef SPDLOG_PISCEM_HEADER_ONLY
#    include <spdlog_piscem/details/log_msg.h>
#endif

#include <spdlog_piscem/details/os.h>

namespace spdlog_piscem {
namespace details {

SPDLOG_PISCEM_INLINE log_msg::log_msg(spdlog_piscem::log_clock::time_point log_time, spdlog_piscem::source_loc loc, string_view_t a_logger_name,
    spdlog_piscem::level::level_enum lvl, spdlog_piscem::string_view_t msg)
    : logger_name(a_logger_name)
    , level(lvl)
    , time(log_time)
#ifndef SPDLOG_PISCEM_NO_THREAD_ID
    , thread_id(os::thread_id())
#endif
    , source(loc)
    , payload(msg)
{}

SPDLOG_PISCEM_INLINE log_msg::log_msg(
    spdlog_piscem::source_loc loc, string_view_t a_logger_name, spdlog_piscem::level::level_enum lvl, spdlog_piscem::string_view_t msg)
    : log_msg(os::now(), loc, a_logger_name, lvl, msg)
{}

SPDLOG_PISCEM_INLINE log_msg::log_msg(string_view_t a_logger_name, spdlog_piscem::level::level_enum lvl, spdlog_piscem::string_view_t msg)
    : log_msg(os::now(), source_loc{}, a_logger_name, lvl, msg)
{}

} // namespace details
} // namespace spdlog_piscem
