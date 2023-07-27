// Copyright(c) 2015-present, Gabi Melman & spdlog contributors.
// Distributed under the MIT License (http://opensource.org/licenses/MIT)

#pragma once

#ifndef SPDLOG_PISCEM_HEADER_ONLY
#    include <spdlog_piscem/common.h>
#endif

#include <algorithm>
#include <iterator>

namespace spdlog_piscem {
namespace level {

#if __cplusplus >= 201703L
constexpr
#endif
    static string_view_t level_string_views[] SPDLOG_PISCEM_LEVEL_NAMES;

static const char *short_level_names[] SPDLOG_PISCEM_SHORT_LEVEL_NAMES;

SPDLOG_PISCEM_INLINE const string_view_t &to_string_view(spdlog_piscem::level::level_enum l) SPDLOG_PISCEM_NOEXCEPT
{
    return level_string_views[l];
}

SPDLOG_PISCEM_INLINE const char *to_short_c_str(spdlog_piscem::level::level_enum l) SPDLOG_PISCEM_NOEXCEPT
{
    return short_level_names[l];
}

SPDLOG_PISCEM_INLINE spdlog_piscem::level::level_enum from_str(const std::string &name) SPDLOG_PISCEM_NOEXCEPT
{
    auto it = std::find(std::begin(level_string_views), std::end(level_string_views), name);
    if (it != std::end(level_string_views))
        return static_cast<level::level_enum>(std::distance(std::begin(level_string_views), it));

    // check also for "warn" and "err" before giving up..
    if (name == "warn")
    {
        return level::warn;
    }
    if (name == "err")
    {
        return level::err;
    }
    return level::off;
}
} // namespace level

SPDLOG_PISCEM_INLINE spdlog_ex::spdlog_ex(std::string msg)
    : msg_(std::move(msg))
{}

SPDLOG_PISCEM_INLINE spdlog_ex::spdlog_ex(const std::string &msg, int last_errno)
{
#ifdef SPDLOG_PISCEM_USE_STD_FORMAT
    msg_ = std::system_error(std::error_code(last_errno, std::generic_category()), msg).what();
#else
    memory_buf_t outbuf;
    fmt::format_system_error(outbuf, last_errno, msg.c_str());
    msg_ = fmt::to_string(outbuf);
#endif
}

SPDLOG_PISCEM_INLINE const char *spdlog_ex::what() const SPDLOG_PISCEM_NOEXCEPT
{
    return msg_.c_str();
}

SPDLOG_PISCEM_INLINE void throw_spdlog_ex(const std::string &msg, int last_errno)
{
    SPDLOG_PISCEM_THROW(spdlog_ex(msg, last_errno));
}

SPDLOG_PISCEM_INLINE void throw_spdlog_ex(std::string msg)
{
    SPDLOG_PISCEM_THROW(spdlog_ex(std::move(msg)));
}

} // namespace spdlog_piscem
