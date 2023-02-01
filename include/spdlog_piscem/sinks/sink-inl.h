// Copyright(c) 2015-present, Gabi Melman & spdlog contributors.
// Distributed under the MIT License (http://opensource.org/licenses/MIT)

#pragma once

#ifndef SPDLOG_PISCEM_HEADER_ONLY
#    include <spdlog_piscem/sinks/sink.h>
#endif

#include <spdlog_piscem/common.h>

SPDLOG_PISCEM_INLINE bool spdlog_piscem::sinks::sink::should_log(spdlog_piscem::level::level_enum msg_level) const
{
    return msg_level >= level_.load(std::memory_order_relaxed);
}

SPDLOG_PISCEM_INLINE void spdlog_piscem::sinks::sink::set_level(level::level_enum log_level)
{
    level_.store(log_level, std::memory_order_relaxed);
}

SPDLOG_PISCEM_INLINE spdlog_piscem::level::level_enum spdlog_piscem::sinks::sink::level() const
{
    return static_cast<spdlog_piscem::level::level_enum>(level_.load(std::memory_order_relaxed));
}
