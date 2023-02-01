// Copyright(c) 2015-present, Gabi Melman & spdlog contributors.
// Distributed under the MIT License (http://opensource.org/licenses/MIT)

#pragma once

#ifndef SPDLOG_PISCEM_HEADER_ONLY
#    include <spdlog_piscem/async_logger.h>
#endif

#include <spdlog_piscem/sinks/sink.h>
#include <spdlog_piscem/details/thread_pool.h>

#include <memory>
#include <string>

SPDLOG_PISCEM_INLINE spdlog_piscem::async_logger::async_logger(
    std::string logger_name, sinks_init_list sinks_list, std::weak_ptr<details::thread_pool> tp, async_overflow_policy overflow_policy)
    : async_logger(std::move(logger_name), sinks_list.begin(), sinks_list.end(), std::move(tp), overflow_policy)
{}

SPDLOG_PISCEM_INLINE spdlog_piscem::async_logger::async_logger(
    std::string logger_name, sink_ptr single_sink, std::weak_ptr<details::thread_pool> tp, async_overflow_policy overflow_policy)
    : async_logger(std::move(logger_name), {std::move(single_sink)}, std::move(tp), overflow_policy)
{}

// send the log message to the thread pool
SPDLOG_PISCEM_INLINE void spdlog_piscem::async_logger::sink_it_(const details::log_msg &msg)
{
    if (auto pool_ptr = thread_pool_.lock())
    {
        pool_ptr->post_log(shared_from_this(), msg, overflow_policy_);
    }
    else
    {
        throw_spdlog_ex("async log: thread pool doesn't exist anymore");
    }
}

// send flush request to the thread pool
SPDLOG_PISCEM_INLINE void spdlog_piscem::async_logger::flush_()
{
    if (auto pool_ptr = thread_pool_.lock())
    {
        pool_ptr->post_flush(shared_from_this(), overflow_policy_);
    }
    else
    {
        throw_spdlog_ex("async flush: thread pool doesn't exist anymore");
    }
}

//
// backend functions - called from the thread pool to do the actual job
//
SPDLOG_PISCEM_INLINE void spdlog_piscem::async_logger::backend_sink_it_(const details::log_msg &msg)
{
    for (auto &sink : sinks_)
    {
        if (sink->should_log(msg.level))
        {
            SPDLOG_PISCEM_TRY
            {
                sink->log(msg);
            }
            SPDLOG_PISCEM_LOGGER_CATCH(msg.source)
        }
    }

    if (should_flush_(msg))
    {
        backend_flush_();
    }
}

SPDLOG_PISCEM_INLINE void spdlog_piscem::async_logger::backend_flush_()
{
    for (auto &sink : sinks_)
    {
        SPDLOG_PISCEM_TRY
        {
            sink->flush();
        }
        SPDLOG_PISCEM_LOGGER_CATCH(source_loc())
    }
}

SPDLOG_PISCEM_INLINE std::shared_ptr<spdlog_piscem::logger> spdlog_piscem::async_logger::clone(std::string new_name)
{
    auto cloned = std::make_shared<spdlog_piscem::async_logger>(*this);
    cloned->name_ = std::move(new_name);
    return cloned;
}
