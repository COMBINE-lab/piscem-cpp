// Copyright(c) 2015-present, Gabi Melman & spdlog contributors.
// Distributed under the MIT License (http://opensource.org/licenses/MIT)

#pragma once

#ifndef SPDLOG_PISCEM_HEADER_ONLY
#    include <spdlog_piscem/details/periodic_worker.h>
#endif

namespace spdlog_piscem {
namespace details {

// stop the worker thread and join it
SPDLOG_PISCEM_INLINE periodic_worker::~periodic_worker()
{
    if (worker_thread_.joinable())
    {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            active_ = false;
        }
        cv_.notify_one();
        worker_thread_.join();
    }
}

} // namespace details
} // namespace spdlog_piscem
