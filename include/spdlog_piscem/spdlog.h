// Copyright(c) 2015-present, Gabi Melman & spdlog contributors.
// Distributed under the MIT License (http://opensource.org/licenses/MIT)

// spdlog main header file.
// see example.cpp for usage example

#ifndef SPDLOG_PISCEM_H
#define SPDLOG_PISCEM_H

#pragma once

#include <spdlog_piscem/common.h>
#include <spdlog_piscem/details/registry.h>
#include <spdlog_piscem/logger.h>
#include <spdlog_piscem/version.h>
#include <spdlog_piscem/details/synchronous_factory.h>

#include <chrono>
#include <functional>
#include <memory>
#include <string>

namespace spdlog_piscem {

using default_factory = synchronous_factory;

// Create and register a logger with a templated sink type
// The logger's level, formatter and flush level will be set according the
// global settings.
//
// Example:
//   spdlog_piscem::create<daily_file_sink_st>("logger_name", "dailylog_filename", 11, 59);
template<typename Sink, typename... SinkArgs>
inline std::shared_ptr<spdlog_piscem::logger> create(std::string logger_name, SinkArgs &&... sink_args)
{
    return default_factory::create<Sink>(std::move(logger_name), std::forward<SinkArgs>(sink_args)...);
}

// Initialize and register a logger,
// formatter and flush level will be set according the global settings.
//
// Useful for initializing manually created loggers with the global settings.
//
// Example:
//   auto mylogger = std::make_shared<spdlog_piscem::logger>("mylogger", ...);
//   spdlog_piscem::initialize_logger(mylogger);
SPDLOG_PISCEM_API void initialize_logger(std::shared_ptr<logger> logger);

// Return an existing logger or nullptr if a logger with such name doesn't
// exist.
// example: spdlog_piscem::get("my_logger")->info("hello {}", "world");
SPDLOG_PISCEM_API std::shared_ptr<logger> get(const std::string &name);

// Set global formatter. Each sink in each logger will get a clone of this object
SPDLOG_PISCEM_API void set_formatter(std::unique_ptr<spdlog_piscem::formatter> formatter);

// Set global format string.
// example: spdlog_piscem::set_pattern("%Y-%m-%d %H:%M:%S.%e %l : %v");
SPDLOG_PISCEM_API void set_pattern(std::string pattern, pattern_time_type time_type = pattern_time_type::local);

// enable global backtrace support
SPDLOG_PISCEM_API void enable_backtrace(size_t n_messages);

// disable global backtrace support
SPDLOG_PISCEM_API void disable_backtrace();

// call dump backtrace on default logger
SPDLOG_PISCEM_API void dump_backtrace();

// Get global logging level
SPDLOG_PISCEM_API level::level_enum get_level();

// Set global logging level
SPDLOG_PISCEM_API void set_level(level::level_enum log_level);

// Determine whether the default logger should log messages with a certain level
SPDLOG_PISCEM_API bool should_log(level::level_enum lvl);

// Set global flush level
SPDLOG_PISCEM_API void flush_on(level::level_enum log_level);

// Start/Restart a periodic flusher thread
// Warning: Use only if all your loggers are thread safe!
template<typename Rep, typename Period>
inline void flush_every(std::chrono::duration<Rep, Period> interval)
{
    details::registry::instance().flush_every(interval);
}

// Set global error handler
SPDLOG_PISCEM_API void set_error_handler(void (*handler)(const std::string &msg));

// Register the given logger with the given name
SPDLOG_PISCEM_API void register_logger(std::shared_ptr<logger> logger);

// Apply a user defined function on all registered loggers
// Example:
// spdlog_piscem::apply_all([&](std::shared_ptr<spdlog_piscem::logger> l) {l->flush();});
SPDLOG_PISCEM_API void apply_all(const std::function<void(std::shared_ptr<logger>)> &fun);

// Drop the reference to the given logger
SPDLOG_PISCEM_API void drop(const std::string &name);

// Drop all references from the registry
SPDLOG_PISCEM_API void drop_all();

// stop any running threads started by spdlog and clean registry loggers
SPDLOG_PISCEM_API void shutdown();

// Automatic registration of loggers when using spdlog_piscem::create() or spdlog_piscem::create_async
SPDLOG_PISCEM_API void set_automatic_registration(bool automatic_registration);

// API for using default logger (stdout_color_mt),
// e.g: spdlog_piscem::info("Message {}", 1);
//
// The default logger object can be accessed using the spdlog_piscem::default_logger():
// For example, to add another sink to it:
// spdlog_piscem::default_logger()->sinks().push_back(some_sink);
//
// The default logger can replaced using spdlog_piscem::set_default_logger(new_logger).
// For example, to replace it with a file logger.
//
// IMPORTANT:
// The default API is thread safe (for _mt loggers), but:
// set_default_logger() *should not* be used concurrently with the default API.
// e.g do not call set_default_logger() from one thread while calling spdlog_piscem::info() from another.

SPDLOG_PISCEM_API std::shared_ptr<spdlog_piscem::logger> default_logger();

SPDLOG_PISCEM_API spdlog_piscem::logger *default_logger_raw();

SPDLOG_PISCEM_API void set_default_logger(std::shared_ptr<spdlog_piscem::logger> default_logger);

template<typename... Args>
inline void log(source_loc source, level::level_enum lvl, format_string_t<Args...> fmt, Args &&... args)
{
    default_logger_raw()->log(source, lvl, fmt, std::forward<Args>(args)...);
}

template<typename... Args>
inline void log(level::level_enum lvl, format_string_t<Args...> fmt, Args &&... args)
{
    default_logger_raw()->log(source_loc{}, lvl, fmt, std::forward<Args>(args)...);
}

template<typename... Args>
inline void trace(format_string_t<Args...> fmt, Args &&... args)
{
    default_logger_raw()->trace(fmt, std::forward<Args>(args)...);
}

template<typename... Args>
inline void debug(format_string_t<Args...> fmt, Args &&... args)
{
    default_logger_raw()->debug(fmt, std::forward<Args>(args)...);
}

template<typename... Args>
inline void info(format_string_t<Args...> fmt, Args &&... args)
{
    default_logger_raw()->info(fmt, std::forward<Args>(args)...);
}

template<typename... Args>
inline void warn(format_string_t<Args...> fmt, Args &&... args)
{
    default_logger_raw()->warn(fmt, std::forward<Args>(args)...);
}

template<typename... Args>
inline void error(format_string_t<Args...> fmt, Args &&... args)
{
    default_logger_raw()->error(fmt, std::forward<Args>(args)...);
}

template<typename... Args>
inline void critical(format_string_t<Args...> fmt, Args &&... args)
{
    default_logger_raw()->critical(fmt, std::forward<Args>(args)...);
}

template<typename T>
inline void log(source_loc source, level::level_enum lvl, const T &msg)
{
    default_logger_raw()->log(source, lvl, msg);
}

template<typename T>
inline void log(level::level_enum lvl, const T &msg)
{
    default_logger_raw()->log(lvl, msg);
}

#ifdef SPDLOG_PISCEM_WCHAR_TO_UTF8_SUPPORT
template<typename... Args>
inline void log(source_loc source, level::level_enum lvl, wformat_string_t<Args...> fmt, Args &&... args)
{
    default_logger_raw()->log(source, lvl, fmt, std::forward<Args>(args)...);
}

template<typename... Args>
inline void log(level::level_enum lvl, wformat_string_t<Args...> fmt, Args &&... args)
{
    default_logger_raw()->log(source_loc{}, lvl, fmt, std::forward<Args>(args)...);
}

template<typename... Args>
inline void trace(wformat_string_t<Args...> fmt, Args &&... args)
{
    default_logger_raw()->trace(fmt, std::forward<Args>(args)...);
}

template<typename... Args>
inline void debug(wformat_string_t<Args...> fmt, Args &&... args)
{
    default_logger_raw()->debug(fmt, std::forward<Args>(args)...);
}

template<typename... Args>
inline void info(wformat_string_t<Args...> fmt, Args &&... args)
{
    default_logger_raw()->info(fmt, std::forward<Args>(args)...);
}

template<typename... Args>
inline void warn(wformat_string_t<Args...> fmt, Args &&... args)
{
    default_logger_raw()->warn(fmt, std::forward<Args>(args)...);
}

template<typename... Args>
inline void error(wformat_string_t<Args...> fmt, Args &&... args)
{
    default_logger_raw()->error(fmt, std::forward<Args>(args)...);
}

template<typename... Args>
inline void critical(wformat_string_t<Args...> fmt, Args &&... args)
{
    default_logger_raw()->critical(fmt, std::forward<Args>(args)...);
}
#endif

template<typename T>
inline void trace(const T &msg)
{
    default_logger_raw()->trace(msg);
}

template<typename T>
inline void debug(const T &msg)
{
    default_logger_raw()->debug(msg);
}

template<typename T>
inline void info(const T &msg)
{
    default_logger_raw()->info(msg);
}

template<typename T>
inline void warn(const T &msg)
{
    default_logger_raw()->warn(msg);
}

template<typename T>
inline void error(const T &msg)
{
    default_logger_raw()->error(msg);
}

template<typename T>
inline void critical(const T &msg)
{
    default_logger_raw()->critical(msg);
}

} // namespace spdlog_piscem

//
// enable/disable log calls at compile time according to global level.
//
// define SPDLOG_PISCEM_ACTIVE_LEVEL to one of those (before including spdlog.h):
// SPDLOG_PISCEM_LEVEL_TRACE,
// SPDLOG_PISCEM_LEVEL_DEBUG,
// SPDLOG_PISCEM_LEVEL_INFO,
// SPDLOG_PISCEM_LEVEL_WARN,
// SPDLOG_PISCEM_LEVEL_ERROR,
// SPDLOG_PISCEM_LEVEL_CRITICAL,
// SPDLOG_PISCEM_LEVEL_OFF
//

#ifndef SPDLOG_PISCEM_NO_SOURCE_LOC
#    define SPDLOG_PISCEM_LOGGER_CALL(logger, level, ...)                                                                                         \
        (logger)->log(spdlog_piscem::source_loc{__FILE__, __LINE__, SPDLOG_PISCEM_FUNCTION}, level, __VA_ARGS__)
#else
#    define SPDLOG_PISCEM_LOGGER_CALL(logger, level, ...) (logger)->log(spdlog_piscem::source_loc{}, level, __VA_ARGS__)
#endif

#if SPDLOG_PISCEM_ACTIVE_LEVEL <= SPDLOG_PISCEM_LEVEL_TRACE
#    define SPDLOG_PISCEM_LOGGER_TRACE(logger, ...) SPDLOG_PISCEM_LOGGER_CALL(logger, spdlog_piscem::level::trace, __VA_ARGS__)
#    define SPDLOG_PISCEM_TRACE(...) SPDLOG_PISCEM_LOGGER_TRACE(spdlog_piscem::default_logger_raw(), __VA_ARGS__)
#else
#    define SPDLOG_PISCEM_LOGGER_TRACE(logger, ...) (void)0
#    define SPDLOG_PISCEM_TRACE(...) (void)0
#endif

#if SPDLOG_PISCEM_ACTIVE_LEVEL <= SPDLOG_PISCEM_LEVEL_DEBUG
#    define SPDLOG_PISCEM_LOGGER_DEBUG(logger, ...) SPDLOG_PISCEM_LOGGER_CALL(logger, spdlog_piscem::level::debug, __VA_ARGS__)
#    define SPDLOG_PISCEM_DEBUG(...) SPDLOG_PISCEM_LOGGER_DEBUG(spdlog_piscem::default_logger_raw(), __VA_ARGS__)
#else
#    define SPDLOG_PISCEM_LOGGER_DEBUG(logger, ...) (void)0
#    define SPDLOG_PISCEM_DEBUG(...) (void)0
#endif

#if SPDLOG_PISCEM_ACTIVE_LEVEL <= SPDLOG_PISCEM_LEVEL_INFO
#    define SPDLOG_PISCEM_LOGGER_INFO(logger, ...) SPDLOG_PISCEM_LOGGER_CALL(logger, spdlog_piscem::level::info, __VA_ARGS__)
#    define SPDLOG_PISCEM_INFO(...) SPDLOG_PISCEM_LOGGER_INFO(spdlog_piscem::default_logger_raw(), __VA_ARGS__)
#else
#    define SPDLOG_PISCEM_LOGGER_INFO(logger, ...) (void)0
#    define SPDLOG_PISCEM_INFO(...) (void)0
#endif

#if SPDLOG_PISCEM_ACTIVE_LEVEL <= SPDLOG_PISCEM_LEVEL_WARN
#    define SPDLOG_PISCEM_LOGGER_WARN(logger, ...) SPDLOG_PISCEM_LOGGER_CALL(logger, spdlog_piscem::level::warn, __VA_ARGS__)
#    define SPDLOG_PISCEM_WARN(...) SPDLOG_PISCEM_LOGGER_WARN(spdlog_piscem::default_logger_raw(), __VA_ARGS__)
#else
#    define SPDLOG_PISCEM_LOGGER_WARN(logger, ...) (void)0
#    define SPDLOG_PISCEM_WARN(...) (void)0
#endif

#if SPDLOG_PISCEM_ACTIVE_LEVEL <= SPDLOG_PISCEM_LEVEL_ERROR
#    define SPDLOG_PISCEM_LOGGER_ERROR(logger, ...) SPDLOG_PISCEM_LOGGER_CALL(logger, spdlog_piscem::level::err, __VA_ARGS__)
#    define SPDLOG_PISCEM_ERROR(...) SPDLOG_PISCEM_LOGGER_ERROR(spdlog_piscem::default_logger_raw(), __VA_ARGS__)
#else
#    define SPDLOG_PISCEM_LOGGER_ERROR(logger, ...) (void)0
#    define SPDLOG_PISCEM_ERROR(...) (void)0
#endif

#if SPDLOG_PISCEM_ACTIVE_LEVEL <= SPDLOG_PISCEM_LEVEL_CRITICAL
#    define SPDLOG_PISCEM_LOGGER_CRITICAL(logger, ...) SPDLOG_PISCEM_LOGGER_CALL(logger, spdlog_piscem::level::critical, __VA_ARGS__)
#    define SPDLOG_PISCEM_CRITICAL(...) SPDLOG_PISCEM_LOGGER_CRITICAL(spdlog_piscem::default_logger_raw(), __VA_ARGS__)
#else
#    define SPDLOG_PISCEM_LOGGER_CRITICAL(logger, ...) (void)0
#    define SPDLOG_PISCEM_CRITICAL(...) (void)0
#endif

#ifdef SPDLOG_PISCEM_HEADER_ONLY
#    include "spdlog-inl.h"
#endif

#endif // SPDLOG_PISCEM_H
