// Copyright(c) 2015-present, Gabi Melman & spdlog contributors.
// Distributed under the MIT License (http://opensource.org/licenses/MIT)

#pragma once

#include <spdlog_piscem/common.h>
#include <ctime> // std::time_t

namespace spdlog_piscem {
namespace details {
namespace os {

SPDLOG_PISCEM_API spdlog_piscem::log_clock::time_point now() SPDLOG_PISCEM_NOEXCEPT;

SPDLOG_PISCEM_API std::tm localtime(const std::time_t &time_tt) SPDLOG_PISCEM_NOEXCEPT;

SPDLOG_PISCEM_API std::tm localtime() SPDLOG_PISCEM_NOEXCEPT;

SPDLOG_PISCEM_API std::tm gmtime(const std::time_t &time_tt) SPDLOG_PISCEM_NOEXCEPT;

SPDLOG_PISCEM_API std::tm gmtime() SPDLOG_PISCEM_NOEXCEPT;

// eol definition
#if !defined(SPDLOG_PISCEM_EOL)
#    ifdef _WIN32
#        define SPDLOG_PISCEM_EOL "\r\n"
#    else
#        define SPDLOG_PISCEM_EOL "\n"
#    endif
#endif

SPDLOG_PISCEM_CONSTEXPR static const char *default_eol = SPDLOG_PISCEM_EOL;

// folder separator
#if !defined(SPDLOG_PISCEM_FOLDER_SEPS)
#    ifdef _WIN32
#        define SPDLOG_PISCEM_FOLDER_SEPS "\\/"
#    else
#        define SPDLOG_PISCEM_FOLDER_SEPS "/"
#    endif
#endif

SPDLOG_PISCEM_CONSTEXPR static const char folder_seps[] = SPDLOG_PISCEM_FOLDER_SEPS;
SPDLOG_PISCEM_CONSTEXPR static const filename_t::value_type folder_seps_filename[] = SPDLOG_PISCEM_FILENAME_T(SPDLOG_PISCEM_FOLDER_SEPS);

// fopen_s on non windows for writing
SPDLOG_PISCEM_API bool fopen_s(FILE **fp, const filename_t &filename, const filename_t &mode);

// Remove filename. return 0 on success
SPDLOG_PISCEM_API int remove(const filename_t &filename) SPDLOG_PISCEM_NOEXCEPT;

// Remove file if exists. return 0 on success
// Note: Non atomic (might return failure to delete if concurrently deleted by other process/thread)
SPDLOG_PISCEM_API int remove_if_exists(const filename_t &filename) SPDLOG_PISCEM_NOEXCEPT;

SPDLOG_PISCEM_API int rename(const filename_t &filename1, const filename_t &filename2) SPDLOG_PISCEM_NOEXCEPT;

// Return if file exists.
SPDLOG_PISCEM_API bool path_exists(const filename_t &filename) SPDLOG_PISCEM_NOEXCEPT;

// Return file size according to open FILE* object
SPDLOG_PISCEM_API size_t filesize(FILE *f);

// Return utc offset in minutes or throw spdlog_ex on failure
SPDLOG_PISCEM_API int utc_minutes_offset(const std::tm &tm = details::os::localtime());

// Return current thread id as size_t
// It exists because the std::this_thread::get_id() is much slower(especially
// under VS 2013)
SPDLOG_PISCEM_API size_t _thread_id() SPDLOG_PISCEM_NOEXCEPT;

// Return current thread id as size_t (from thread local storage)
SPDLOG_PISCEM_API size_t thread_id() SPDLOG_PISCEM_NOEXCEPT;

// This is avoid msvc issue in sleep_for that happens if the clock changes.
// See https://github.com/gabime/spdlog_piscem/issues/609
SPDLOG_PISCEM_API void sleep_for_millis(unsigned int milliseconds) SPDLOG_PISCEM_NOEXCEPT;

SPDLOG_PISCEM_API std::string filename_to_str(const filename_t &filename);

SPDLOG_PISCEM_API int pid() SPDLOG_PISCEM_NOEXCEPT;

// Determine if the terminal supports colors
// Source: https://github.com/agauniyal/rang/
SPDLOG_PISCEM_API bool is_color_terminal() SPDLOG_PISCEM_NOEXCEPT;

// Determine if the terminal attached
// Source: https://github.com/agauniyal/rang/
SPDLOG_PISCEM_API bool in_terminal(FILE *file) SPDLOG_PISCEM_NOEXCEPT;

#if (defined(SPDLOG_PISCEM_WCHAR_TO_UTF8_SUPPORT) || defined(SPDLOG_PISCEM_WCHAR_FILENAMES)) && defined(_WIN32)
SPDLOG_PISCEM_API void wstr_to_utf8buf(wstring_view_t wstr, memory_buf_t &target);

SPDLOG_PISCEM_API void utf8_to_wstrbuf(string_view_t str, wmemory_buf_t &target);
#endif

// Return directory name from given path or empty string
// "abc/file" => "abc"
// "abc/" => "abc"
// "abc" => ""
// "abc///" => "abc//"
SPDLOG_PISCEM_API filename_t dir_name(const filename_t &path);

// Create a dir from the given path.
// Return true if succeeded or if this dir already exists.
SPDLOG_PISCEM_API bool create_dir(const filename_t &path);

// non thread safe, cross platform getenv/getenv_s
// return empty string if field not found
SPDLOG_PISCEM_API std::string getenv(const char *field);

} // namespace os
} // namespace details
} // namespace spdlog_piscem

#ifdef SPDLOG_PISCEM_HEADER_ONLY
#    include "os-inl.h"
#endif
