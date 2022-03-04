cmake_minimum_required(VERSION 3.5)
project(SSHASH)

set(CMAKE_CXX_STANDARD 17)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif ()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif ()

MESSAGE(STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if (UNIX)

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-missing-braces")

  # Flags for PTHash:
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mbmi2 -msse4.2") # for hardware popcount and pdep
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")        # for multi-threading

  if (SSHASH_USE_SANITIZERS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
  endif()

endif()

find_package(ZLIB REQUIRED)

add_executable(build src/build.cpp)
target_link_libraries(build z)

add_executable(bench src/bench.cpp)

add_executable(check src/check.cpp)
target_link_libraries(check z)

add_executable(query src/query.cpp)
target_link_libraries(query z)

add_executable(build_contig_table src/build_contig_table.cpp)
target_include_directories(build_contig_table PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(build_contig_table z)

add_executable(test_load src/test_load.cpp)
target_include_directories(test_load PUBLIC ${CMAKE_SOURCE_DIR}/include ${ZLIB_INCLUDE_DIRS})
target_link_libraries(test_load ZLIB::ZLIB)

add_executable(mapper src/mapper.cpp) 
target_include_directories(mapper PUBLIC ${CMAKE_SOURCE_DIR}/include ${ZLIB_INCLUDE_DIRS})
target_link_libraries(mapper ZLIB::ZLIB)
