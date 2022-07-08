cmake_minimum_required(VERSION 3.5)
project(SSHASH)

# Set path for modules required to search for existing packages in the system.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

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
find_package(Threads REQUIRED)

# Prepare the `jemalloc` library. It provides scalable concurrency support and better avoidance
# of fragmentation.
set(MALLOC_LIB "")
set(JEMALLOC_MIN_VERSION "5.2.1")
find_package(Jemalloc ${JEMALLOC_MIN_VERSION})
if(JEMALLOC_FOUND)
    message("Found jemalloc (v${JEMALLOC_VERSION}) in the system")
    set(MALLOC_LIB ${JEMALLOC_LIBRARIES})
else()
  message("Could not find jemalloc (>= v${JEMALLOC_MIN_VERSION}). Please consider installing this.")
endif()

set(Z_LIB_SOURCES
  include/gz/zip_stream.cpp
)

set(SSHASH_SOURCES
  include/dictionary.cpp
  include/info.cpp
)

# Create a static lib
add_library(sshash_static STATIC
  ${Z_LIB_SOURCES}
  ${SSHASH_SOURCES}
  include/builder/build.cpp
)
target_include_directories(sshash_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(build src/build.cpp)
target_link_libraries(build sshash_static z ${MALLOC_LIB})
target_include_directories(build PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(bench src/bench.cpp ${SSHASH_SOURCES})

add_executable(check src/check.cpp ${SSHASH_SOURCES} ${Z_LIB_SOURCES})
target_link_libraries(check z)

add_executable(query src/query.cpp ${SSHASH_SOURCES} ${Z_LIB_SOURCES})
target_link_libraries(query z)

#add_executable(build_contig_table src/build_contig_table.cpp)
#target_include_directories(build_contig_table PUBLIC ${CMAKE_SOURCE_DIR}/include)
#target_link_libraries(build_contig_table z)

add_executable(test_load src/test_load.cpp)
target_include_directories(test_load PUBLIC ${CMAKE_SOURCE_DIR}/include ${ZLIB_INCLUDE_DIRS})
target_link_libraries(test_load ZLIB::ZLIB sshash_static)

add_executable(pesc src/pesc.cpp) 
target_include_directories(pesc PUBLIC ${CMAKE_SOURCE_DIR}/include ${ZLIB_INCLUDE_DIRS})
target_link_libraries(pesc ZLIB::ZLIB Threads::Threads sshash_static ${MALLOC_LIB}) 

add_executable(evaluator src/index_evaluator.cpp) 
target_include_directories(evaluator PUBLIC ${CMAKE_SOURCE_DIR}/include ${ZLIB_INCLUDE_DIRS})
target_link_libraries(evaluator ZLIB::ZLIB Threads::Threads sshash_static) 

add_executable(test_parse src/test_parse_geo.cpp)

