# 3.15 introuced many nice features and seems 
# a reasonable minimum version worth supporting.
cmake_minimum_required(VERSION 3.15)
project(SSHASH)
include_directories(AFTER external/sshash/src external/sshash/include external/sshash external/sshash/external/pthash external/sshash/external/pthash/external/essentials/include)
add_subdirectory(${CMAKE_SOURCE_DIR}/external/libradicl)
# Set path for modules required to search for existing packages in the system.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
add_definitions(-DSSHASH_USE_TRADITIONAL_NUCLEOTIDE_ENCODING)

set(CMAKE_CXX_STANDARD 17)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif ()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif ()

MESSAGE(STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE})

# if CONDA_BUILD is defined AND true, then set a local variable for this
if (CONDA_BUILD)
  # it's already true, we don't have to do anything
else()
  if ($ENV{CONDA_BUILD})
    # if it wasn't passed in but is an environment variable
    set(CONDA_BUILD TRUE)
  endif()
endif()

if (CONDA_BUILD)
  message(STATUS "Building in the context of CONDA_BUILD = TRUE")
else()
  message(STATUS "Building in the context of CONDA_BUILD = FALSE")
endif()

set(CUSTOM_INSTRUCTION_FLAGS "")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
MESSAGE(STATUS "Compiling for processor: " ${CMAKE_HOST_SYSTEM_PROCESSOR})
if (UNIX AND (CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "x86_64"))
  
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.2") # for vectorized instructions
  list(APPEND CUSTOM_INSTRUCTION_FLAGS "-msse4.2")

  if (NOT CONDA_BUILD) 
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
    list(APPEND CUSTOM_INSTRUCTION_FLAGS "-march=native")
  endif()

  if(NOT NO_BMI2)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mbmi2") # for hardware popcount and pdep
    list(APPEND CUSTOM_INSTRUCTION_FLAGS "-mbmi2")
  endif()
  MESSAGE(STATUS "Compiling with custom instruction flags: ${CUSTOM_INSTRUCTION_FLAGS}")
endif()
## figure out what we might want / need to put here if we have a non x86_64 arch.

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DL1_CACHE_LINE_SIZE=64")
if (UNIX)

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-missing-braces")

  # Flags for PTHash:
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")        # for multi-threading

  if (SSHASH_USE_SANITIZERS)
    message(STATUS "Using sanitizers. Compiling with flags: -fsanitize=address -fno-omit-frame-pointer")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
  endif()

  if (CONDA_BUILD)
    # because of the kind of ridiculous RAM limitations on the conda build jobs
    add_compile_options("SHELL:--param ggc-min-expand=20" "SHELL:--param ggc-min-heapsize=8192")
  endif()

endif()

#add_subdirectory(external/zlib-cloudflare)
include(ExternalProject)
ExternalProject_Add(zlib-cloudflare
  GIT_REPOSITORY https://github.com/cloudflare/zlib
  GIT_TAG 7aa510344e06fecd6fe09195ac22e9a424ceb660
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/zlib-cloudflare
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/zlib-cloudflare
)

ExternalProject_Add_Step(zlib-cloudflare
  rename_lib
  COMMAND cp ${CMAKE_BINARY_DIR}/zlib-cloudflare/lib/libz.a ${CMAKE_BINARY_DIR}/zlib-cloudflare/lib/libzcf.a
  COMMENT "copying cloudflare libz.a to libzcf.a"
  DEPENDEES install
)

## can't reply on find_package because the ExternalProject_Add
# steps won't complete by compile time, when the find_package 
# looks for the ZLIB
#set(ZLIB_ROOT ${CMAKE_BINARY_DIR}/zlib-cloudflare)
#find_package(ZLIB REQUIRED)


# install path to the lib / include directories of the 
# cloudflare libz
include_directories(${CMAKE_BINARY_DIR}/zlib-cloudflare/include)
link_directories(${CMAKE_BINARY_DIR}/zlib-cloudflare/lib)

# direct path to the cloudflare libz
set(CLOUDFLARE_ZLIB ${CMAKE_BINARY_DIR}/zlib-cloudflare/lib/libzcf.a)


set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED) # for pthreads

#find_package(GPerfTools REQUIRED) ## gperftools
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
  external/sshash/include/gz/zip_stream.cpp
)

set(SSHASH_SOURCES
  external/sshash/include/dictionary.cpp
  external/sshash/include/info.cpp
)

set(LIBRADICL_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/external/libradicl/include)

# Create a static lib
add_library(sshash_static STATIC
  ${Z_LIB_SOURCES}
  ${SSHASH_SOURCES}
  external/sshash/include/builder/build.cpp
)
target_include_directories(sshash_static PRIVATE  ${LIBRADICL_INCLUDE})


add_library(pesc_static STATIC
  src/pesc_bulk.cpp src/pesc_sc.cpp src/pesc_sc_atac.cpp src/hit_searcher.cpp src/FastxParser.cpp
)
target_include_directories(pesc_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include ${ZLIB_INCLUDE_DIRS} ${LIBRADICL_INCLUDE})

add_library(build_static STATIC
  src/build.cpp
  src/poison_table_builder.cpp)
target_include_directories(build_static PRIVATE  ${ZLIB_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/include)

#add_executable(build src/build.cpp)
#target_link_libraries(build sshash_static z ${MALLOC_LIB})
#target_include_directories(build PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include ${LIBRADICL_INCLUDE})

#add_executable(bench src/bench.cpp ${SSHASH_SOURCES})
#target_include_directories(bench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include ${LIBRADICL_INCLUDE})

add_executable(check src/check.cpp src/FastxParser.cpp)
target_link_libraries(check ${CLOUDFLARE_ZLIB} sshash_static Threads::Threads)
target_include_directories(check PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include ${ZLIB_INCLUDE_DIRS} ${LIBRADICL_INCLUDE})

#add_executable(query src/query.cpp ${SSHASH_SOURCES} ${Z_LIB_SOURCES})
#target_link_libraries(query z Threads::Threads)
#target_include_directories(query PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include ${LIBRADICL_INCLUDE})

#add_executable(test_load src/test_load.cpp)
#target_include_directories(test_load PUBLIC ${CMAKE_SOURCE_DIR}/include ${ZLIB_INCLUDE_DIRS} ${LIBRADICL_INCLUDE})
#target_link_libraries(test_load ${CLOUDFLARE_ZLIB} sshash_static Threads::Threads)

add_executable(pesc-sc src/pesc_sc_runner.cpp) 
target_include_directories(pesc-sc PUBLIC ${CMAKE_SOURCE_DIR}/include ${ZLIB_INCLUDE_DIRS} ${LIBRADICL_INCLUDE})
target_link_libraries(pesc-sc pesc_static ${CLOUDFLARE_ZLIB} Threads::Threads sshash_static ${MALLOC_LIB})

add_executable(pesc-bulk src/pesc_bulk_runner.cpp) 
target_include_directories(pesc-bulk PUBLIC ${CMAKE_SOURCE_DIR}/include ${ZLIB_INCLUDE_DIRS} ${LIBRADICL_INCLUDE})
target_link_libraries(pesc-bulk pesc_static ${CLOUDFLARE_ZLIB} Threads::Threads sshash_static ${MALLOC_LIB})

add_executable(pesc-sc-atac src/pesc_sc_atac_runner.cpp) 
target_include_directories(pesc-sc-atac PUBLIC ${LIBRADICL_INCLUDE} ${CMAKE_SOURCE_DIR}/include ${ZLIB_INCLUDE_DIRS})
target_link_libraries(pesc-sc-atac pesc_static radicl ${CLOUDFLARE_ZLIB} Threads::Threads sshash_static ${MALLOC_LIB})
#target_link_directories(pesc-sc-atac PUBLIC /fs/cbcb-lab/rob/students/noor/miniforge3/lib)
add_dependencies(pesc-sc-atac radicl)

# add_executable(check-overlap src/check_overlap.cpp) 
# target_include_directories(check-overlap PUBLIC ${CMAKE_SOURCE_DIR}/include ${ZLIB_INCLUDE_DIRS})
# target_link_libraries(check-overlap pesc_static ZLIB::ZLIB Threads::Threads sshash_static ${MALLOC_LIB})



#add_executable(mapper src/mapper.cpp src/hit_searcher.cpp) 
#target_include_directories(mapper PUBLIC ${CMAKE_SOURCE_DIR}/include ${ZLIB_INCLUDE_DIRS})
#target_link_libraries(mapper ${CLOUDFLARE_ZLIB} Threads::Threads sshash_static ${MALLOC_LIB}) 

add_executable(build src/build_runner.cpp) 
target_include_directories(build PUBLIC ${CMAKE_SOURCE_DIR}/include ${ZLIB_INCLUDE_DIRS} ${LIBRADICL_INCLUDE})
target_link_libraries(build Threads::Threads build_static sshash_static ${CLOUDFLARE_ZLIB} ${MALLOC_LIB}) 

add_executable(evaluator src/index_evaluator.cpp) 
target_include_directories(evaluator PUBLIC ${CMAKE_SOURCE_DIR}/include ${ZLIB_INCLUDE_DIRS} ${LIBRADICL_INCLUDE})
target_link_libraries(evaluator ${CLOUDFLARE_ZLIB} Threads::Threads sshash_static) 

add_executable(build-poison-table src/poison_table_build_runner.cpp) 
target_include_directories(build-poison-table PUBLIC ${CMAKE_SOURCE_DIR}/include ${ZLIB_INCLUDE_DIRS})
target_link_libraries(build-poison-table Threads::Threads build_static sshash_static pesc_static ${CLOUDFLARE_ZLIB} ${MALLOC_LIB})

add_executable(poison_filter src/poison_read_filter.cpp) 
target_include_directories(poison_filter PUBLIC ${CMAKE_SOURCE_DIR}/include ${ZLIB_INCLUDE_DIRS})
target_link_libraries(poison_filter Threads::Threads sshash_static pesc_static ${CLOUDFLARE_ZLIB}) 

## tests
add_subdirectory(external/doctest)

# Make test executable
add_executable(tests tests/test_geom_extractor.cpp)
target_compile_features(tests PRIVATE cxx_std_17)
target_include_directories(tests PRIVATE doctest)
target_link_libraries(tests PRIVATE doctest pesc_static)


## depend on the cloudflare zlib library
add_dependencies(check zlib-cloudflare)
add_dependencies(evaluator zlib-cloudflare)
add_dependencies(build zlib-cloudflare)
add_dependencies(pesc-sc zlib-cloudflare)
add_dependencies(pesc-sc-atac zlib-cloudflare)
add_dependencies(pesc-bulk zlib-cloudflare)
add_dependencies(build-poison-table zlib-cloudflare)
add_dependencies(poison_filter zlib-cloudflare)


#add_executable(test_parse src/test_parse_geo.cpp)
#target_include_directories(evaluator PRIVATE ${CMAKE_SOURCE_DIR}/include)

install(TARGETS build build-poison-table pesc-sc pesc-bulk pesc-sc-atac
  CONFIGURATIONS Release
  RUNTIME DESTINATION Release/bin)

install(TARGETS build build-poison-table pesc-sc pesc-bulk pesc-sc-atac
        CONFIGURATIONS Debug
        RUNTIME DESTINATION Debug/bin)

install(TARGETS pesc_static build_static sshash_static
  CONFIGURATIONS Debug 
  RUNTIME DESTINATION Debug/lib)

install(FILES ${CLOUDFLARE_ZLIB}
  CONFIGURATIONS Debug
  RUNTIME DESTINATION Debug/lib)

install(TARGETS pesc_static build_static sshash_static
  CONFIGURATIONS Release
  RUNTIME DESTINATION Release/lib)

install(FILES ${CLOUDFLARE_ZLIB}
  CONFIGURATIONS Release
  RUNTIME DESTINATION Release/lib)
