
set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)


set(SRC
        Header.cpp
        Tags.cpp
        Buffer.cpp
        Byte_Array.cpp
        RAD_Writer.cpp
)


add_library(libradicl STATIC ${SRC})


target_include_directories(libradicl PRIVATE ${INCLUDE_DIR})


set(COMPILE_FLAGS ${WARNING_FLAGS} ${OPTIMIZE_FLAGS})


target_compile_options(libradicl PRIVATE ${COMPILE_FLAGS})


install(TARGETS libradicl
        DESTINATION lib)


add_executable(test test.cpp)

target_include_directories(test PRIVATE ${INCLUDE_DIR})
target_compile_options(test PRIVATE ${COMPILE_FLAGS})
target_compile_options(test PRIVATE ${OPTIMIZE_FLAGS})
target_link_libraries(test libradicl)
