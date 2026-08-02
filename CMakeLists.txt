cmake_minimum_required(VERSION 3.14)

project(snrk LANGUAGES CXX)

option(BUILD_SHARED_LIBS OFF)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(PROJECT_SOURCES
    include/snrk.h
    include/circuit.h src/circuit.cpp
    include/setup.h src/setup.cpp
    include/polynom.h src/polynom.cpp
    include/proof.h src/proof.cpp
    include/types.h src/types.cpp
    include/threader.h src/threader.cpp
    include/params.h
)

add_library(snrk ${PROJECT_SOURCES})

target_include_directories(snrk PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

target_compile_options(snrk PRIVATE
    $<$<CONFIG:Release>:-O3 -march=native>
)

target_link_libraries(snrk PUBLIC mcl crypto)

install(TARGETS snrk
    LIBRARY DESTINATION lib
)

install(DIRECTORY include/ DESTINATION include/snrk)
