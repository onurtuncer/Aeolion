cmake_minimum_required(VERSION 3.16)
project(aeolion CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
if(NOT MSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "-O2")
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

option(AEOLION_BUILD_TESTS "Build the test suite" ON)
option(AEOLION_BUILD_VIEWER "Build the OpenGL viewer application" ON)
option(AEOLION_WARNINGS "Build with warnings enabled" ON)

# ---------------------------------------------------------------- library --
# Header-only: include/Aeolion/<Namespace>/*.h holds Math, Geometry, BEMT,
# and DragEstimate; Solver lives in its own top-level solver/include/
# (mirroring viewer/'s layout) since it's the one module other modules build
# on top of. This INTERFACE target hands executables and tests both include
# paths, the C++23 requirement, warnings, and the LAPACK backend via a
# single target_link_libraries(aeolion).
add_library(aeolion INTERFACE)
target_include_directories(aeolion INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/solver/include
)
target_compile_features(aeolion INTERFACE cxx_std_23)

include(CompilerWarnings)
aeolion_enable_warnings(aeolion)

include(LapackBackend)
aeolion_link_lapack(aeolion)
find_package(nlohmann_json CONFIG REQUIRED)
target_link_libraries(aeolion INTERFACE nlohmann_json::nlohmann_json)

# ---------------------------------------------------- executables & tests --
add_subdirectory(src)

if(AEOLION_BUILD_VIEWER)
    add_subdirectory(viewer)
endif()

if(AEOLION_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
