cmake_minimum_required(VERSION 3.16)
project(aeolion CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
if(NOT MSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "-O2")
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

option(AEOLION_BUILD_TESTS "Build the test suite" ON)
option(AEOLION_WARNINGS "Build with warnings enabled" ON)

# ---------------------------------------------------------------- library --
# Header-only: everything lives under include/Aeolion/*.h. This INTERFACE
# target hands executables and tests the include path, the C++20 requirement,
# warnings, and the LAPACK backend via a single target_link_libraries(aeolion).
add_library(aeolion INTERFACE)
target_include_directories(aeolion INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(aeolion INTERFACE cxx_std_20)

include(CompilerWarnings)
aeolion_enable_warnings(aeolion)

include(LapackBackend)
aeolion_link_lapack(aeolion)
find_package(nlohmann_json CONFIG REQUIRED)
target_link_libraries(aeolion INTERFACE nlohmann_json::nlohmann_json)

# ---------------------------------------------------- executables & tests --
add_subdirectory(src)

if(AEOLION_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
