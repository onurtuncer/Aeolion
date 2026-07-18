cmake_minimum_required(VERSION 3.16)
project(aeolion CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
if(NOT MSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "-O2")
endif()

option(AEOLION_BUILD_TESTS "Build the test suite" ON)
option(AEOLION_WARNINGS "Build with -Wall -Wextra" ON)

# ---------------------------------------------------------------- library --
# Header-only: everything lives under include/Aeolion/*.h. This target
# exists so executables (and tests) just do target_link_libraries(... aeolion)
# and get the include path and the LAPACK backend for free.
add_library(aeolion INTERFACE)
target_include_directories(aeolion INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(aeolion INTERFACE cxx_std_17)

if(AEOLION_WARNINGS)
    if(MSVC)
        target_compile_options(aeolion INTERFACE /W4)
    else()
        target_compile_options(aeolion INTERFACE -Wall -Wextra)
    endif()
endif()

# The dense linear solve calls LAPACK's dgetrf/dgetrs through their Fortran
# ABI directly -- no LAPACKE C header needed, so a LAPACK *library* is the
# only hard dependency. CMake's FindLAPACK locates OpenBLAS, reference
# LAPACK, MKL, Accelerate, etc. and pulls in the matching BLAS.
find_package(LAPACK REQUIRED)
target_link_libraries(aeolion INTERFACE LAPACK::LAPACK)
message(STATUS "aeolion: using LAPACK (${LAPACK_LIBRARIES})")

# -------------------------------------------------------------- executables
add_executable(vlm_demo src/Main.cpp)
target_link_libraries(vlm_demo PRIVATE aeolion)

add_executable(geometry_vlm src/GeometryVlm.cpp)
target_link_libraries(geometry_vlm PRIVATE aeolion)

add_executable(airfoil_polar src/AirfoilPolar.cpp)
target_link_libraries(airfoil_polar PRIVATE aeolion)

# ------------------------------------------------------------------- tests
if(AEOLION_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
