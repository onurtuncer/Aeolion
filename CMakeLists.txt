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

# enable_testing() and aeolion_add_test() need to be in place before any
# add_subdirectory() below, since solver/ and panelbuilder/ each register
# their OWN tests/ subdirectory rather than a single central one -- every
# module is expected to own its build target and its tests together.
if(AEOLION_BUILD_TESTS)
  enable_testing()
  include(AeolionTest)
endif()

# ------------------------------------------------------------ core library --
# Header-only: include/Aeolion/<Namespace>/*.h holds Math, Lattice, Geometry,
# and DragEstimate. Solver lives in its own top-level solver/include/
# folder instead, because it is the module every other one builds on top
# of. This INTERFACE target hands consumers every header-only include
# path, the C++23 requirement, warnings, and the LAPACK backend.
add_library(aeolion INTERFACE)
target_include_directories(aeolion INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include
                                             ${CMAKE_CURRENT_SOURCE_DIR}/solver/include)
target_compile_features(aeolion INTERFACE cxx_std_23)

include(CompilerWarnings)
aeolion_enable_warnings(aeolion)

include(LapackBackend)
aeolion_link_lapack(aeolion)
find_package(nlohmann_json CONFIG REQUIRED)
target_link_libraries(aeolion INTERFACE nlohmann_json::nlohmann_json)

# The propeller BEMT solver was excised into its own project
# (https://github.com/onurtuncer/BEMT) -- it is a momentum method, not a
# panel method, and Aeolion carries no dependency on it. The handoff
# contract's propulsion_bemt block (Geometry::PropulsionSpec) stays as
# schema vocabulary for whatever propulsion method consumes it next.

# --------------------------------------------------- header-only modules --
# Owns nothing but its own tests/ subdirectory -- see solver/CMakeLists.txt.
add_subdirectory(solver)

# ------------------------------------------------------ compiled libraries --
add_subdirectory(logger)
add_subdirectory(panelbuilder)

# ---------------------------------------------------- executables & tests --
add_subdirectory(app)

if(AEOLION_BUILD_VIEWER)
  add_subdirectory(viewer)
endif()

if(AEOLION_BUILD_TESTS)
  # The core header-only library's (Math/Lattice/Geometry/DragEstimate)
  # own tests -- everything else lives beside the module it tests.
  add_subdirectory(tests)
endif()
