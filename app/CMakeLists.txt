# Driver executables. Each links the header-only `aeolion` interface library,
# which carries the include path, C++23 requirement, warnings, and LAPACK.

add_executable(solver_demo main.cpp)
target_link_libraries(solver_demo PRIVATE aeolion)

add_executable(aeolion_geometry GeometryContractCLI.cpp)
target_link_libraries(aeolion_geometry PRIVATE aeolion)

# Also needs the compiled BEMT solver, not just the header-only core.
add_executable(aeolion_prop PropellerCLI.cpp)
target_link_libraries(aeolion_prop PRIVATE aeolion aeolion_bemt)
