# Driver executables. Each links the header-only `aeolion` interface library,
# which carries the include path, C++23 requirement, warnings, and LAPACK,
# plus the logger every application reports through.

add_executable(solver_demo main.cpp)
target_link_libraries(solver_demo PRIVATE aeolion aeolion_logger)

add_executable(aeolion_geometry GeometryContractCLI.cpp)
target_link_libraries(aeolion_geometry PRIVATE aeolion aeolion_logger)
