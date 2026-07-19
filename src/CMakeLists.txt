# Driver executables. Each links the header-only `aeolion` interface library,
# which carries the include path, C++20 requirement, warnings, and LAPACK.

add_executable(vlm_demo main.cpp)
target_link_libraries(vlm_demo PRIVATE aeolion)

add_executable(aeolion_geometry GeometryContractCLI.cpp)
target_link_libraries(aeolion_geometry PRIVATE aeolion)
