# Driver executables. Each links the header-only `aeolion` interface library,
# which carries the include path, C++20 requirement, warnings, and LAPACK.

add_executable(vlm_demo main.cpp)
target_link_libraries(vlm_demo PRIVATE aeolion)

add_executable(geometry_vlm GeometryVLM.cpp)
target_link_libraries(geometry_vlm PRIVATE aeolion)
