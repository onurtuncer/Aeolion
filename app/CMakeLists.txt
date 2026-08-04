# Driver executables. Each links the header-only `aeolion` interface library,
# which carries the include path, C++23 requirement, warnings, and LAPACK,
# plus the logger every application reports through.

add_executable(solver_demo main.cpp)
target_link_libraries(solver_demo PRIVATE aeolion aeolion_logger)

add_executable(aeolion_geometry GeometryContractCLI.cpp)
target_link_libraries(aeolion_geometry PRIVATE aeolion aeolion_logger)

# Exports the paper's reference ducted-fan case (lattice + solved
# circulation) as JSON for the publication figure renderer in
# papers/aiaa-scitech/figures/.
add_executable(aeolion_lattice_export LatticeFigureExport.cpp)
target_link_libraries(aeolion_lattice_export PRIVATE aeolion aeolion_panelbuilder)
