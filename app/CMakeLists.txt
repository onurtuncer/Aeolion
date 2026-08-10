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

# Runs the Journal of Aircraft paper's alpha/beta application sweep on the
# real geometry handoff (wing attachment line + fuselage skin-flow
# topology) and exports it as JSON for the figure renderer in
# papers/journal-of-aircraft/figures/.
add_executable(aeolion_attachment_sweep AttachmentSweepExport.cpp)
target_link_libraries(aeolion_attachment_sweep PRIVATE aeolion aeolion_panelbuilder)

# Exports the aft-fan induction paper's figure data (axis verification,
# chordwise gradient across the span, meridional field slice) from the
# vortex-cylinder model in Solver/DiskInduction.h.
add_executable(aeolion_fan_induction FanInductionExport.cpp)
target_link_libraries(aeolion_fan_induction PRIVATE aeolion aeolion_panelbuilder)

# Phase 0 of the post-separation study: the Level-2 coupled solve pushed
# over alpha in [-4, 90] x beta in [0, 30] with the current analytic
# deep-stall blend, exporting the flat-plate convergence metrics as JSON.
add_executable(aeolion_poststall_sweep PostStallSweepExport.cpp)
target_link_libraries(aeolion_poststall_sweep PRIVATE aeolion aeolion_panelbuilder)
