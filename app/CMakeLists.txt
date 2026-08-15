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

# Section-level validation of the anchored post-stall model against the
# Sheldahl-Klimas NACA 0015 table (Part II's declared acceptance test):
# CST-fits the section, computes f(alpha) through the same Hess-Smith +
# attachment-march pipeline as the airframe, exports the model polar.
add_executable(aeolion_section_validation SectionValidationExport.cpp)
target_link_libraries(aeolion_section_validation PRIVATE aeolion)

# The aeroCD0 table of the DAVE-ML flight model: the parasite drag of the
# body and duct (never the wing -- its profile drag is already inside the
# coupled solve), component buildup plus a slender-body crossflow branch
# so the estimate stays meaningful to 90 degrees incidence.
add_executable(aeolion_parasite_drag ParasiteDragExport.cpp)
target_link_libraries(aeolion_parasite_drag PRIVATE aeolion aeolion_panelbuilder)

# The propulsion side of the DAVE-ML flight model (models/README.md): the
# contract's ducted propulsor swept over advance ratio and the three vane
# command modes, through the two-way rotor-vane coupling, exported as JSON
# for the DAVE-ML assembler.
add_executable(aeolion_propulsion_map PropulsionMapExport.cpp)
target_link_libraries(aeolion_propulsion_map PRIVATE aeolion aeolion_panelbuilder)

# The 3-D unsteady cross-check: the particle-wake solve at the declared
# attitudes (alpha 30/60/90 x beta 0/15) on the handoff wing -- full-span
# and coarse, built from the contract's planform law (see the driver's
# header for the measured reasons) -- mean and RMS loads against the
# quasi-steady map's cycle means.
add_executable(aeolion_particle_crosscheck ParticleCrosscheckExport.cpp)
target_link_libraries(aeolion_particle_crosscheck PRIVATE aeolion)
find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
  target_link_libraries(aeolion_particle_crosscheck PRIVATE OpenMP::OpenMP_CXX)
endif()
