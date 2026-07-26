# Aeolion BEMT: hover-safe blade-element momentum theory for a propeller,
# plus the slipstream field downstream control vanes read.
#
# STATIC, same rationale as panelbuilder/CMakeLists.txt: no ABI-stability
# or plugin boundary to justify a DLL's __declspec plumbing.
#
# Separated from the header-only `aeolion` target because it is a
# self-contained solver with real implementation weight (see BEMT.h),
# not a thin data-vocabulary header. PropGeometry.h stays in the
# header-only core precisely because it IS just data -- Geometry's
# HandoffContract builds a PropGeometry without ever calling into this
# solver.

add_library(aeolion_bemt STATIC
    src/BEMT.cpp
)

target_include_directories(aeolion_bemt PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# PUBLIC: anything using the solver also uses PropGeometry (core Geometry
# vocabulary), Math::Vec3, and the warnings/standard/LAPACK aeolion already
# carries -- inherited automatically since aeolion_enable_warnings() put
# those flags on aeolion's own INTERFACE.
target_link_libraries(aeolion_bemt PUBLIC aeolion)

if(AEOLION_BUILD_TESTS)
    add_subdirectory(tests)
endif()
