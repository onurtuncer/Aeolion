# Aeolion PanelBuilder: builds the solver's lattice from a parsed geometry
# handoff -- cambered CST surface, spanwise and chordwise discretization,
# deflected control surfaces.
#
# STATIC rather than SHARED on purpose. A DLL here would need
# __declspec(dllexport/dllimport) plumbing on every exported symbol and a
# runtime deployment story, and would buy nothing: there is no ABI-stability
# requirement, no plugin boundary, and no separate release cadence. The one
# real benefit of a shared build -- not recompiling consumers -- does not
# apply when the type vocabulary (Lattice::Panel, Geometry::HandoffContract)
# is header-only anyway.
#
# Separated from the header-only `aeolion` target because it is the one
# component with a compiled part; consumers link it in addition to aeolion.

add_library(aeolion_panelbuilder STATIC
    src/PanelBuilder.cpp
)

target_include_directories(aeolion_panelbuilder PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# PUBLIC: anything using the builder also uses Lattice::Panel and the
# contract types it takes and returns.
target_link_libraries(aeolion_panelbuilder PUBLIC aeolion)
