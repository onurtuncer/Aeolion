// Geometry/MeshTopology.h
//
// Requested lattice discretization. This is the producer's recommendation for
// how finely to panel the surface, not a property of the geometry itself.
#pragma once

namespace Aeolion::Geometry {

// The schema's wake_model is validated at parse time but not stored: the
// solver implements exactly one wake (frozen, trailing legs along +x), so
// there is nothing for a consumer to branch on. A contract asking for a
// relaxed wake is REJECTED rather than quietly solved with a frozen one --
// see the parser. Store a wake model here only when more than one exists.
/** Requested lattice discretization. */
struct MeshTopology {
    int ChordwisePanels = 0;
    int SpanwisePanelsPerSection = 0; ///< Panels between each adjacent station pair.
};

} // namespace Aeolion::Geometry
