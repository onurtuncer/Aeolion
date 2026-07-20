// Geometry/MeshTopology.h
//
// Requested lattice discretization. This is the producer's recommendation for
// how finely to panel the surface, not a property of the geometry itself.
#pragma once

namespace Aeolion::Geometry {

enum class WakeModel { Frozen, Relaxed };

struct MeshTopology {
    int ChordwisePanels = 0;
    int SpanwisePanelsPerSection = 0; // panels between each adjacent station pair
    WakeModel Wake = WakeModel::Frozen;
};

} // namespace Aeolion::Geometry
