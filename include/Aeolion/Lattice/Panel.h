// Lattice/Panel.h
//
// One horseshoe-vortex panel: the atomic geometric unit of the lattice.
//
// This is shared vocabulary, not solver-owned: PanelBuilder produces
// panels, Solver consumes them, and the viewer draws them. It therefore
// lives in its own small module rather than inside any one of them, the
// same way Math::Vec3 does -- and, like Vec3, it is re-exported into
// Aeolion::Solver at the bottom of this file so Solver::Panel keeps
// resolving for code that thinks of it as part of the solver's API.
//
// A lattice may be a single spanwise row (one panel per span station, the
// flat wing BuildWing() produces) or several chordwise rows stacked on a
// cambered surface (what PanelBuilder builds from a geometry handoff). In
// the second case the panel's quarter-chord/three-quarter-chord points are
// those of ITS OWN chordwise strip of the section, not of the whole chord,
// and StripIndex ties the stack back together for per-span reporting.
#pragma once
#include <string>
#include "Aeolion/Math/Vec3.h"

namespace Aeolion::Lattice {

using Math::Vec3;

struct Panel {
    Vec3 A, B;              // bound vortex endpoints (panel quarter-chord), A.y < B.y
    Vec3 ControlPoint;      // panel three-quarter-chord, mid-span
    Vec3 Normal;            // unit outward normal at control point (twist + dihedral + camber slope)
    Vec3 TrailDirA, TrailDirB; // unit direction of trailing legs from A and from B (downstream)
    // Area of the actual (cambered, dihedral-tilted) panel surface, and the
    // area of its projection onto the reference plane. They are equal only
    // for a flat, un-dihedralled panel.
    //
    // These are NOT interchangeable. Force and moment coefficients are
    // normalized by PLANFORM area by universal convention -- CL = L/(q*S)
    // with S projected -- so every coefficient in the solver divides by
    // PlanformArea. Area is the true wetted-side geometry, which is what a
    // skin-friction / viscous buildup wants (see DragEstimate). Normalizing
    // a coefficient by Area instead would silently redefine CL and break
    // comparison against every published lift curve.
    double Area = 0.0;
    double PlanformArea = 0.0;
    double SpanwiseWidth = 0.0;
    std::string Surface;    // which lifting surface this panel belongs to (e.g. "wing", "htail")

    // Panels sharing a StripIndex are the chordwise stack covering ONE
    // spanwise station, and the solver reports them as a single
    // StationResult (summed circulation and lift over the section chord).
    // Negative means "this panel is its own station" -- the single-row
    // case, where per-panel and per-station are the same thing.
    int StripIndex = -1;
};

} // namespace Aeolion::Lattice

// Re-export into the solver's namespace (see file header).
namespace Aeolion::Solver {
    using Lattice::Panel;
} // namespace Aeolion::Solver
