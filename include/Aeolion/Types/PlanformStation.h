// Types/PlanformStation.h
//
// One spanwise station of the handoff contract's planform. Stations are
// ordered by increasing Eta (0 at the root, 1 at the tip) and the quantities
// between them are interpolated linearly.
#pragma once

namespace Aeolion::Geometry {

struct PlanformStation {
    double Eta = 0.0;                  // normalized semi-span position, 0..1
    double Chord = 0.0;                // [m]
    double TwistDeg = 0.0;             // local geometric twist (JSON "twist")
    double SweepQuarterChordDeg = 0.0; // JSON "sweep_qc"
    double DihedralDeg = 0.0;          // JSON "dihedral"
};

} // namespace Aeolion::Geometry
