// Solver/SourceInfluence.h
//
// Velocity induced by a constant-strength quadrilateral SOURCE panel --
// the body-surface counterpart to Solver.h's horseshoe Biot-Savart kernel.
//
// The in-plane components are the classical Hess & Smith flat-panel result
// (Katz & Plotkin, "Low-Speed Aerodynamics" 10.4.1): transform into a local
// frame where the panel lies in z = 0 and sum a logarithm per edge.
//
// The NORMAL component is deliberately not computed by the textbook's
// per-edge arctangent form. That form carries a slope (y2-y1)/(x2-x1) which
// blows up for any edge parallel to the local y axis -- a case that arises
// constantly on a body of revolution, where panel edges align with the
// axis. Instead the normal component is the signed SOLID ANGLE the panel
// subtends at the field point, which is the same quantity and has a
// singularity-free closed form:
//
//     w = Omega / (4 pi),   Omega = solid angle of the panel
//
// computed per triangle by Van Oosterom & Strackee's formula
//
//     tan(Omega/2) = a . (b x c) / (|a||b||c| + (a.b)|c| + (a.c)|b| + (b.c)|a|)
//
// evaluated with a two-argument arctangent so the quadrant is never lost.
// On the sheet Omega -> +/- 2 pi and w -> +/- 1/2, which is exactly the
// analytic source-sheet jump, recovered rather than special-cased.
//
// The remaining care is the logarithm: r1 + r2 - d vanishes for a field
// point lying on an edge's own line. That divergence is integrable and
// physically harmless -- it is the rim of a sheet, not a real singularity --
// so it is floored rather than allowed to put an infinity into the
// influence matrix.
//
// Panel corners are assumed near-planar and are projected into the local
// plane: exact for a flat quad, second-order for the slightly warped quads
// a body of revolution produces.

#pragma once

#include "Aeolion/Lattice/SourcePanel.h"
#include "Aeolion/Math/Constants.h"
#include "Aeolion/Math/Vec3.h"

#include <array>
#include <cmath>
#include <numbers>

namespace Aeolion::Solver {

using Lattice::SourcePanel;

// --- tolerances -------------------------------------------------------------
inline constexpr double SourceEdgeLengthEps = 1e-12;    // a zero-length edge contributes nothing
inline constexpr double SourceLogArgumentFloor = 1e-12; // keeps the on-edge logarithm finite
inline constexpr double SourceSolidAngleEps = 1e-30;    // degenerate triangle guard
inline constexpr double InvFourPiSource = 1.0 / (4.0 * std::numbers::pi);

// Self-induced normal velocity of a constant-strength source sheet on the
// sheet itself. Outward, and exactly half the strength.
inline constexpr double SourceSelfInfluence = 0.5;

/** The panel's local orthonormal frame: E1/E2 span the panel plane, E3 is the outward normal. */
struct SourcePanelFrame {
    Vec3 Origin; ///< Panel centroid.
    Vec3 E1, E2, E3;
};

/** Build a SourcePanelFrame for panel. */
[[nodiscard]] inline SourcePanelFrame MakeSourcePanelFrame(const SourcePanel& panel) {
    SourcePanelFrame frame;
    frame.Origin = panel.ControlPoint;
    frame.E3 = panel.Normal.Normalized();

    // Any in-plane direction will do; take the first edge and orthogonalize
    // it against the normal so the frame stays exactly orthonormal even when
    // the quad is slightly warped.
    Vec3 along = panel.Corners[1] - panel.Corners[0];
    along = along - frame.E3 * Dot(along, frame.E3);
    if (along.Norm() < SourceEdgeLengthEps) {
        along = panel.Corners[2] - panel.Corners[0];
        along = along - frame.E3 * Dot(along, frame.E3);
    }
    frame.E1 = along.Normalized();
    frame.E2 = Cross(frame.E3, frame.E1);
    return frame;
}

namespace Detail {

// Signed solid angle subtended at the origin by the triangle (a, b, c),
// Van Oosterom & Strackee. The sign follows the winding of a->b->c.
[[nodiscard]] inline double TriangleSolidAngle(const Vec3& a, const Vec3& b, const Vec3& c) {
    const double la = a.Norm();
    const double lb = b.Norm();
    const double lc = c.Norm();
    const double numerator = Dot(a, Cross(b, c));
    const double denominator =
        la * lb * lc + Dot(a, b) * lc + Dot(a, c) * lb + Dot(b, c) * la;
    if (std::fabs(numerator) < SourceSolidAngleEps && std::fabs(denominator) < SourceSolidAngleEps)
        return 0.0;
    return 2.0 * std::atan2(numerator, denominator);
}

} // namespace Detail

/**
 * Velocity at `point` induced by `panel` carrying UNIT source strength.
 * Multiply by the panel's sigma for the physical velocity.
 */
[[nodiscard]] inline Vec3 SourcePanelVelocity(const Vec3& point, const SourcePanel& panel) {
    const SourcePanelFrame frame = MakeSourcePanelFrame(panel);

    const Vec3 offset = point - frame.Origin;
    const double x = Dot(offset, frame.E1);
    const double y = Dot(offset, frame.E2);
    const double z = Dot(offset, frame.E3);

    std::array<double, 4> cornerX{};
    std::array<double, 4> cornerY{};
    for (std::size_t k = 0; k < 4; ++k) {
        const Vec3 local = panel.Corners[k] - frame.Origin;
        cornerX[k] = Dot(local, frame.E1);
        cornerY[k] = Dot(local, frame.E2);
    }

    // --- in-plane components: one logarithm per edge ------------------------
    double u = 0.0;
    double v = 0.0;
    for (std::size_t k = 0; k < 4; ++k) {
        const std::size_t next = (k + 1) % 4;
        const double x1 = cornerX[k], y1 = cornerY[k];
        const double x2 = cornerX[next], y2 = cornerY[next];

        const double dx = x2 - x1;
        const double dy = y2 - y1;
        const double d = std::sqrt(dx * dx + dy * dy);
        if (d < SourceEdgeLengthEps) continue;

        const double r1 = std::sqrt((x - x1) * (x - x1) + (y - y1) * (y - y1) + z * z);
        const double r2 = std::sqrt((x - x2) * (x - x2) + (y - y2) * (y - y2) + z * z);

        // ln((r1+r2+d)/(r1+r2-d)), i.e. the log of a quantity >= 1, so the
        // term is positive and a source pushes fluid outward. The small
        // factor is the one that vanishes for a field point on the edge's
        // own line, so it is the one that gets floored.
        double small = r1 + r2 - d;
        if (small < SourceLogArgumentFloor) small = SourceLogArgumentFloor;
        const double logTerm = std::log((r1 + r2 + d) / small);

        // Edge length in the denominator, not the edge slope, so an edge
        // parallel to either local axis needs no special case.
        u += (dy / d) * logTerm;
        v -= (dx / d) * logTerm;
    }

    // --- normal component: the solid angle the quad subtends ----------------
    // Split into two triangles about corner 0. Vectors run from the FIELD
    // POINT to the corners, so the solid angle is the one seen from there.
    const Vec3 r0 = panel.Corners[0] - point;
    const Vec3 r1 = panel.Corners[1] - point;
    const Vec3 r2 = panel.Corners[2] - point;
    const Vec3 r3 = panel.Corners[3] - point;
    const double solidAngle =
        Detail::TriangleSolidAngle(r0, r1, r2) + Detail::TriangleSolidAngle(r0, r2, r3);

    // Winding is right-handed about Normal, so a point on the +Normal side
    // sees a negative solid angle by the convention above; flip it so the
    // induced normal velocity is outward, as a source must be.
    const double w = -solidAngle;

    return frame.E1 * (u * InvFourPiSource) + frame.E2 * (v * InvFourPiSource) +
           frame.E3 * (w * InvFourPiSource);
}

/**
 * Normal-velocity influence coefficient: the normal component, at
 * `target`'s control point, of the velocity induced by unit source strength
 * on `source`. One entry of the body block of the influence matrix.
 */
[[nodiscard]] inline double SourceNormalInfluence(const SourcePanel& target, const SourcePanel& source) {
    // A panel's own control point lies exactly on its sheet, where the solid
    // angle is indeterminate. The answer there is analytic.
    if (&target == &source) return SourceSelfInfluence;
    return Dot(SourcePanelVelocity(target.ControlPoint, source), target.Normal);
}

} // namespace Aeolion::Solver
