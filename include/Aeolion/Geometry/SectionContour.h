// Geometry/SectionContour.h
//
// The THICK closed contour of a CST section -- both surfaces joined into one
// loop -- as opposed to CstSurface.h's mean line.
//
// --- why the toolkit suddenly needs thickness --------------------------------
// Everything the vortex lattice does needs only the camber line, which is
// why CstSurface.h stops there: a lattice is a zero-thickness sheet, and the
// mean line is the whole of what it models. But a zero-thickness sheet has
// no stagnation point. The flow does not stop anywhere on it; its leading
// edge carries a square-root velocity singularity instead. So the moment
// the question becomes "where does the boundary layer start on this
// section", the camber line has nothing to say and the real contour does.
// The nose radius IS the answer's leading term -- a sharper nose puts the
// attachment point closer to the leading edge at the same incidence -- and
// a mean line has no nose radius at all.
//
// --- ordering, and why it is the one it is -----------------------------------
// Points run TRAILING EDGE -> LOWER surface -> LEADING EDGE -> UPPER surface
// -> trailing edge: one closed loop, traversed CLOCKWISE in the section's
// own (chordwise, upward) axes, so the interior lies to the right of the
// direction of travel. That is the orientation a Hess-Smith panel method
// wants, because it makes the outward normal of a panel whose tangent is at
// angle theta simply (-sin theta, cos theta) with no per-panel sign test.
// The first and last panels are then the two adjacent to the trailing edge,
// which is exactly the pair a Kutta condition is written on.
//
// --- spacing ------------------------------------------------------------------
// Cosine in psi, clustering at BOTH ends. The trailing-edge clustering is
// the usual accuracy argument; the leading-edge clustering is not optional
// here. The stagnation point sits within a nose radius or so of the leading
// edge -- a per-cent of chord on a typical section -- and uniform spacing
// would place the entire quantity of interest inside the first panel.
//
// Everything is normalized by chord: multiply by the local chord for metres.

#pragma once

#include "Aeolion/Geometry/AirfoilSection.h"
#include "Aeolion/Geometry/CstSurface.h"
#include "Aeolion/Math/Constants.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <vector>

namespace Aeolion::Geometry {

inline constexpr int DefaultContourPointsPerSurface = 80; ///< Per surface, so ~160 panels a section.
inline constexpr int MinContourPointsPerSurface = 8;

/**
 * A closed section contour in the section's own 2-D axes: Psi aft from the
 * leading edge, Zeta upward, both as fractions of chord.
 *
 * Ordered trailing edge -> lower -> leading edge -> upper -> trailing edge,
 * clockwise, with the leading-edge point appearing once. Point 0 and the
 * last point are the same trailing-edge location; the loop is closed
 * explicitly rather than implied, so a panel method can walk it without a
 * wrap-around special case.
 */
struct SectionContour {
    std::vector<double> Psi;  ///< x/c, aft from the leading edge.
    std::vector<double> Zeta; ///< z/c, positive toward the suction side.

    /**
     * Leading-edge radius as a fraction of chord.
     *
     * Free from the CST parameterization rather than fitted from the
     * points: with the class exponent N1 = 1/2 the surface behaves like
     * sqrt(psi) * A_0 near the nose, and a curve z = A_0 sqrt(psi) is a
     * parabola of radius A_0^2 / 2 at its vertex. Upper and lower give the
     * same answer for any sanely fitted section; the mean of the two is
     * taken so a slightly mismatched fit degrades rather than picks a side.
     */
    double LeadingEdgeRadius = 0.0;

    [[nodiscard]] std::size_t Count() const { return Psi.size(); }
    [[nodiscard]] bool Valid() const { return Psi.size() >= 4 && Psi.size() == Zeta.size(); }
};

/**
 * Leading-edge radius of one CST surface, r_LE/c = A_0^2 / 2 (see
 * SectionContour::LeadingEdgeRadius).
 */
[[nodiscard]] inline double SurfaceLeadingEdgeRadius(const std::vector<double>& coefficients) {
    if (coefficients.empty()) return 0.0;
    return Math::Half * coefficients.front() * coefficients.front();
}

/**
 * Build the closed contour of one CST section. `pointsPerSurface` counts the
 * points on each surface between the leading and trailing edges.
 */
[[nodiscard]] inline SectionContour BuildSectionContour(
    const AirfoilSection& section, int pointsPerSurface = DefaultContourPointsPerSurface) {
    SectionContour contour;
    const int n = std::max(pointsPerSurface, MinContourPointsPerSurface);

    // Cosine spacing, leading edge at psi = 0.
    const auto psiAt = [n](int k) {
        return Math::Half * (1.0 - std::cos(std::numbers::pi * static_cast<double>(k) / n));
    };

    contour.Psi.reserve(static_cast<std::size_t>(2 * n + 1));
    contour.Zeta.reserve(static_cast<std::size_t>(2 * n + 1));

    // Trailing edge -> leading edge along the LOWER surface.
    for (int k = n; k >= 0; --k) {
        const double psi = psiAt(k);
        contour.Psi.push_back(psi);
        contour.Zeta.push_back(SurfaceOrdinate(section.CoefficientsLower, psi));
    }
    // Leading edge -> trailing edge along the UPPER surface, skipping the
    // leading-edge point already placed.
    for (int k = 1; k <= n; ++k) {
        const double psi = psiAt(k);
        contour.Psi.push_back(psi);
        contour.Zeta.push_back(SurfaceOrdinate(section.CoefficientsUpper, psi));
    }

    // The class function drives both surfaces to zero at psi = 1 (the
    // producer removed the trailing-edge ramp before fitting, see
    // AirfoilSection.h), so the loop closes on itself and the sharp
    // trailing edge a Kutta condition assumes is real, not imposed here.
    const double upper = SurfaceLeadingEdgeRadius(section.CoefficientsUpper);
    const double lower = SurfaceLeadingEdgeRadius(section.CoefficientsLower);
    const int surfaces = (section.CoefficientsUpper.empty() ? 0 : 1) +
                         (section.CoefficientsLower.empty() ? 0 : 1);
    contour.LeadingEdgeRadius = (surfaces > 0) ? (upper + lower) / surfaces : 0.0;
    return contour;
}

/**
 * The contour at an arbitrary span station, interpolated between bracketing
 * sections the same way CamberAt does -- on the EVALUATED ordinates, because
 * adjacent sections may carry different CST orders and blending coefficients
 * of unequal-order Bernstein bases would be meaningless.
 */
[[nodiscard]] inline SectionContour ContourAt(const std::vector<AirfoilSection>& sections, double eta,
                                              int pointsPerSurface = DefaultContourPointsPerSurface) {
    SectionContour contour;
    if (sections.empty()) return contour; // no section data: no thickness to speak of

    const Detail::SectionBracket bracket = Detail::BracketByEta(sections, eta);
    const SectionContour low = BuildSectionContour(sections[bracket.Low], pointsPerSurface);
    const SectionContour high = BuildSectionContour(sections[bracket.High], pointsPerSurface);
    if (!low.Valid() || low.Count() != high.Count()) return low;

    contour.Psi = low.Psi; // identical spacing by construction
    contour.Zeta.resize(low.Count());
    for (std::size_t k = 0; k < low.Count(); ++k)
        contour.Zeta[k] = low.Zeta[k] + (high.Zeta[k] - low.Zeta[k]) * bracket.Weight;
    contour.LeadingEdgeRadius =
        low.LeadingEdgeRadius + (high.LeadingEdgeRadius - low.LeadingEdgeRadius) * bracket.Weight;
    return contour;
}

/**
 * Transform a contour into the plane NORMAL to a swept leading edge.
 *
 * Sections in the handoff are streamwise: they are cut at a span fraction,
 * along the flight direction. The attachment line of a swept wing does not
 * live in that plane. Infinite-swept-wing theory poses the problem in the
 * plane perpendicular to the leading edge, where the chord shortens by
 * cos(sweep) while the ordinates do not, so the section gets THICKER and its
 * nose blunter in proportion:
 *
 *     psi_n = psi,   zeta_n = zeta / cos(Lambda),   r_LE,n = r_LE / cos^2(Lambda)
 *
 * once both are renormalized by the shortened normal chord c_n = c cos(Lambda).
 * Taking the streamwise section instead would put the attachment point too
 * close to the leading edge on a swept wing, and would get the sign of its
 * movement with sideslip wrong on one side of the aircraft -- which is the
 * whole point of doing this in the normal plane.
 */
[[nodiscard]] inline SectionContour ToLeadingEdgeNormal(const SectionContour& streamwise,
                                                        double sweepRad) {
    SectionContour normal = streamwise;
    const double cosSweep = std::cos(sweepRad);
    if (!(std::fabs(cosSweep) > Math::Tiny)) return normal; // fully swept: no normal plane left
    const double stretch = 1.0 / cosSweep;
    for (double& zeta : normal.Zeta) zeta *= stretch;
    normal.LeadingEdgeRadius = streamwise.LeadingEdgeRadius * stretch * stretch;
    return normal;
}

/** Maximum thickness of a contour as a fraction of chord, upper minus lower. */
[[nodiscard]] inline double MaxThickness(const SectionContour& contour) {
    if (!contour.Valid()) return 0.0;
    // The loop is symmetric about its midpoint index: point k from the start
    // (lower) pairs with point Count()-1-k from the end (upper) at the same psi.
    double thickest = 0.0;
    const std::size_t count = contour.Count();
    for (std::size_t k = 0; k < count / 2; ++k) {
        const std::size_t mirror = count - 1 - k;
        if (std::fabs(contour.Psi[k] - contour.Psi[mirror]) > Math::Tiny) continue;
        thickest = std::max(thickest, contour.Zeta[mirror] - contour.Zeta[k]);
    }
    return thickest;
}

} // namespace Aeolion::Geometry
