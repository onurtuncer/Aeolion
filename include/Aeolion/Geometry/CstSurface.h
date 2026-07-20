// Geometry/CstSurface.h
//
// Evaluation of the Class-Shape Transformation sections carried by the
// handoff (see AirfoilSection.h for the schema side and what the
// coefficients mean). AirfoilSection.h deliberately keeps the coefficients
// as fitted; this is where they get turned into geometry.
//
//     y(psi) = psi^N1 * (1 - psi)^N2 * sum_i A_i * B_{i,n}(psi)
//
// with psi = x/c measured aft from the leading edge, B_{i,n} the Bernstein
// basis of order n = (coefficient count - 1), and N1/N2 fixed by the
// schema at 0.5 / 1.0.
//
// A vortex lattice is a THIN lifting surface: it models the camber line,
// not thickness. So what the solver needs from a section is the mean line
//
//     camber(psi) = (y_upper(psi) + y_lower(psi)) / 2
//
// and its slope d(camber)/d(psi), which sets the surface normal a panel's
// flow-tangency condition is applied against. Both are normalized by chord
// -- multiply by the local chord to get metres.
//
// Sign convention: 'upper' is the +z (suction) side in the solver's axes
// (x aft, y right, z up), so positive camber bows upward and a cambered
// section lifts at zero angle of attack. The contract's own
// aetherion_body_frd frame applies to VECTOR fields like hinge axes (see
// ControlSurface.h), not to a 2D section's own chordwise ordinates.
//
// Because the producer subtracts the trailing-edge ramp before fitting
// (AirfoilSection.h), the class function's (1-psi)^N2 factor already
// drives both surfaces -- and therefore the camber line -- to exactly zero
// at psi=1. Nothing here re-adds a trailing-edge offset.

#pragma once

#include "Aeolion/Geometry/AirfoilSection.h"
#include "Aeolion/Math/Constants.h"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace Aeolion::Geometry {

// The class function's psi^(N1-1) derivative term is singular at the
// leading edge (N1 = 0.5), which is a real property of a round-nosed
// section, not a numerical defect. Lattice control points never sit at
// psi=0, but clamping keeps a caller that asks for the endpoint from
// getting an infinity back.
inline constexpr double CstEvaluationEpsilon = 1e-9;

namespace Detail {

[[nodiscard]] inline double BinomialCoefficient(std::size_t n, std::size_t i) {
    double result = 1.0;
    for (std::size_t k = 0; k < i; ++k)
        result = result * static_cast<double>(n - k) / static_cast<double>(k + 1);
    return result;
}

[[nodiscard]] inline double Bernstein(std::size_t i, std::size_t n, double psi) {
    return BinomialCoefficient(n, i) * std::pow(psi, static_cast<double>(i)) *
           std::pow(1.0 - psi, static_cast<double>(n - i));
}

[[nodiscard]] inline double BernsteinDerivative(std::size_t i, std::size_t n, double psi) {
    // d/dpsi [ C(n,i) psi^i (1-psi)^(n-i) ]; the i=0 and i=n edge terms
    // drop out rather than raising psi/(1-psi) to a negative power.
    double lower = 0.0;
    double upper = 0.0;
    if (i > 0)
        lower = static_cast<double>(i) * std::pow(psi, static_cast<double>(i) - 1.0) *
                std::pow(1.0 - psi, static_cast<double>(n - i));
    if (n > i)
        upper = static_cast<double>(n - i) * std::pow(psi, static_cast<double>(i)) *
                std::pow(1.0 - psi, static_cast<double>(n - i) - 1.0);
    return BinomialCoefficient(n, i) * (lower - upper);
}

[[nodiscard]] inline double ClassFunction(double psi) {
    return std::pow(psi, CstN1) * std::pow(1.0 - psi, CstN2);
}

[[nodiscard]] inline double ClassFunctionDerivative(double psi) {
    return CstN1 * std::pow(psi, CstN1 - 1.0) * std::pow(1.0 - psi, CstN2) -
           CstN2 * std::pow(psi, CstN1) * std::pow(1.0 - psi, CstN2 - 1.0);
}

[[nodiscard]] inline double ShapeFunction(const std::vector<double>& coefficients, double psi) {
    const std::size_t n = coefficients.size() - 1;
    double sum = 0.0;
    for (std::size_t i = 0; i < coefficients.size(); ++i) sum += coefficients[i] * Bernstein(i, n, psi);
    return sum;
}

[[nodiscard]] inline double ShapeFunctionDerivative(const std::vector<double>& coefficients, double psi) {
    const std::size_t n = coefficients.size() - 1;
    double sum = 0.0;
    for (std::size_t i = 0; i < coefficients.size(); ++i)
        sum += coefficients[i] * BernsteinDerivative(i, n, psi);
    return sum;
}

[[nodiscard]] inline double ClampPsi(double psi) {
    if (psi < CstEvaluationEpsilon) return CstEvaluationEpsilon;
    if (psi > 1.0 - CstEvaluationEpsilon) return 1.0 - CstEvaluationEpsilon;
    return psi;
}

} // namespace Detail

// --- one surface of one section -------------------------------------------
[[nodiscard]] inline double SurfaceOrdinate(const std::vector<double>& coefficients, double psi) {
    if (coefficients.empty()) return 0.0;
    const double p = Detail::ClampPsi(psi);
    return Detail::ClassFunction(p) * Detail::ShapeFunction(coefficients, p);
}

[[nodiscard]] inline double SurfaceSlope(const std::vector<double>& coefficients, double psi) {
    if (coefficients.empty()) return 0.0;
    const double p = Detail::ClampPsi(psi);
    return Detail::ClassFunctionDerivative(p) * Detail::ShapeFunction(coefficients, p) +
           Detail::ClassFunction(p) * Detail::ShapeFunctionDerivative(coefficients, p);
}

// --- camber line of one section -------------------------------------------
[[nodiscard]] inline double SectionCamber(const AirfoilSection& section, double psi) {
    return Math::Half * (SurfaceOrdinate(section.CoefficientsUpper, psi) +
                         SurfaceOrdinate(section.CoefficientsLower, psi));
}

[[nodiscard]] inline double SectionCamberSlope(const AirfoilSection& section, double psi) {
    return Math::Half * (SurfaceSlope(section.CoefficientsUpper, psi) +
                         SurfaceSlope(section.CoefficientsLower, psi));
}

// --- camber line at an arbitrary span station ------------------------------
// Sections are ordered by strictly increasing Eta and span [0,1] (the parser
// enforces both, see HandoffContract's RequireSpanningEtaSequence), so a
// linear search for the bracketing pair is well-defined.
//
// Adjacent sections may legitimately carry different CST orders, so this
// interpolates the EVALUATED ordinate/slope between the two bracketing
// sections rather than their coefficient vectors -- blending coefficients of
// unequal-order Bernstein bases would be meaningless.
namespace Detail {

struct SectionBracket {
    std::size_t Low = 0;
    std::size_t High = 0;
    double Weight = 0.0; // 0 at Low, 1 at High
};

[[nodiscard]] inline SectionBracket BracketByEta(const std::vector<AirfoilSection>& sections, double eta) {
    SectionBracket bracket;
    if (sections.size() == 1) return bracket; // Low == High == 0, weight 0
    for (std::size_t i = 0; i + 1 < sections.size(); ++i) {
        if (eta <= sections[i + 1].Eta || i + 2 == sections.size()) {
            bracket.Low = i;
            bracket.High = i + 1;
            const double span = sections[i + 1].Eta - sections[i].Eta;
            bracket.Weight = (span > 0.0) ? (eta - sections[i].Eta) / span : 0.0;
            if (bracket.Weight < 0.0) bracket.Weight = 0.0;
            if (bracket.Weight > 1.0) bracket.Weight = 1.0;
            return bracket;
        }
    }
    return bracket;
}

} // namespace Detail

[[nodiscard]] inline double CamberAt(const std::vector<AirfoilSection>& sections, double eta, double psi) {
    if (sections.empty()) return 0.0; // no section data -> flat plate
    const Detail::SectionBracket bracket = Detail::BracketByEta(sections, eta);
    const double low = SectionCamber(sections[bracket.Low], psi);
    const double high = SectionCamber(sections[bracket.High], psi);
    return low + (high - low) * bracket.Weight;
}

[[nodiscard]] inline double CamberSlopeAt(const std::vector<AirfoilSection>& sections, double eta, double psi) {
    if (sections.empty()) return 0.0;
    const Detail::SectionBracket bracket = Detail::BracketByEta(sections, eta);
    const double low = SectionCamberSlope(sections[bracket.Low], psi);
    const double high = SectionCamberSlope(sections[bracket.High], psi);
    return low + (high - low) * bracket.Weight;
}

// --- arc length along the camber line --------------------------------------
// Length of the mean line between two chordwise fractions, as a multiple of
// chord (multiply by the local chord for metres). This is what makes a
// panel's area the area of the actual curved surface rather than of its
// flat projection.
//
// Integrated in u = sqrt(psi) rather than directly in psi. A CAMBERED CST
// section's mean line has an inverse-square-root slope singularity at the
// round leading edge: the class function contributes psi^(N1-1) with
// N1 = 1/2, and while the upper and lower surfaces' leading-edge slopes
// cancel for a symmetric section, they do not once the section is cambered.
// The length integral stays finite (psi^-1/2 is integrable), but quadrature
// in psi would converge terribly against it. Under psi = u^2 the integrand
// sqrt(1 + s^2) dpsi becomes 2u*sqrt(1 + s^2) du, and since s ~ c/u near the
// nose that tends to the perfectly smooth 2*sqrt(u^2 + c^2).
//
// Composite Simpson in u is then exact for the flat case (integrand 2u) and
// converges fast for every real section.
inline constexpr int CamberArcQuadratureIntervals = 16; // even, for Simpson

[[nodiscard]] inline double CamberArcLengthFraction(const std::vector<AirfoilSection>& sections, double eta,
                                                    double psiStart, double psiEnd) {
    if (sections.empty()) return psiEnd - psiStart; // flat plate: arc length is the chord fraction
    const double uStart = std::sqrt(psiStart);
    const double uEnd = std::sqrt(psiEnd);
    if (!(uEnd > uStart)) return 0.0;

    const auto integrand = [&](double u) {
        const double slope = CamberSlopeAt(sections, eta, u * u);
        return 2.0 * u * std::sqrt(1.0 + slope * slope);
    };

    const double step = (uEnd - uStart) / CamberArcQuadratureIntervals;
    double sum = integrand(uStart) + integrand(uEnd);
    for (int i = 1; i < CamberArcQuadratureIntervals; ++i)
        sum += ((i % 2 == 1) ? 4.0 : 2.0) * integrand(uStart + i * step);
    return sum * step / 3.0;
}

} // namespace Aeolion::Geometry
