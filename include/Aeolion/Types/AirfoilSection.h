// Types/AirfoilSection.h
//
// A section shape attached to a spanwise station, carried as Class-Shape
// Transformation (CST) coefficients rather than point coordinates.
//
// The CST surface is
//     y(psi) = psi^N1 * (1 - psi)^N2 * sum_i A_i * B_{i,n}(psi)
// with psi = x/c, B_{i,n} the Bernstein basis, and the class exponents fixed
// by convention at N1 = 0.5 / N2 = 1.0 -- the round-leading-edge,
// sharp-trailing-edge airfoil topology.
//
// The producer subtracts the linear psi * y_TE trailing-edge ramp before
// fitting, so these coefficients describe the SHARP-TE equivalent section.
// That is what a VLM camber surface consumes; there is deliberately no
// trailing-edge thickness field in the schema.
//
// The coefficients are the intended optimization design variables, so they
// are kept as fitted -- nothing here bakes them down into coordinates.
#pragma once

#include <vector>

namespace Aeolion::Geometry {

enum class SectionParameterization { CST };

// Class-function exponents implied by SectionParameterization::CST.
inline constexpr double CstN1 = 0.5;
inline constexpr double CstN2 = 1.0;

struct AirfoilSection {
    double Eta = 0.0; // normalized semi-span position, 0..1
    SectionParameterization Parameterization = SectionParameterization::CST;
    std::vector<double> CoefficientsUpper;
    std::vector<double> CoefficientsLower;
};

} // namespace Aeolion::Geometry
