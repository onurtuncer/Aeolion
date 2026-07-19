// Math/Constants.h
//
// Shared numeric constants and helpers, so the rest of the toolkit never
// hard-codes a bare literal. Anything with a broader-than-one-module
// meaning lives here; module-specific tuning/empirical constants live in a
// constants block at the top of that module.

#pragma once
#include <numbers>

namespace Aeolion::Math {

// --- dimensionless factors that recur in formulae -------------------------
inline constexpr double Half           = 0.5;
inline constexpr double Two            = 2.0;
inline constexpr double QuarterChord   = 0.25; // x/c of the bound vortex line
inline constexpr double ThreeQuarterChord = 0.75; // x/c of the control point

// --- angle conversion -----------------------------------------------------
inline constexpr double DegreesPerHalfTurn = 180.0;
inline constexpr double SecondsPerMinute   = 60.0;

constexpr double DegToRad(double deg) { return deg * std::numbers::pi / DegreesPerHalfTurn; }
constexpr double RadToDeg(double rad) { return rad * DegreesPerHalfTurn / std::numbers::pi; }

// --- generic numerical guards ---------------------------------------------
inline constexpr double Tiny       = 1e-9;   // "practically zero" length/area guard
inline constexpr double UnitClampLo = -1.0;  // acos/asin argument clamp
inline constexpr double UnitClampHi = 1.0;

} // namespace Aeolion::Math
