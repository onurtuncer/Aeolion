// Math/Vec3.h
//
// The one math primitive the toolkit is built on: a 3D double vector plus
// the handful of free functions (dot, cross, axis rotations) that operate
// on it. Lives in namespace Aeolion::Math.
//
// For source-compatibility the names are also re-exported into
// Aeolion::Solver, so existing Solver::Vec3 / Solver::Dot / Solver::Cross
// references keep resolving -- Solver is the umbrella that these math types
// feed into.

#pragma once
#include <cmath>

namespace Aeolion::Math {

struct Vec3 {
    double x = 0.0, y = 0.0, z = 0.0;

    Vec3() = default;
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    [[nodiscard]] Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    [[nodiscard]] Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    [[nodiscard]] Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
    [[nodiscard]] Vec3 operator-() const { return {-x, -y, -z}; }

    [[nodiscard]] double Norm() const { return std::sqrt(x * x + y * y + z * z); }

    [[nodiscard]] Vec3 Normalized() const {
        double n = Norm();
        if (n < 1e-14) return {0, 0, 0};
        return {x / n, y / n, z / n};
    }
};

[[nodiscard]] inline double Dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

[[nodiscard]] inline Vec3 Cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

[[nodiscard]] inline Vec3 RotateAboutX(const Vec3& v, double angleRad) {
    double c = std::cos(angleRad), s = std::sin(angleRad);
    return {v.x, c * v.y - s * v.z, s * v.y + c * v.z};
}

// Rotate about a local spanwise axis (approximated as global y) — used to
// apply geometric twist to the panel normal vector.
[[nodiscard]] inline Vec3 RotateAboutY(const Vec3& v, double angleRad) {
    double c = std::cos(angleRad), s = std::sin(angleRad);
    return {c * v.x + s * v.z, v.y, -s * v.x + c * v.z};
}

} // namespace Aeolion::Math

// Re-export the math vocabulary into Aeolion::Solver (see file header).
namespace Aeolion::Solver {
    using Math::Vec3;
    using Math::Dot;
    using Math::Cross;
    using Math::RotateAboutX;
    using Math::RotateAboutY;
} // namespace Aeolion::Solver
