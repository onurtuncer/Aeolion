// Solver/ParticleTree.h
//
// Phase B of the particle-wake program: a Barnes-Hut octree over vector
// vortex particles, evaluating the regularized Biot-Savart velocity AND
// its gradient (the stretching term needs it) in O(N log N) instead of
// O(N^2). The Phase-A verdict was that configuration means fail for want
// of RESOLUTION, not physics -- this is the component that buys the
// resolution.
//
// Design, deliberately plain:
//
//   - The tree is rebuilt per time step over the particles' CURRENT
//     positions and serves every evaluation of that step (both RK2
//     stages sample the same source positions; only the targets move).
//   - Cells carry the MONOPOLE only: total vector strength at the
//     strength-weighted centroid. A far cell therefore evaluates exactly
//     like one big particle through the same kernels the direct path
//     uses -- no second code path for the physics, only for the
//     traversal. Accuracy comes from the opening criterion, and the
//     equivalence test pins it against direct summation.
//   - A cell is accepted when it is small seen from the target
//     (size/distance < Theta) AND the target is well clear of the
//     cell's largest smoothing core -- spreading cores must never be
//     approximated across, or diffusion would leak through the
//     multipole error.
//   - Leaves fall back to the exact pairwise kernel with symmetrized
//     cores, self-exclusion by particle index.
//
// The evaluator is read-only after Build, so callers may fan targets
// over threads freely.

#pragma once

#include "Aeolion/Math/Constants.h"
#include "Aeolion/Math/Vec3.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <numbers>
#include <vector>

namespace Aeolion::Solver {

inline constexpr double PtDefaultTheta = 0.5; ///< Opening criterion, size/distance.
inline constexpr int    PtLeafSize = 24;      ///< Direct summation below this population.
inline constexpr double PtCoreClearance = 9.0;///< Accept cells only beyond clearance^ (1/2) x core.

namespace Detail {

/** The particle fields the tree needs; mirrors WakeParticle's layout needs. */
struct TreeSource {
    Vec3 X{0, 0, 0};
    Vec3 Alpha{0, 0, 0};
    double Core2 = 0.0;
    int Original = -1; ///< Caller's index, for self-exclusion.
};

struct TreeCell {
    Vec3 Min{0, 0, 0}, Max{0, 0, 0};
    Vec3 Centroid{0, 0, 0};
    Vec3 AlphaSum{0, 0, 0};
    double MaxCore2 = 0.0;
    double Size2 = 0.0; ///< Squared diagonal, for the opening test.
    int Begin = 0, End = 0;
    int FirstChild = -1; ///< Eight consecutive children, or -1 for a leaf.
};

inline void TreeKernel(const Vec3& r, const Vec3& alpha, double pair2, Vec3& u, Vec3 grad[3],
                       bool wantGrad) {
    const double rho2 = Dot(r, r) + pair2;
    const double inv3 = 1.0 / (4.0 * std::numbers::pi * rho2 * std::sqrt(rho2));
    const Vec3 axr = Cross(alpha, r);
    u = u + axr * inv3;
    if (!wantGrad) return;
    const double inv5 = 3.0 * inv3 / rho2;
    const Vec3 unit[3] = {Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1)};
    for (int l = 0; l < 3; ++l) {
        const double rl = (l == 0) ? r.x : (l == 1) ? r.y : r.z;
        grad[l] = grad[l] + Cross(alpha, unit[l]) * inv3 - axr * (inv5 * rl);
    }
}

} // namespace Detail

class ParticleTree {
public:
    template <typename ParticleRange>
    void Build(const ParticleRange& particles) {
        const std::size_t np = particles.size();
        Sources.resize(np);
        for (std::size_t i = 0; i < np; ++i) {
            Sources[i].X = particles[i].X;
            Sources[i].Alpha = particles[i].Alpha;
            Sources[i].Core2 = particles[i].Core2;
            Sources[i].Original = static_cast<int>(i);
        }
        Cells.clear();
        if (np == 0) return;
        Cells.reserve(2 * np / PtLeafSize + 8);
        Detail::TreeCell root;
        root.Min = root.Max = Sources[0].X;
        for (const Detail::TreeSource& s : Sources) {
            root.Min.x = std::min(root.Min.x, s.X.x);
            root.Min.y = std::min(root.Min.y, s.X.y);
            root.Min.z = std::min(root.Min.z, s.X.z);
            root.Max.x = std::max(root.Max.x, s.X.x);
            root.Max.y = std::max(root.Max.y, s.X.y);
            root.Max.z = std::max(root.Max.z, s.X.z);
        }
        root.Begin = 0;
        root.End = static_cast<int>(np);
        Cells.push_back(root);
        Subdivide(0);
    }

    /**
     * Velocity (and, when grad is non-null, its gradient) induced at
     * `at` by every source except `excludeOriginal` (-1 excludes none).
     * `targetCore2` enters the symmetrized pair smoothing on direct pairs.
     */
    void Evaluate(const Vec3& at, double targetCore2, int excludeOriginal, Vec3& u,
                  Vec3* grad) const {
        if (Cells.empty()) return;
        EvaluateCell(0, at, targetCore2, excludeOriginal, u, grad);
    }

    double Theta = PtDefaultTheta;

private:
    std::vector<Detail::TreeSource> Sources;
    std::vector<Detail::TreeCell> Cells;

    void Subdivide(std::size_t cellIndex) {
        Detail::TreeCell& cell = Cells[cellIndex];
        // Moments first (every cell carries them, leaf or not).
        double weight = 0.0;
        Vec3 centroid(0, 0, 0);
        cell.AlphaSum = Vec3(0, 0, 0);
        cell.MaxCore2 = 0.0;
        for (int i = cell.Begin; i < cell.End; ++i) {
            const Detail::TreeSource& s = Sources[static_cast<std::size_t>(i)];
            const double w = std::max(s.Alpha.Norm(), Math::Tiny);
            weight += w;
            centroid = centroid + s.X * w;
            cell.AlphaSum = cell.AlphaSum + s.Alpha;
            cell.MaxCore2 = std::max(cell.MaxCore2, s.Core2);
        }
        cell.Centroid = centroid * (1.0 / weight);
        cell.Size2 = Dot(cell.Max - cell.Min, cell.Max - cell.Min);
        if (cell.End - cell.Begin <= PtLeafSize) return; // leaf

        const Vec3 mid = (cell.Min + cell.Max) * 0.5;
        // In-place octant partition: three successive binary partitions.
        std::array<int, 9> bounds{};
        bounds[0] = cell.Begin;
        bounds[8] = cell.End;
        const auto part = [&](int lo, int hi, auto pred) {
            int a = lo, b = hi;
            while (a < b) {
                if (pred(Sources[static_cast<std::size_t>(a)]))
                    ++a;
                else
                    std::swap(Sources[static_cast<std::size_t>(a)],
                              Sources[static_cast<std::size_t>(--b)]);
            }
            return a;
        };
        bounds[4] = part(bounds[0], bounds[8],
                         [&](const Detail::TreeSource& s) { return s.X.x < mid.x; });
        bounds[2] = part(bounds[0], bounds[4],
                         [&](const Detail::TreeSource& s) { return s.X.y < mid.y; });
        bounds[6] = part(bounds[4], bounds[8],
                         [&](const Detail::TreeSource& s) { return s.X.y < mid.y; });
        bounds[1] = part(bounds[0], bounds[2],
                         [&](const Detail::TreeSource& s) { return s.X.z < mid.z; });
        bounds[3] = part(bounds[2], bounds[4],
                         [&](const Detail::TreeSource& s) { return s.X.z < mid.z; });
        bounds[5] = part(bounds[4], bounds[6],
                         [&](const Detail::TreeSource& s) { return s.X.z < mid.z; });
        bounds[7] = part(bounds[6], bounds[8],
                         [&](const Detail::TreeSource& s) { return s.X.z < mid.z; });

        const int firstChild = static_cast<int>(Cells.size());
        Cells[cellIndex].FirstChild = firstChild;
        for (int o = 0; o < 8; ++o) {
            Detail::TreeCell child;
            child.Begin = bounds[static_cast<std::size_t>(o)];
            child.End = bounds[static_cast<std::size_t>(o) + 1];
            child.Min = Vec3((o & 4) ? mid.x : Cells[cellIndex].Min.x,
                            (o & 2) ? mid.y : Cells[cellIndex].Min.y,
                            (o & 1) ? mid.z : Cells[cellIndex].Min.z);
            child.Max = Vec3((o & 4) ? Cells[cellIndex].Max.x : mid.x,
                            (o & 2) ? Cells[cellIndex].Max.y : mid.y,
                            (o & 1) ? Cells[cellIndex].Max.z : mid.z);
            Cells.push_back(child);
        }
        for (int o = 0; o < 8; ++o) {
            const std::size_t ci = static_cast<std::size_t>(firstChild + o);
            if (Cells[ci].End > Cells[ci].Begin) Subdivide(ci);
        }
    }

    void EvaluateCell(std::size_t cellIndex, const Vec3& at, double targetCore2,
                      int excludeOriginal, Vec3& u, Vec3* grad) const {
        const Detail::TreeCell& cell = Cells[cellIndex];
        if (cell.End <= cell.Begin) return;
        const Vec3 r = at - cell.Centroid;
        const double d2 = Dot(r, r);
        const bool farEnough = cell.Size2 < Theta * Theta * d2 &&
                               d2 > PtCoreClearance * std::max(cell.MaxCore2, targetCore2);
        if (cell.FirstChild < 0 || farEnough) {
            if (farEnough) {
                Detail::TreeKernel(r, cell.AlphaSum,
                                   0.5 * (targetCore2 + cell.MaxCore2), u, grad,
                                   grad != nullptr);
                return;
            }
            for (int i = cell.Begin; i < cell.End; ++i) {
                const Detail::TreeSource& s = Sources[static_cast<std::size_t>(i)];
                if (s.Original == excludeOriginal) continue;
                Detail::TreeKernel(at - s.X, s.Alpha, 0.5 * (targetCore2 + s.Core2), u, grad,
                                   grad != nullptr);
            }
            return;
        }
        for (int o = 0; o < 8; ++o)
            EvaluateCell(static_cast<std::size_t>(cell.FirstChild + o), at, targetCore2,
                         excludeOriginal, u, grad);
    }
};

} // namespace Aeolion::Solver
