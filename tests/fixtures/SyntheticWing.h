// SyntheticWing.h -- test fixture: builds an in-memory thin-shell wing
// mesh (biconvex symmetric section, small nonzero thickness) for a given
// planform, so MeshSlice.h has a real upper/lower skin to slice. Shared
// by tests that need a mesh-derived wing to compare against the parametric
// VLM::BuildWing() path.
#pragma once
#include <numbers>
#include "Aeolion/Math/Vec3.h"
#include "Aeolion/Mesh.h"
#include <cmath>
#include <vector>

inline Aeolion::MeshIO::PartMesh BuildSyntheticWingMesh(double span, double rootChord, double tipChord,
                                                double sweepDeg, double dihedralDeg, double twistTipDeg,
                                                int nSpanSeg, int nChordSeg, double thicknessRatio,
                                                const std::string& name = "wing") {
    using Aeolion::VLM::Vec3; using Aeolion::VLM::Cross;
    Aeolion::MeshIO::PartMesh part;
    part.Name = name;
    double halfSpan = span / 2.0;
    double sweep = sweepDeg * std::numbers::pi / 180.0;
    double dihedral = dihedralDeg * std::numbers::pi / 180.0;

    auto chordAt = [&](double yAbs) { double eta = yAbs / halfSpan; return rootChord + (tipChord - rootChord) * eta; };
    auto qcX = [&](double y) { return std::fabs(y) * std::tan(sweep); };
    auto zOf = [&](double y) { return std::fabs(y) * std::tan(dihedral); };
    auto twistOf = [&](double y) { return (twistTipDeg * std::numbers::pi / 180.0) * (std::fabs(y) / halfSpan); };

    int NS = nSpanSeg, NC = nChordSeg;
    std::vector<std::vector<Vec3>> upper(NS + 1, std::vector<Vec3>(NC + 1));
    std::vector<std::vector<Vec3>> lower(NS + 1, std::vector<Vec3>(NC + 1));

    for (int k = 0; k <= NS; ++k) {
        double y = -halfSpan + span * k / NS;
        double c = chordAt(std::fabs(y));
        double xLE = qcX(y) - 0.25 * c;
        double eps = twistOf(y);
        double zc = zOf(y);
        for (int j = 0; j <= NC; ++j) {
            double u = (double)j / NC;
            double xLocal = u * c;
            double halfThick = 0.5 * thicknessRatio * c * 4.0 * u * (1.0 - u); // parabolic biconvex
            double ct = std::cos(eps), st = std::sin(eps);
            double xu = xLocal * ct + halfThick * st, zu = -xLocal * st + halfThick * ct;
            double xl = xLocal * ct - halfThick * st, zl = -xLocal * st - halfThick * ct;
            upper[k][j] = Vec3(xLE + xu, y, zc + zu);
            lower[k][j] = Vec3(xLE + xl, y, zc + zl);
        }
    }

    auto addTri = [&](const Vec3& a, const Vec3& b, const Vec3& c) {
        Vec3 n = Cross(b - a, c - a).Normalized();
        part.Tris.push_back({a, b, c, n});
    };
    for (int k = 0; k < NS; ++k) {
        for (int j = 0; j < NC; ++j) {
            addTri(upper[k][j], upper[k][j + 1], upper[k + 1][j]);
            addTri(upper[k][j + 1], upper[k + 1][j + 1], upper[k + 1][j]);
            addTri(lower[k][j], lower[k + 1][j], lower[k][j + 1]);
            addTri(lower[k][j + 1], lower[k + 1][j], lower[k + 1][j + 1]);
        }
    }
    return part;
}
