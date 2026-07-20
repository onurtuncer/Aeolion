// TestSourcePanel.cpp -- validates the constant-strength source panel kernel
// against results that are known in closed form, not against itself:
//
//   - Far field: at range, a panel of area A and unit strength must look
//     like a point source of strength A, i.e. radial flow of magnitude
//     A/(4*pi*r^2). This pins the overall scale factor AND the sign.
//   - On the sheet: the normal velocity just off a panel's own face must be
//     +/- 1/2, the analytic source-sheet jump, approached from both sides.
//   - Sphere: a source-panelled sphere in uniform flow must reproduce
//     potential flow's exact answer, surface speed = (3/2) U sin(theta),
//     equivalently Cp = 1 - (9/4) sin^2(theta). This is the real test: it
//     exercises the kernel, the influence matrix and the solve together,
//     and it is sensitive to any error in the solid-angle sign convention.
#include "Aeolion/Solver/SourceInfluence.h"
#include "Aeolion/Solver/Solver.h"

#include <cmath>
#include <iostream>
#include <numbers>
#include <string>
#include <vector>

using Aeolion::Lattice::SourcePanel;
using Aeolion::Math::Vec3;

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } } while (0)

namespace {

// A flat unit square in the z = 0 plane, normal +z, wound right-handed.
SourcePanel UnitSquare() {
    SourcePanel panel;
    panel.Corners = {Vec3(-0.5, -0.5, 0.0), Vec3(0.5, -0.5, 0.0), Vec3(0.5, 0.5, 0.0), Vec3(-0.5, 0.5, 0.0)};
    panel.ControlPoint = Vec3(0, 0, 0);
    panel.Normal = Vec3(0, 0, 1);
    panel.Area = 1.0;
    return panel;
}

void TestFarFieldLooksLikeAPointSource() {
    const SourcePanel panel = UnitSquare();

    // Far enough that the panel's shape is irrelevant.
    for (double r : {20.0, 50.0, 200.0}) {
        for (const Vec3& direction : {Vec3(0, 0, 1), Vec3(1, 0, 0), Vec3(0.577, 0.577, 0.577)}) {
            const Vec3 unit = direction.Normalized();
            const Vec3 point = unit * r;
            const Vec3 velocity = Aeolion::Solver::SourcePanelVelocity(point, panel);

            const double expected = panel.Area / (4.0 * std::numbers::pi * r * r);
            const double radial = Aeolion::Math::Dot(velocity, unit);

            // A source pushes fluid OUTWARD: radial component positive.
            CHECK(radial > 0.0, "far-field velocity must point away from the source panel, r=" + std::to_string(r));
            CHECK(std::fabs(radial - expected) < 0.02 * expected,
                  "far field should match a point source of strength A at r=" + std::to_string(r) +
                      ": got " + std::to_string(radial) + ", expected " + std::to_string(expected));
            // And it must be essentially purely radial out there.
            const double tangential = (velocity - unit * radial).Norm();
            CHECK(tangential < 0.02 * expected, "far field should be radial at r=" + std::to_string(r));
        }
    }
}

void TestSheetJump() {
    const SourcePanel panel = UnitSquare();

    // Approach the face from both sides. The normal velocity must tend to
    // +/- 1/2 -- the analytic jump across a constant-strength source sheet.
    for (double offset : {1e-3, 1e-4, 1e-5}) {
        const double above =
            Aeolion::Math::Dot(Aeolion::Solver::SourcePanelVelocity(Vec3(0, 0, offset), panel), panel.Normal);
        const double below =
            Aeolion::Math::Dot(Aeolion::Solver::SourcePanelVelocity(Vec3(0, 0, -offset), panel), panel.Normal);

        CHECK(std::fabs(above - 0.5) < 0.01,
              "normal velocity just above the sheet should be +1/2, got " + std::to_string(above));
        CHECK(std::fabs(below + 0.5) < 0.01,
              "normal velocity just below the sheet should be -1/2, got " + std::to_string(below));
    }

    // The coefficient used on the matrix diagonal must agree with that limit.
    CHECK(std::fabs(Aeolion::Solver::SourceNormalInfluence(panel, panel) - 0.5) < 1e-15,
          "the self-influence coefficient must be exactly 1/2");
}

// --- sphere in uniform flow -------------------------------------------------
std::vector<SourcePanel> BuildSphere(double radius, int polarPanels, int azimuthPanels) {
    const auto position = [radius](double theta, double phi) {
        // theta measured from +x, which is the freestream direction below.
        return Vec3(radius * std::cos(theta), radius * std::sin(theta) * std::cos(phi),
                    radius * std::sin(theta) * std::sin(phi));
    };

    std::vector<SourcePanel> panels;
    panels.reserve(static_cast<std::size_t>(polarPanels * azimuthPanels));

    for (int i = 0; i < polarPanels; ++i) {
        const double theta0 = std::numbers::pi * i / polarPanels;
        const double theta1 = std::numbers::pi * (i + 1) / polarPanels;
        for (int j = 0; j < azimuthPanels; ++j) {
            const double phi0 = 2.0 * std::numbers::pi * j / azimuthPanels;
            const double phi1 = 2.0 * std::numbers::pi * (j + 1) / azimuthPanels;

            SourcePanel panel;
            panel.Corners = {position(theta0, phi0), position(theta1, phi0), position(theta1, phi1),
                             position(theta0, phi1)};

            Vec3 centroid(0, 0, 0);
            for (const Vec3& corner : panel.Corners) centroid = centroid + corner * 0.25;
            panel.ControlPoint = centroid;

            // Outward normal from the quad itself, with the winding fixed to
            // match it -- the kernel's solid angle is winding-sensitive.
            Vec3 normal = Aeolion::Math::Cross(panel.Corners[2] - panel.Corners[0],
                                               panel.Corners[3] - panel.Corners[1]);
            if (Aeolion::Math::Dot(normal, centroid) < 0.0) {
                std::swap(panel.Corners[1], panel.Corners[3]);
                normal = normal * -1.0;
            }
            panel.Normal = normal.Normalized();
            panel.Area = 0.5 * Aeolion::Math::Cross(panel.Corners[2] - panel.Corners[0],
                                                    panel.Corners[3] - panel.Corners[1])
                                   .Norm();
            // Push the control point onto the true sphere surface rather than
            // the chord of the panel, so the comparison below is against the
            // analytic solution at the same place.
            panel.ControlPoint = centroid.Normalized() * radius;
            panels.push_back(panel);
        }
    }
    return panels;
}

void TestSphereMatchesPotentialFlow() {
    const double radius = 1.0;
    const double freestream = 1.0;
    const int polarPanels = 24;
    const int azimuthPanels = 24;

    const std::vector<SourcePanel> panels = BuildSphere(radius, polarPanels, azimuthPanels);
    const int n = static_cast<int>(panels.size());
    const Vec3 Vinf(freestream, 0, 0);

    // Influence matrix and the tangency right-hand side.
    std::vector<double> storage(static_cast<std::size_t>(n) * n, 0.0);
    Aeolion::Solver::DenseMatrixView matrix(storage.data(), n, n);
    std::vector<double> rhs(static_cast<std::size_t>(n), 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j)
            matrix[i, j] = Aeolion::Solver::SourceNormalInfluence(panels[static_cast<std::size_t>(i)],
                                                                  panels[static_cast<std::size_t>(j)]);
        rhs[static_cast<std::size_t>(i)] =
            -Aeolion::Math::Dot(Vinf, panels[static_cast<std::size_t>(i)].Normal);
    }

    const auto factorization = Aeolion::Solver::LuFactorize(matrix);
    CHECK(!factorization.NearSingular, "the sphere influence matrix should factorize cleanly");
    const std::vector<double> sigma = Aeolion::Solver::LuSolve(factorization, rhs);

    // Surface velocity, and how well tangency actually came out.
    double worstTangency = 0.0;
    double worstSpeedError = 0.0;
    for (int i = 0; i < n; ++i) {
        const SourcePanel& panel = panels[static_cast<std::size_t>(i)];
        // The panel's own contribution at its centroid is purely normal: the
        // in-plane part cancels by symmetry.
        Vec3 velocity = Vinf + panel.Normal * (sigma[static_cast<std::size_t>(i)] * 0.5);
        for (int j = 0; j < n; ++j) {
            if (j == i) continue;
            velocity = velocity + Aeolion::Solver::SourcePanelVelocity(panel.ControlPoint,
                                                                       panels[static_cast<std::size_t>(j)]) *
                                      sigma[static_cast<std::size_t>(j)];
        }

        const double normalComponent = Aeolion::Math::Dot(velocity, panel.Normal);
        worstTangency = std::max(worstTangency, std::fabs(normalComponent));

        // Exact potential flow over a sphere: surface speed = 3/2 U sin(theta).
        const double cosTheta = panel.ControlPoint.x / radius;
        const double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
        const double exact = 1.5 * freestream * sinTheta;

        const double tangentialSpeed = (velocity - panel.Normal * normalComponent).Norm();
        // Skip the two polar caps, where the quads degenerate to triangles
        // and the exact speed goes to zero (relative error is meaningless).
        if (sinTheta < 0.2) continue;
        worstSpeedError = std::max(worstSpeedError, std::fabs(tangentialSpeed - exact) / exact);
    }

    std::cout << "sphere: " << n << " panels, worst |V.n| = " << worstTangency
              << ", worst surface-speed error = " << worstSpeedError * 100.0 << "%\n";

    CHECK(worstTangency < 1e-8, "the solve should enforce tangency to round-off, got " +
                                    std::to_string(worstTangency));
    CHECK(worstSpeedError < 0.03, "surface speed should match 3/2 U sin(theta) within 3%, worst was " +
                                      std::to_string(worstSpeedError * 100.0) + "%");
}

} // namespace

int main() {
    TestFarFieldLooksLikeAPointSource();
    TestSheetJump();
    TestSphereMatchesPotentialFlow();

    if (failures == 0) { std::cout << "PASS: TestSourcePanel\n"; return 0; }
    std::cerr << failures << " check(s) failed in TestSourcePanel\n";
    return 1;
}
