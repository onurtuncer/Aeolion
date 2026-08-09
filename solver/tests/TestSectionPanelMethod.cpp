// TestSectionPanelMethod.cpp -- validates the 2-D Hess-Smith section solve
// and the stagnation point it exists to find.
//
//   - Circular cylinder. Fed a circle instead of an airfoil, the method must
//     reproduce the exact potential-flow answer |V| = 2U sin(theta) and
//     Cp = 1 - 4 sin^2(theta), with zero circulation. This tests the
//     influence coefficients, the assembly and the solve together, against
//     a result that owes nothing to airfoil theory. The Kutta condition is
//     satisfied trivially at zero incidence, since the two panels flanking
//     the start point are mirror images.
//
//   - Symmetric section at zero incidence. The stagnation point must sit
//     exactly on the leading edge and the lift must vanish -- a symmetry
//     the solve either has or does not.
//
//   - Lift slope. A thin symmetric section must give cl close to 2 pi alpha,
//     a little above it for thickness.
//
//   - Stagnation-point movement, against matched asymptotics. Raising
//     incidence must walk the stagnation point aft along the LOWER surface,
//     monotonically, and it must do so on the scale that inner/outer
//     matching predicts:
//
//         s_stag/c  ~  sqrt(2 r_LE/c) * alpha_e.
//
//     The naive reading -- that a round nose behaves like a cylinder in a
//     stream at angle alpha, giving s ~ r_LE * alpha -- is WRONG, by an
//     order of magnitude on a typical section, and it is worth saying why.
//     A cylinder sits in a uniform stream; an airfoil nose sits in the
//     leading-edge singularity of the outer thin-airfoil solution, where
//     u ~ alpha sqrt(c/s). Setting that against the parabolic nose's own
//     surface speed sqrt(s/r_LE) gives s ~ alpha sqrt(r_LE c), so the offset
//     scales with the SQUARE ROOT of the nose radius. The test asserts the
//     collapse across a twelvefold range of nose radius, which is the part
//     no fitted constant could fake.
//
//   - Boundary-layer runs. The two runs must start at the stagnation point,
//     between them cover the whole contour, and the upper one must be
//     LONGER than the bare upper surface, because it wraps around the nose
//     from a stagnation point that lies on the lower side.
#include "Aeolion/Geometry/SectionContour.h"
#include "Aeolion/Solver/SectionPanelMethod.h"

#include <cmath>
#include <iostream>
#include <numbers>
#include <string>
#include <vector>

using Aeolion::Geometry::AirfoilSection;
using Aeolion::Geometry::SectionContour;
namespace S = Aeolion::Solver;
namespace M = Aeolion::Math;

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } } while (0)

namespace {

constexpr double Pi = std::numbers::pi;

// A circle of diameter 1 spanning psi in [0,1], traversed clockwise from
// psi = 1 -- the same ordering BuildSectionContour produces, so the panel
// method cannot tell it is not an airfoil.
SectionContour CircleContour(int points) {
    SectionContour contour;
    for (int k = 0; k <= points; ++k) {
        const double phi = 2.0 * Pi * k / points;
        contour.Psi.push_back(0.5 + 0.5 * std::cos(phi));
        contour.Zeta.push_back(-0.5 * std::sin(phi));
    }
    contour.LeadingEdgeRadius = 0.5;
    return contour;
}

// A symmetric CST section: equal and opposite coefficients on the two
// surfaces. The first coefficient sets the nose radius, r_LE/c = A0^2/2.
AirfoilSection SymmetricSection(double scale) {
    AirfoilSection section;
    section.Eta = 0.0;
    section.CoefficientsUpper = {scale, scale, scale, scale};
    section.CoefficientsLower = {-scale, -scale, -scale, -scale};
    return section;
}

void TestCylinderIsExact() {
    const SectionContour contour = CircleContour(240);
    const S::SectionSolution solution = S::SolveSectionContour(contour, 0.0);
    CHECK(solution.Valid, "the cylinder solve must succeed");
    if (!solution.Valid) return;

    CHECK(std::fabs(solution.cl) < 1e-6,
          "a cylinder at zero incidence carries no lift, got cl=" + std::to_string(solution.cl));

    double worstSpeed = 0.0, worstCp = 0.0;
    for (std::size_t k = 0; k < solution.Ue.size(); ++k) {
        // Polar angle from the FRONT stagnation point, which is at psi = 0.
        const double dx = solution.Psi[k] - 0.5, dz = solution.Zeta[k];
        const double theta = std::atan2(std::hypot(dz, 0.0), -dx); // angle from -x
        const double expectedSpeed = 2.0 * std::sin(theta);
        const double expectedCp = 1.0 - 4.0 * std::sin(theta) * std::sin(theta);
        worstSpeed = std::max(worstSpeed, std::fabs(std::fabs(solution.Ue[k]) - expectedSpeed));
        worstCp = std::max(worstCp, std::fabs(solution.Cp[k] - expectedCp));
    }
    CHECK(worstSpeed < 0.01, "cylinder surface speed must be 2 sin(theta), worst error " +
                                 std::to_string(worstSpeed));
    CHECK(worstCp < 0.02, "cylinder Cp must be 1 - 4 sin^2(theta), worst error " + std::to_string(worstCp));

    // The stagnation point sits at the leading edge of the circle.
    CHECK(solution.StagnationFound, "the cylinder's forward stagnation point must be located");
    CHECK(std::fabs(solution.StagnationPsi) < 0.01,
          "the cylinder's stagnation point is at psi = 0, got " + std::to_string(solution.StagnationPsi));

    // Near the stagnation point of a cylinder of radius a, U_e = 2 U s/a, so
    // the strain rate is 2/a = 4 per chord for a = 1/2.
    const double expectedStrain = 4.0;
    CHECK(std::fabs(solution.StagnationStrain - expectedStrain) / expectedStrain < 0.05,
          "cylinder stagnation strain must be 2U/a, got " + std::to_string(solution.StagnationStrain) +
              " vs " + std::to_string(expectedStrain));
}

void TestSymmetricSectionAtZeroIncidence() {
    const SectionContour contour = Aeolion::Geometry::BuildSectionContour(SymmetricSection(0.17), 100);
    const S::SectionSolution solution = S::SolveSectionContour(contour, 0.0);
    CHECK(solution.Valid, "the symmetric section solve must succeed");
    if (!solution.Valid) return;

    CHECK(std::fabs(solution.cl) < 1e-6,
          "a symmetric section at zero incidence carries no lift, got cl=" + std::to_string(solution.cl));
    CHECK(solution.StagnationFound, "the stagnation point must be located");
    CHECK(std::fabs(solution.StagnationPsi) < 1e-3,
          "at zero incidence a symmetric section stagnates ON the leading edge, got psi=" +
              std::to_string(solution.StagnationPsi));
}

void TestLiftSlope() {
    const SectionContour contour = Aeolion::Geometry::BuildSectionContour(SymmetricSection(0.17), 100);
    const double thickness = Aeolion::Geometry::MaxThickness(contour);

    for (double alphaDeg : {2.0, 5.0, 8.0}) {
        const S::SectionSolution solution = S::SolveSectionContour(contour, M::DegToRad(alphaDeg));
        const double thin = 2.0 * Pi * M::DegToRad(alphaDeg);
        const double ratio = solution.cl / thin;
        // Thickness raises the slope above 2 pi by roughly 0.77 t/c.
        CHECK(ratio > 1.0 && ratio < 1.0 + 4.0 * thickness,
              "cl must sit just above thin-airfoil at alpha=" + std::to_string(alphaDeg) +
                  " deg: cl/2*pi*alpha = " + std::to_string(ratio) + " with t/c=" +
                  std::to_string(thickness));
    }
}

void TestStagnationMovesAftWithIncidence() {
    const SectionContour contour = Aeolion::Geometry::BuildSectionContour(SymmetricSection(0.17), 120);
    const double noseRadius = contour.LeadingEdgeRadius;
    CHECK(std::fabs(noseRadius - 0.5 * 0.17 * 0.17) < 1e-12,
          "the CST nose radius must be A0^2/2, got " + std::to_string(noseRadius));

    double previousArc = -1.0;
    for (double alphaDeg : {0.5, 1.0, 2.0, 4.0, 6.0, 8.0}) {
        const double alpha = M::DegToRad(alphaDeg);
        const S::SectionSolution solution = S::SolveSectionContour(contour, alpha);
        const std::string label = "alpha=" + std::to_string(alphaDeg) + " deg";
        CHECK(solution.StagnationFound, "a stagnation point must exist at " + label);
        if (!solution.StagnationFound) continue;

        CHECK(solution.StagnationOnLower,
              "at positive incidence the stagnation point is on the LOWER surface, " + label);
        CHECK(solution.StagnationZeta < 0.0,
              "the stagnation point must lie below the chord line at " + label);

        // Arc from the leading edge, positive going back along the lower
        // surface. The leading edge sits at the contour's midpoint arc.
        const double leadingEdgeArc = 0.5 * solution.Perimeter;
        const double offset = leadingEdgeArc - solution.StagnationArc;
        CHECK(offset > previousArc,
              "the stagnation point must move monotonically aft with incidence at " + label);
        previousArc = offset;

        // Matched asymptotics: s_stag ~ sqrt(2 r_LE) * alpha on a unit chord.
        const double estimate = std::sqrt(2.0 * noseRadius) * alpha;
        if (alphaDeg <= 2.0) {
            const double ratio = offset / estimate;
            CHECK(ratio > 0.85 && ratio < 1.15,
                  "the stagnation offset must follow sqrt(2 r_LE)*alpha at " + label + ": got " +
                      std::to_string(offset) + " vs " + std::to_string(estimate) + " (ratio " +
                      std::to_string(ratio) + ")");
        }

        CHECK(solution.StagnationStrain > 0.0, "the stagnation strain rate must be positive at " + label);
    }
}

// The sharp form of the same statement: the offset must COLLAPSE onto
// sqrt(r_LE)*alpha across nose radii. A wrong scaling law (r_LE*alpha, say)
// cannot survive a twelvefold change in nose radius no matter what constant
// is put in front of it.
void TestStagnationOffsetScalesWithSqrtNoseRadius() {
    const double alpha = M::DegToRad(1.5);
    double firstCollapsed = 0.0;

    for (double scale : {0.10, 0.17, 0.24, 0.34}) {
        const SectionContour contour = Aeolion::Geometry::BuildSectionContour(SymmetricSection(scale), 200);
        const S::SectionSolution solution = S::SolveSectionContour(contour, alpha);
        CHECK(solution.StagnationFound, "a stagnation point must exist for A0=" + std::to_string(scale));
        if (!solution.StagnationFound) continue;

        const double offset = 0.5 * solution.Perimeter - solution.StagnationArc;
        const double collapsed = offset / (std::sqrt(contour.LeadingEdgeRadius) * alpha);

        // The collapsed value is sqrt(2) = 1.414 to leading order.
        CHECK(collapsed > 1.15 && collapsed < 1.65,
              "offset/(sqrt(r_LE)*alpha) must be near sqrt(2) for A0=" + std::to_string(scale) +
                  ", got " + std::to_string(collapsed));

        if (firstCollapsed == 0.0) {
            firstCollapsed = collapsed;
        } else {
            const double spread = std::fabs(collapsed - firstCollapsed) / firstCollapsed;
            CHECK(spread < 0.12, "the collapse must hold across nose radii: A0=" + std::to_string(scale) +
                                     " gives " + std::to_string(collapsed) + " against " +
                                     std::to_string(firstCollapsed));
        }
    }
}

void TestBoundaryLayerRunsCoverTheContour() {
    const SectionContour contour = Aeolion::Geometry::BuildSectionContour(SymmetricSection(0.17), 100);
    const S::SectionSolution solution = S::SolveSectionContour(contour, M::DegToRad(6.0));
    CHECK(solution.StagnationFound, "runs need a stagnation point");
    if (!solution.StagnationFound) return;

    const S::SurfaceRun upper = solution.UpperRun();
    const S::SurfaceRun lower = solution.LowerRun();
    CHECK(upper.Count() > 10 && lower.Count() > 10, "both runs must carry stations");
    CHECK(upper.Count() + lower.Count() == solution.S.size(),
          "between them the two runs must cover every panel exactly once");

    // Both start at the stagnation point and run to a trailing edge, so
    // together they span the whole perimeter.
    const double covered = upper.Length() + lower.Length();
    CHECK(std::fabs(covered - solution.Perimeter) < 0.05 * solution.Perimeter,
          "the two runs together must span the perimeter: " + std::to_string(covered) + " vs " +
              std::to_string(solution.Perimeter));

    // The upper run wraps around the nose from a stagnation point on the
    // lower surface, so it is longer than half the perimeter.
    CHECK(upper.Length() > 0.5 * solution.Perimeter,
          "the upper run must include the wrap around the nose, got " + std::to_string(upper.Length()) +
              " against half-perimeter " + std::to_string(0.5 * solution.Perimeter));

    // Arc must increase and speeds be non-negative on both runs.
    for (const S::SurfaceRun* run : {&upper, &lower}) {
        for (std::size_t k = 1; k < run->Count(); ++k)
            CHECK(run->S[k] > run->S[k - 1], "run arc length must increase monotonically");
        for (double speed : run->Ue) CHECK(speed >= 0.0, "run edge speed must be a magnitude");
    }

    // The suction peak belongs on the upper run, near its start.
    std::size_t peak = 0;
    for (std::size_t k = 1; k < upper.Count(); ++k)
        if (upper.Ue[k] > upper.Ue[peak]) peak = k;
    CHECK(upper.S[peak] < 0.25 * upper.Length(),
          "the suction peak must sit in the first quarter of the upper run");
}

void TestHiemenzInitialThickness() {
    // theta_0 = sqrt(0.075 nu / (dUe/ds)). A 1 m chord at 30 m/s in air with
    // a strain rate of 20 per chord gives a layer of order tens of microns.
    const double theta = S::StagnationMomentumThickness(20.0, 1.0, 30.0, 1.46e-5);
    const double expected = std::sqrt(0.075 * 1.46e-5 / (20.0 * 30.0 / 1.0));
    CHECK(std::fabs(theta - expected) < 1e-12, "Hiemenz initial thickness formula");
    CHECK(theta > 0.0 && theta < 1e-3, "the initial momentum thickness must be small but nonzero, got " +
                                           std::to_string(theta));
    CHECK(S::StagnationMomentumThickness(0.0, 1.0, 30.0, 1.46e-5) == 0.0,
          "a zero strain rate must not produce an infinity");
}

} // namespace

int main() {
    TestCylinderIsExact();
    TestSymmetricSectionAtZeroIncidence();
    TestLiftSlope();
    TestStagnationMovesAftWithIncidence();
    TestStagnationOffsetScalesWithSqrtNoseRadius();
    TestBoundaryLayerRunsCoverTheContour();
    TestHiemenzInitialThickness();

    if (failures == 0) std::cout << "PASS: TestSectionPanelMethod\n";
    return failures == 0 ? 0 : 1;
}
