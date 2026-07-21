// TestAirframe.cpp -- the whole chain on the real airframe: parse the 1.4.0
// handoff, build the wing lattice and the fuselage panels, and solve them
// coupled. This is the first place every piece meets, so it is also where
// integration problems that no single-component test can see should show up.
#include "Aeolion/PanelBuilder/PanelBuilder.h"
#include "Aeolion/Solver/Solver.h"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <numbers>
#include <string>
#include <vector>

using Aeolion::Math::Vec3;

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } } while (0)

namespace {

namespace PB = Aeolion::PanelBuilder;

std::string FixturePath() {
    return std::string(AEOLION_TEST_DATA_DIR) + "/AeolionGeometryHandoff-1.5.0.json";
}

// The solver frame is x aft and the contract frame x forward, so a
// solver-frame position is negated before asking the contract's own radius
// law (Geometry::RadiusAt).
double BodyRadiusAtSolverX(const Aeolion::Geometry::BodyGeometry& body, double solverX) {
    return Aeolion::Geometry::RadiusAt(body, -solverX);
}

void TestCoupledAirframe() {
    const auto contract = Aeolion::Geometry::LoadHandoff(FixturePath());
    PB::LatticeBuilder builder(contract);
    const auto wing = builder.Build();
    const auto body = builder.BuildBody();

    CHECK(!wing.empty(), "the fixture should produce a wing lattice");
    CHECK(!body.empty(), "the fixture should produce body panels");
    std::cout << "airframe: " << wing.size() << " wing panels, " << body.size() << " body panels\n";

    // --- where the two components actually sit ------------------------------
    double bodyNose = 1e30, bodyTail = -1e30, bodyMaxRadius = 0.0;
    for (const auto& panel : body)
        for (const Vec3& corner : panel.Corners) {
            bodyNose = std::min(bodyNose, corner.x);
            bodyTail = std::max(bodyTail, corner.x);
            bodyMaxRadius = std::max(bodyMaxRadius, std::sqrt(corner.y * corner.y + corner.z * corner.z));
        }
    double wingFore = 1e30, wingAft = -1e30;
    for (const auto& panel : wing) {
        wingFore = std::min({wingFore, panel.A.x, panel.B.x});
        wingAft = std::max({wingAft, panel.A.x, panel.B.x});
    }
    std::cout << "  body x [" << bodyNose << ", " << bodyTail << "], max radius " << bodyMaxRadius << "\n"
              << "  wing bound vortices x [" << wingFore << ", " << wingAft << "]\n";

    // --- INTEGRATION FINDING: the contract states no relative placement -----
    // The planform gives span/chord/sweep/twist and the body gives its own
    // profile, but nothing positions one against the other. Both therefore
    // start at the origin, which puts the wing at the fuselage NOSE. Count
    // how much of the wing ends up buried inside the body, because a lifting
    // panel whose control point is inside a solid volume is not a physically
    // meaningful boundary condition.
    int buried = 0;
    for (const auto& panel : wing) {
        const double radius = std::sqrt(panel.ControlPoint.y * panel.ControlPoint.y +
                                        panel.ControlPoint.z * panel.ControlPoint.z);
        if (radius < BodyRadiusAtSolverX(contract.Body, panel.ControlPoint.x)) ++buried;
    }
    std::cout << "  wing control points inside the body: " << buried << " of " << wing.size() << "\n";

    // --- the coupled solve --------------------------------------------------
    // GROSS planform area, from the contract's stations -- not the trimmed
    // panel sum. CL is conventionally referred to the wing extended through
    // the fuselage to the centreline, so letting the trim shrink the
    // reference would rescale every coefficient by whatever was cut.
    const double area = builder.GrossPlanformArea();
    Aeolion::Solver::ReferenceGeometry ref;
    ref.Area = area; ref.Span = contract.Span; ref.Chord = area / contract.Span;
    std::cout << "  gross planform area " << area << " m^2, trimmed at eta " << builder.TrimEta() << "\n";
    Aeolion::Solver::FreestreamConditions fc;
    fc.Vinf = 25.0; fc.rho = 1.225;
    const double trail = 50.0 * contract.Span;

    const auto wingOnly = Aeolion::Solver::Prepare(Aeolion::Solver::PanelSystem{wing, {}}, trail);
    const auto coupled = Aeolion::Solver::Prepare(Aeolion::Solver::PanelSystem{wing, body}, trail);

    CHECK(!coupled.Factorization.NearSingular, "the coupled airframe matrix must factorize");

    // CONDITIONING, and why the obvious metric lies here. Prepare() reports
    // MinPivotRatio 0.86 for the wing alone and 0.0064 coupled, which looks
    // like the coupling destroys the system. It does not. That ratio divides
    // every pivot by the largest entry anywhere in A, which is only
    // meaningful if the entries are commensurable -- and they are not: a
    // horseshoe's normal influence per unit circulation carries units of
    // 1/length, while a source panel's per unit strength is dimensionless.
    //
    // Three hypotheses were tested and rejected before landing on that.
    // Trimming the wing panels buried in the body moved the ratio to 0.0062
    // (no help). Body-alone scores 0.96 and the body's lateral surface 0.99,
    // so neither component is individually ill-conditioned. Moving the body
    // out of the wing's wake recovered only a factor of two. What does
    // explain it is the block magnitudes, measured below, and the residual,
    // which is the honest test and is at machine precision.
    double vortexColumnMax = 0.0, sourceColumnMax = 0.0;
    const Aeolion::Solver::PanelSystem system{wing, body};
    const int n = system.UnknownCount();
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            const double a = std::fabs(system.NormalInfluence(i, j, trail));
            if (system.IsVortex(j)) vortexColumnMax = std::max(vortexColumnMax, a);
            else sourceColumnMax = std::max(sourceColumnMax, a);
        }
    std::cout << "  pivot ratio: wing-only " << wingOnly.Factorization.MinPivotRatio << ", coupled "
              << coupled.Factorization.MinPivotRatio << "\n"
              << "  largest |A|: vortex columns " << vortexColumnMax << ", source columns "
              << sourceColumnMax << " (ratio " << vortexColumnMax / sourceColumnMax << ")\n";

    // 0.86 / 149 is 0.0058, which is the "collapsed" ratio almost exactly --
    // the metric is measuring the unit mismatch, not the matrix.
    CHECK(std::fabs(wingOnly.Factorization.MinPivotRatio / (vortexColumnMax / sourceColumnMax) -
                    coupled.Factorization.MinPivotRatio) < 0.3 * coupled.Factorization.MinPivotRatio,
          "the coupled pivot ratio should be explained by the block magnitude mismatch");

    // The real measure: does the solution satisfy its own equations?
    fc.alphaDeg = 6.0;
    const double alphaRad = fc.alphaDeg * std::numbers::pi / 180.0;
    const Vec3 Vinf(fc.Vinf * std::cos(alphaRad), 0.0, fc.Vinf * std::sin(alphaRad));
    std::vector<double> rhs(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i)
        rhs[static_cast<std::size_t>(i)] =
            system.PrescribedNormalVelocity(i) - Aeolion::Math::Dot(Vinf, system.Normal(i));

    const std::vector<double> x = Aeolion::Solver::LuSolve(coupled.Factorization, rhs);
    double residual = 0.0, scale = 0.0;
    for (int i = 0; i < n; ++i) {
        double row = 0.0;
        for (int j = 0; j < n; ++j) row += system.NormalInfluence(i, j, trail) * x[static_cast<std::size_t>(j)];
        residual = std::max(residual, std::fabs(row - rhs[static_cast<std::size_t>(i)]));
        scale = std::max(scale, std::fabs(rhs[static_cast<std::size_t>(i)]));
    }
    std::cout << "  coupled solve relative residual " << residual / scale << "\n";
    CHECK(residual / scale < 1e-10,
          "the coupled solve must satisfy its own equations; relative residual " +
              std::to_string(residual / scale));

    // --- carryover: does the body pick up the load the trim removed? -------
    // Trimming deletes the wing's root panels, so a trimmed wing ALONE must
    // lose lift. The body is what carries that load in reality, so the
    // coupled answer should come back to roughly what the untrimmed wing
    // gave. That round trip is the real test of the wing-body coupling: it
    // is not a comparison against the builder's own arithmetic, and neither
    // half of it passes on its own.
    PB::LatticeOptions untrimmedOptions;
    untrimmedOptions.TrimWingAtBody = false;
    const auto untrimmedWing = PB::LatticeBuilder(contract, untrimmedOptions).Build();
    const auto untrimmed =
        Aeolion::Solver::Prepare(Aeolion::Solver::PanelSystem{untrimmedWing, {}}, trail);

    std::cout << "  alpha  CL untrimmed  CL trimmed  CL trimmed+body   recovery\n";
    double previousDCm = 0.0;
    for (double alphaDeg : {0.0, 4.0, 8.0}) {
        fc.alphaDeg = alphaDeg;
        const auto whole = Aeolion::Solver::SolveWithSystem(untrimmed, fc, ref);
        const auto cut = Aeolion::Solver::SolveWithSystem(wingOnly, fc, ref);
        const auto both = Aeolion::Solver::SolveWithSystem(coupled, fc, ref);

        const double recovery = (both.CL - cut.CL) / (whole.CL - cut.CL);
        std::printf("  %5.1f  %11.5f  %10.5f  %15.5f  %8.1f%%\n", alphaDeg, whole.CL, cut.CL, both.CL,
                    100.0 * recovery);

        CHECK(std::isfinite(both.CL) && std::isfinite(both.Cm), "the coupled solve must be finite");
        CHECK(cut.CL < whole.CL, "trimming the root must cost lift at alpha=" + std::to_string(alphaDeg));

        // The body recovers most but NOT all of it -- measured 64% to 75%
        // over this range, rising with incidence. The shortfall is
        // structural, not a tuning problem: a source distribution cannot
        // carry circulation, and lift IS circulation. Sources enforce
        // tangency and so produce the upwash that loads the exposed wing
        // harder, which is where the recovered part comes from, but the
        // bound vorticity that the real wing carries THROUGH the fuselage
        // has nowhere to live in a source-only body.
        //
        // Closing the remaining quarter needs a carry-through vortex
        // spanning the body at the strength of the innermost exposed strip
        // -- the classical Pitts/Nielsen/Kaattari treatment -- which is a
        // constraint between unknowns rather than another panel. Until then
        // this band records what the model actually does rather than what
        // slender-body theory says it should.
        CHECK(recovery > 0.55 && recovery < 0.90,
              "the body should carry most of the trimmed root load, recovered " +
                  std::to_string(100.0 * recovery) + "% at alpha=" + std::to_string(alphaDeg));

        if (alphaDeg > 0.0)
            CHECK(std::fabs(both.Cm - cut.Cm) > std::fabs(previousDCm),
                  "the body's pitching contribution should grow with incidence");
        previousDCm = both.Cm - cut.Cm;
    }
}

} // namespace

int main() {
    TestCoupledAirframe();

    if (failures == 0) { std::cout << "PASS: TestAirframe\n"; return 0; }
    std::cerr << failures << " check(s) failed in TestAirframe\n";
    return 1;
}
