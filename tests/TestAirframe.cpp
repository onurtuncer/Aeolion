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

std::string FixturePath180() {
    return std::string(AEOLION_TEST_DATA_DIR) + "/AeolionGeometryHandoff-1.8.0.json";
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

    // Scaling the coupled ratio back up by the block mismatch should land
    // near the wing-only figure. This is an order-of-magnitude argument, not
    // an identity -- the two matrices are not the same matrix -- so the band
    // is a factor of three. It is enough to distinguish "the metric is
    // measuring a unit mismatch" from "the matrix is a hundred times worse",
    // which is the question that mattered.
    const double rescaled =
        coupled.Factorization.MinPivotRatio * (vortexColumnMax / sourceColumnMax);
    CHECK(rescaled > wingOnly.Factorization.MinPivotRatio / 3.0 &&
              rescaled < wingOnly.Factorization.MinPivotRatio * 3.0,
          "the coupled pivot ratio should be explained by the block magnitude mismatch: rescaled " +
              std::to_string(rescaled) + " against wing-only " +
              std::to_string(wingOnly.Factorization.MinPivotRatio));

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
    // Three lattices isolate what each mechanism does:
    //   whole  -- untrimmed, wing run through the body (geometrically wrong,
    //             but it is the lift a full-span wing carries)
    //   cut    -- trimmed with NO carry-through: the root load simply lost
    //   joined -- trimmed WITH the carry-through vortex
    // Trimming must cost lift, and the carry-through must give it back. If
    // it did not, the carry-through would be decoration.
    PB::LatticeOptions untrimmedOptions;
    untrimmedOptions.TrimWingAtBody = false;
    const auto untrimmed = Aeolion::Solver::Prepare(
        Aeolion::Solver::PanelSystem{PB::LatticeBuilder(contract, untrimmedOptions).Build(), {}}, trail);

    PB::LatticeOptions noCarryOptions;
    noCarryOptions.CarryThroughLift = false;
    const auto severed = Aeolion::Solver::Prepare(
        Aeolion::Solver::PanelSystem{PB::LatticeBuilder(contract, noCarryOptions).Build(), {}}, trail);

    std::cout << "  alpha    whole      cut    joined   joined+body   carry-through   body\n";
    double previousDCm = 0.0;
    for (double alphaDeg : {0.0, 4.0, 8.0}) {
        fc.alphaDeg = alphaDeg;
        const auto whole = Aeolion::Solver::SolveWithSystem(untrimmed, fc, ref);
        const auto cut = Aeolion::Solver::SolveWithSystem(severed, fc, ref);
        const auto joined = Aeolion::Solver::SolveWithSystem(wingOnly, fc, ref);
        const auto both = Aeolion::Solver::SolveWithSystem(coupled, fc, ref);

        const double restored = (joined.CL - cut.CL) / (whole.CL - cut.CL);
        const double bodyEffect = (both.CL - joined.CL) / joined.CL;
        std::printf("  %5.1f  %8.5f %8.5f  %8.5f  %11.5f  %12.1f%% %+6.1f%%\n", alphaDeg, whole.CL, cut.CL,
                    joined.CL, both.CL, 100.0 * restored, 100.0 * bodyEffect);

        CHECK(std::isfinite(both.CL) && std::isfinite(both.Cm), "the coupled solve must be finite");
        CHECK(cut.CL < whole.CL, "trimming without carry-through must cost lift at alpha=" +
                                     std::to_string(alphaDeg));

        // The carry-through restores the bound circulation the fuselage
        // interrupts, so a trimmed wing carrying it should reproduce the
        // full-span wing almost exactly. It is the same vorticity spanning
        // the same distance -- only the tangency points moved outboard, off
        // the part of the wing that is inside the body.
        CHECK(std::fabs(joined.CL - whole.CL) < 0.02 * whole.CL,
              "the carry-through should restore the full-span lift, got " + std::to_string(joined.CL) +
                  " against " + std::to_string(whole.CL) + " at alpha=" + std::to_string(alphaDeg));
        CHECK(restored > 0.95,
              "the carry-through should account for essentially all the trimmed load, restored " +
                  std::to_string(100.0 * restored) + "%");

        // On top of that, the body's own displacement adds a modest
        // increment -- blockage and upwash over the exposed wing. Modest is
        // the point: the body's frontal area is a few percent of the wing's,
        // so a large swing would mean the coupling is wrong, not strong.
        CHECK(bodyEffect > 0.0 && bodyEffect < 0.20,
              "the body should add a modest increment on top of the carry-through, got " +
                  std::to_string(100.0 * bodyEffect) + "% at alpha=" + std::to_string(alphaDeg));

        if (alphaDeg > 0.0)
            CHECK(std::fabs(both.Cm - joined.Cm) > std::fabs(previousDCm),
                  "the body's pitching contribution should grow with incidence");
        previousDCm = both.Cm - joined.Cm;
    }
}

// The viewer's Application keeps ONE LatticeBuilder alive for the whole
// session and, on every frame a control-surface slider changes,
// ClearDeflections()s + Deflect()s + Build()s it again, re-solving the SAME
// cached fuselage/duct source panels each time (see
// viewer/src/Core/Application.cpp::Resolve()). TestPanelBuilder.cpp already
// proves a single Deflect().Build() moves the WING lattice; this checks the
// pattern the viewer actually runs -- repeated rebuilds of one builder,
// solved coupled with the body and the duct -- because that repetition is
// exactly what was untested before the viewer started doing it every frame.
void TestControlDeflectionMovesCoupledSolve() {
    const auto contract = Aeolion::Geometry::LoadHandoff(FixturePath180());

    std::size_t aileronIndex = contract.ControlSurfaces.size();
    for (std::size_t i = 0; i < contract.ControlSurfaces.size(); ++i)
        if (contract.ControlSurfaces[i].Name == "aileron") { aileronIndex = i; break; }
    CHECK(aileronIndex < contract.ControlSurfaces.size(), "the 1.8.0 fixture should carry an aileron");

    PB::LatticeBuilder builder(contract);
    const auto body = builder.BuildBody();
    const auto duct = builder.BuildDuct();
    CHECK(!body.empty(), "the 1.8.0 fixture should produce body panels");
    CHECK(!duct.empty(), "the 1.8.0 fixture should produce duct panels");
    std::vector<Aeolion::Lattice::SourcePanel> sources = body;
    sources.insert(sources.end(), duct.begin(), duct.end());

    const double area = builder.GrossPlanformArea();
    Aeolion::Solver::ReferenceGeometry ref;
    ref.Area = area; ref.Span = contract.Span; ref.Chord = area / contract.Span;
    Aeolion::Solver::FreestreamConditions fc;
    fc.Vinf = 25.0; fc.rho = 1.225; fc.alphaDeg = 4.0;
    const double trail = 50.0 * contract.Span;

    const auto solve = [&](const std::vector<Aeolion::Lattice::Panel>& wing) {
        const auto prepared =
            Aeolion::Solver::Prepare(Aeolion::Solver::PanelSystem{wing, sources}, trail);
        return Aeolion::Solver::SolveWithSystem(prepared, fc, ref);
    };

    // Three states on the SAME builder instance, in sequence -- exactly what
    // Application::Resolve() runs across consecutive frames as a slider moves.
    const auto neutral = solve(builder.Build());

    builder.Deflect(PB::Antisymmetric(aileronIndex, 10.0));
    const auto deflected = solve(builder.Build());

    std::cout << "control deflection: neutral CL=" << neutral.CL << " Croll=" << neutral.Croll
              << "  deflected CL=" << deflected.CL << " Croll=" << deflected.Croll << "\n";

    // An antisymmetric aileron command must roll the coupled airframe...
    CHECK(std::fabs(deflected.Croll - neutral.Croll) > 1e-3,
          "deflecting the aileron through the cached LatticeBuilder must move the coupled Croll, got neutral=" +
              std::to_string(neutral.Croll) + " deflected=" + std::to_string(deflected.Croll));
    // ...without changing how much lift the coupled airframe makes overall.
    CHECK(std::fabs(deflected.CL - neutral.CL) < 0.02 * std::fabs(neutral.CL) + 1e-6,
          "an antisymmetric aileron deflection should leave total coupled CL essentially unchanged, got neutral=" +
              std::to_string(neutral.CL) + " deflected=" + std::to_string(deflected.CL));

    // ClearDeflections() on the same builder must reproduce the ORIGINAL
    // neutral answer exactly -- state must not leak between Build() calls,
    // since Application::Resolve() calls ClearDeflections() before every
    // Deflect() precisely to prevent that.
    builder.ClearDeflections();
    const auto cleared = solve(builder.Build());
    CHECK(std::fabs(cleared.CL - neutral.CL) < 1e-12,
          "ClearDeflections() must exactly reproduce the neutral CL");
    CHECK(std::fabs(cleared.Croll - neutral.Croll) < 1e-12,
          "ClearDeflections() must exactly reproduce the neutral Croll");
}

} // namespace

int main() {
    TestCoupledAirframe();
    TestControlDeflectionMovesCoupledSolve();

    if (failures == 0) { std::cout << "PASS: TestAirframe\n"; return 0; }
    std::cerr << failures << " check(s) failed in TestAirframe\n";
    return 1;
}
