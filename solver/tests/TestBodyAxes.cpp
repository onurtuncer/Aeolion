// TestBodyAxes.cpp -- pins the solver-to-contract-frame conversion and the
// body-axis rate derivatives (Solver/BodyAxes.h) against textbook values,
// WITH THEIR SIGNS.
//
// This suite exists because of two defects that a casual check passes.
// Both were live in the DAVE-ML aero-map generator and both are recorded
// in BodyAxes.h; the tests below are the regression guard for each.
//
//   1. A rate derivative does not follow the wrench frame rule. Applying
//      "x and z flip" to Cl_p negates it, giving +0.45 for a wing whose
//      textbook roll damping is -0.45. The MAGNITUDE is right, so any
//      check that tests |Cl_p| ~ 0.45 -- or a range like (-0.8, -0.2)
//      applied to the wrong frame -- passes an aircraft with roll
//      ANTI-damping. Only a signed check catches it, and the sign is the
//      whole physical content: damping opposes the rate.
//
//   2. SolveViscousCoupled leaves its SolveResult's coefficient members
//      at zero and reports dimensional forces instead. Reading res.Base.CL
//      yields a fully converged sweep of exactly zero forces.
//      TestCoupledCoefficientsAreRealNumbers is the guard: a lifting
//      configuration must produce lift.
//
// The textbook anchors used here:
//   - roll damping of a plain rectangular wing sits near Cl_p = -0.45
//     (the same anchor TestSolverCore uses in the solver frame);
//   - lifting-line lift slope CL_alpha = 2*pi/(1 + 2/AR);
//   - an unswept wing at small incidence has essentially no yaw-due-to-
//     roll-rate, and its yaw damping is negative but small;
//   - a 180-degree rotation is an involution, and it must leave y alone.
#include "Aeolion/Solver/BodyAxes.h"
#include "Aeolion/Solver/Solver.h"
#include "Aeolion/Solver/ViscousCoupling.h"

#include <cmath>
#include <iostream>
#include <numbers>
#include <vector>

using namespace Aeolion::Solver;
namespace Math = Aeolion::Math;

static int failures = 0;
#define CHECK(cond, msg)                                                                           \
    do {                                                                                           \
        if (!(cond)) {                                                                             \
            std::cerr << "FAIL: " << msg << "\n";                                                  \
            ++failures;                                                                            \
        }                                                                                          \
    } while (0)

namespace {

constexpr double Span = 6.0;
constexpr double Chord = 1.0;
constexpr double Vinf = 25.0;
constexpr double Rho = 1.225;

struct Fixture {
    std::vector<Panel> Panels;
    ReferenceGeometry Ref;
    PreparedSystem Prepared;
};

/** A plain AR = 6 rectangular unswept wing -- the textbook case. */
Fixture MakeWing() {
    WingParams w;
    w.Span = Span;
    w.RootChord = Chord;
    w.TipChord = Chord;
    w.SweepQuarterChordDeg = 0.0;
    w.DihedralDeg = 0.0;
    w.TwistTipDeg = 0.0;
    w.NPanelsSemiSpan = 20;

    Fixture f;
    f.Panels = BuildWing(w);
    for (Panel& p : f.Panels) p.Surface = "wing";
    double area = 0.0;
    for (const Panel& p : f.Panels) area += p.PlanformArea;
    f.Ref.Area = area;
    f.Ref.Span = w.Span;
    f.Ref.Chord = area / w.Span;
    f.Prepared = Prepare(PanelSystem{f.Panels, {}}, 50.0 * w.Span);
    return f;
}

FreestreamConditions Conditions(double alphaDeg) {
    FreestreamConditions fc;
    fc.Vinf = Vinf;
    fc.rho = Rho;
    fc.alphaDeg = alphaDeg;
    return fc;
}

// --- 1. the frame conversion itself -----------------------------------------

void TestFrameConversionIsAProperInvolution() {
    const Vec3 v(1.0, 2.0, 3.0);
    const Vec3 once = FrdFromSolver(v);
    const Vec3 twice = FrdFromSolver(once);

    CHECK(std::fabs(once.x + v.x) < 1e-15, "x must flip, got " << once.x);
    CHECK(std::fabs(once.y - v.y) < 1e-15, "y must be invariant, got " << once.y);
    CHECK(std::fabs(once.z + v.z) < 1e-15, "z must flip, got " << once.z);

    CHECK(std::fabs(twice.x - v.x) < 1e-15 && std::fabs(twice.y - v.y) < 1e-15 &&
              std::fabs(twice.z - v.z) < 1e-15,
          "the conversion must be an involution: applying it twice is the identity");

    // A PROPER rotation preserves handedness, which is what lets a
    // rotation axis (a pseudovector) be converted like an ordinary vector
    // and keep its right-hand rule -- the property every deflection sign
    // in the flight model rests on. Check it as a determinant: the images
    // of the basis triad must still form a right-handed set.
    const Vec3 ex = FrdFromSolver(Vec3(1, 0, 0));
    const Vec3 ey = FrdFromSolver(Vec3(0, 1, 0));
    const Vec3 ez = FrdFromSolver(Vec3(0, 0, 1));
    const double det = Dot(Cross(ex, ey), ez);
    CHECK(std::fabs(det - 1.0) < 1e-15,
          "the conversion must be a proper rotation (det = +1), got " << det);
}

// --- 2. rate derivatives, with signs -----------------------------------------

void TestRollDampingIsNegativeInBodyAxes() {
    const Fixture f = MakeWing();
    const BodyAxisRateDerivatives d =
        ComputeBodyAxisRateDerivatives(f.Prepared, Conditions(2.0), f.Ref);

    std::cout << "AR=6 body axes: Clp=" << d.Clp << " Cnr=" << d.Cnr << " Cmq=" << d.Cmq
              << " Cnp=" << d.Cnp << " CYp=" << d.CYp << "\n";

    // THE regression guard. A wing must damp roll: the sign is the
    // physics, and the failure mode this catches had the right magnitude.
    CHECK(d.Clp < 0.0, "roll damping must be NEGATIVE in body axes -- a positive Cl_p is roll "
                       "ANTI-damping, i.e. an aircraft that diverges in roll. Got "
                           << d.Clp);
    CHECK(d.Clp < -0.2 && d.Clp > -0.8,
          "AR=6 roll damping should sit near the textbook Cl_p = -0.45, got " << d.Clp);

    // Yaw damping likewise opposes yaw rate.
    CHECK(d.Cnr <= 0.0, "yaw damping must not be positive, got " << d.Cnr);

    // An unswept wing at small incidence generates almost no
    // yaw-due-to-roll-rate; a large value means an axis has been crossed.
    CHECK(std::fabs(d.Cnp) < 0.15,
          "an unswept wing at 2 deg should have small Cn_p, got " << d.Cnp);
}

void TestRateDerivativeFrameRuleIsNotTheWrenchRule() {
    // The identity that the aero-map bug violated. Under the 180-degree
    // rotation about y, both the response and the rate flip for the
    // roll/yaw pairs, so their derivatives are INVARIANT -- they must
    // equal the solver frame's own reduced-rate values EXACTLY, not up to
    // a tolerance and not up to a sign.
    const Fixture f = MakeWing();
    const FreestreamConditions fc = Conditions(2.0);

    const BodyAxisRateDerivatives d = ComputeBodyAxisRateDerivatives(f.Prepared, fc, f.Ref);
    const StabilityDerivatives s =
        ComputeDerivatives(f.Panels, fc, f.Ref, 50.0 * Span);

    const auto identical = [](double a, double b) {
        return std::fabs(a - b) <= 1e-9 * std::fabs(b) + 1e-12;
    };

    CHECK(identical(d.Clp, s.Croll_p_nd),
          "Cl_p must be frame INVARIANT (moment and rate both flip), body=" << d.Clp
              << " solver=" << s.Croll_p_nd);
    CHECK(identical(d.Cnr, s.Cn_r_nd),
          "Cn_r must be frame INVARIANT, body=" << d.Cnr << " solver=" << s.Cn_r_nd);
    CHECK(identical(d.Cnp, s.Cn_p_nd),
          "Cn_p must be frame INVARIANT, body=" << d.Cnp << " solver=" << s.Cn_p_nd);
    CHECK(identical(d.Clr, s.Croll_r_nd),
          "Cl_r must be frame INVARIANT, body=" << d.Clr << " solver=" << s.Croll_r_nd);
    CHECK(identical(d.Cmq, s.Cm_q_nd),
          "Cm_q must be frame INVARIANT (both about y), body=" << d.Cmq
              << " solver=" << s.Cm_q_nd);

    // CZ_q, by contrast, DOES flip: CZ = -CL and q is invariant.
    CHECK(identical(d.CZq, -s.CL_q_nd),
          "CZ_q must be the NEGATIVE of CL_q (CZ flips, q does not), body=" << d.CZq
              << " solver CL_q=" << s.CL_q_nd);
}

void TestReducedRateScalingIsNotInverted() {
    // The (2V/length)^2 trap, restated in body axes: the reduced-rate
    // derivative is the per-rad/s one DIVIDED by length/(2V). Inverting it
    // is a factor of ~2000 here, which reads as an aircraft with no
    // damping at all rather than as an obviously broken number.
    const Fixture f = MakeWing();
    const FreestreamConditions fc = Conditions(2.0);
    const BodyAxisRateDerivatives d = ComputeBodyAxisRateDerivatives(f.Prepared, fc, f.Ref);
    const StabilityDerivatives s = ComputeDerivatives(f.Panels, fc, f.Ref, 50.0 * Span);

    // The exact identity, not a magnitude heuristic: the reduced-rate
    // derivative is the per-rad/s one DIVIDED by length/(2V). Writing the
    // multiplication instead scales by (b/2V)^2 = 0.0144 here -- Cl_p
    // would read -0.0065, which looks like an aircraft with almost no
    // roll damping rather than like an obviously broken number.
    const double reduceFactor = Span / (2.0 * Vinf);
    const double expected = s.Croll_p / reduceFactor;
    CHECK(std::fabs(d.Clp - expected) <= 1e-9 * std::fabs(expected) + 1e-12,
          "Cl_p must equal Croll_p / (b/2V), got " << d.Clp << " against " << expected);

    const double invertedForm = s.Croll_p * reduceFactor;
    CHECK(std::fabs(d.Clp - invertedForm) > 0.1 * std::fabs(d.Clp),
          "Cl_p must not equal the INVERTED form Croll_p * (b/2V) = " << invertedForm);
}

// --- 3. the coefficient assembly ---------------------------------------------

void TestCoupledCoefficientsAreRealNumbers() {
    // The second recorded defect: SolveViscousCoupled leaves its
    // SolveResult coefficient members zero, so a consumer reading
    // res.Base.CL gets a converged sweep of exactly nothing. A lifting
    // configuration must produce lift, and in FRD lift is NEGATIVE CZ
    // (z points down).
    const Fixture f = MakeWing();
    const FreestreamConditions fc = Conditions(4.0);

    std::vector<StripSection> strips;
    strips.reserve(f.Panels.size());
    for (const Panel& p : f.Panels) {
        StripSection s;
        s.ChordDir = Vec3(1.0, 0.0, 0.0);
        s.LiftDir = Vec3(0.0, 0.0, 1.0);
        s.Chord = (p.SpanwiseWidth > 0.0) ? p.PlanformArea / p.SpanwiseWidth : Chord;
        s.Width = p.SpanwiseWidth;
        s.Eta = std::fabs((p.A.y + p.B.y) * 0.5) / (0.5 * Span);
        s.Alpha0Deg = 0.0; // flat plate: the fixture is uncambered
        strips.push_back(s);
    }

    ViscousCouplingOptions options;
    options.Relaxation = 0.05;
    options.AndersonDepth = 0;
    options.MaxIterations = 1000;

    const ViscousCoupledResult res = SolveViscousCoupled(
        f.Panels, strips, fc, f.Ref, 50.0 * Span, AnalyticSectionModel{}, options);
    CHECK(res.Converged, "the attached coupled solve should converge, residual " << res.MaxResidual);

    const double q = 0.5 * Rho * Vinf * Vinf;
    const BodyAxisCoefficients c = BodyAxisFromCoupled(res, q, f.Ref);
    std::cout << "AR=6 coupled at 4 deg: CX=" << c.CX << " CZ=" << c.CZ << " Cm=" << c.Cm << "\n";

    CHECK(std::fabs(c.CZ) > 1e-6,
          "a lifting wing must produce a nonzero CZ -- exactly zero means the coefficient "
          "members of the coupled result were read instead of the dimensional forces");

    // Lift acts UP, which is -z in FRD.
    CHECK(c.CZ < 0.0, "a wing at +4 deg must lift, i.e. CZ < 0 in FRD, got " << c.CZ);

    // Against lifting-line theory: CL = 2*pi*alpha/(1 + 2/AR). The
    // coupled solve carries section drag too, so a loose band is right;
    // the point is order and sign, not a third digit.
    const double AR = Span * Span / f.Ref.Area;
    const double alphaRad = Math::DegToRad(4.0);
    const double clTheory = 2.0 * std::numbers::pi * alphaRad / (1.0 + 2.0 / AR);
    const double clBody = -c.CZ;
    CHECK(std::fabs(clBody - clTheory) < 0.25 * clTheory,
          "coupled CL should track lifting-line theory within 25%, got " << clBody << " against "
                                                                        << clTheory);

    // NOT a check that CX < 0. The body-axis AXIAL force of a lifting
    // surface is positive here, and correctly so: CX = -CD cos(a) + CL
    // sin(a), and a thin lifting surface's leading-edge suction outweighs
    // its profile drag at moderate incidence. (The aero map shows CX
    // crossing zero near 5 degrees for exactly this reason.) Testing the
    // axial force for a drag-like sign would pin an artifact.
    //
    // What must be positive is the STREAMWISE drag -- the force resolved
    // along the free-stream direction, which in FRD at beta = 0 is
    // u = (cos a, 0, sin a). Drag opposes it, so CD = -(C . u).
    const double ca = std::cos(alphaRad), sa = std::sin(alphaRad);
    const double cd = -(c.CX * ca + c.CZ * sa);
    std::cout << "  streamwise CD=" << cd << " (axial CX=" << c.CX
              << " is positive by leading-edge suction)\n";
    CHECK(cd > 0.0, "streamwise drag must be positive, got " << cd);
    CHECK(cd < 0.2, "and physically small on an attached AR=6 wing, got " << cd);
}

void TestSymmetryAtZeroSideslip() {
    // A symmetric configuration at zero sideslip carries no lateral
    // wrench. This catches a y-axis sign slip, which the involution test
    // cannot: y is invariant, so an error there survives it.
    const Fixture f = MakeWing();
    const SolveResult r = SolveWithSystem(f.Prepared, Conditions(6.0), f.Ref);

    CHECK(std::fabs(r.CY) < 1e-9, "symmetric wing at beta=0 must carry no side force, got " << r.CY);
    CHECK(std::fabs(r.Croll) < 1e-9, "and no rolling moment, got " << r.Croll);
    CHECK(std::fabs(r.Cn) < 1e-9, "and no yawing moment, got " << r.Cn);
}

} // namespace

int main() {
    TestFrameConversionIsAProperInvolution();
    TestRollDampingIsNegativeInBodyAxes();
    TestRateDerivativeFrameRuleIsNotTheWrenchRule();
    TestReducedRateScalingIsNotInverted();
    TestCoupledCoefficientsAreRealNumbers();
    TestSymmetryAtZeroSideslip();

    if (failures == 0) {
        std::cout << "PASS: TestBodyAxes\n";
        return 0;
    }
    std::cerr << failures << " check(s) failed in TestBodyAxes\n";
    return 1;
}
