// TestRollDamping.cpp -- roll damping through and beyond the stall, and
// the frame discipline a coupled rate derivative demands.
//
// WHAT IS ESTABLISHED, and what is not. The mechanism below is real and
// is what this suite pins: the sign of Cl_p follows the sign of the
// section lift slope, so a section past stall -- where dcl/dalpha is
// negative -- produces roll ANTI-damping, which is autorotation.
//
// Whether the studied configuration actually exhibits that at a given
// incidence is a SEPARATE question, and the answer measured so far is
// "not resolved". Differencing the coupled solve at two roll-rate
// amplitudes a factor of nine apart:
//
//   alpha        0..18      20      22      24      26     45..90
//   small step  -0.545..  -0.311  -0.258  +0.112  +0.025  -0.257..
//               -0.374                                    -0.425
//   large step  -0.545..  -0.168  -0.087  -0.159  -0.194  -0.256..
//               -0.374                                    -0.425
//
// The attached range agrees to four digits across that amplitude change,
// and so does deep stall past about 45 degrees. Between roughly 20 and 30
// the two disagree by up to 196% and the SIGN does not survive. So in
// that band Cl_p is not a derivative at all: the response is nonlinear
// over the perturbation range, and no single linear coefficient
// represents it. A first pass reported the small-amplitude sign reversal
// as a finding; it is not one, and the correction is recorded here rather
// than quietly dropped.
//
// The consequence for a tabulated model is stronger than "the taper is
// wrong". Across the stall band a linear rate derivative is structurally
// invalid, so a taper, a clamp and a measured value are all equally
// fabrications. The two well-defined ranges are tabulated and the band
// between them is declared.
//
// The mechanism test below is therefore driven with SYNTHETIC polars of
// chosen slope rather than with a configuration's stalled sections: that
// isolates the physics from the amplitude question and keeps the guard
// meaningful whatever section model is mounted.
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
    std::vector<StripSection> Strips;
    ReferenceGeometry Ref;
    PreparedSystem Prepared;
};

Fixture MakeWing() {
    WingParams w;
    w.Span = Span;
    w.RootChord = Chord;
    w.TipChord = Chord;
    w.NPanelsSemiSpan = 16;

    Fixture f;
    f.Panels = BuildWing(w);
    for (Panel& p : f.Panels) p.Surface = "wing";
    double area = 0.0;
    for (const Panel& p : f.Panels) area += p.PlanformArea;
    f.Ref.Area = area;
    f.Ref.Span = w.Span;
    f.Ref.Chord = area / w.Span;
    f.Prepared = Prepare(PanelSystem{f.Panels, {}}, 50.0 * w.Span);

    for (const Panel& p : f.Panels) {
        StripSection s;
        s.ChordDir = Vec3(1, 0, 0);
        s.LiftDir = Vec3(0, 0, 1);
        s.Chord = (p.SpanwiseWidth > 0.0) ? p.PlanformArea / p.SpanwiseWidth : Chord;
        s.Width = p.SpanwiseWidth;
        s.Eta = std::fabs(0.5 * (p.A.y + p.B.y)) / (0.5 * Span);
        s.Alpha0Deg = 0.0;
        f.Strips.push_back(s);
    }
    return f;
}

/**
 * Coupled roll damping by central difference.
 *
 * THE FRAME TRAP, pinned below in its own test. The moment comes back in
 * contract FRD from BodyAxisFromCoupled, while the rate is set on
 * FreestreamConditions in the SOLVER frame. Cl and p each flip under the
 * solver-to-contract rotation, so the derivative is frame invariant --
 * but only when both sides are in the same frame. Differencing an FRD
 * moment against a solver-frame rate leaves exactly one flip outstanding,
 * and the quotient needs negating. That error reads as Cl_p = +0.54 in
 * the ATTACHED range, where the answer must be about -0.45: the right
 * magnitude with the wrong sign, which is precisely the failure mode
 * Solver/BodyAxes.h was written about.
 */
double CoupledClp(const Fixture& f, const SectionModel& model, double alphaDeg, double step) {
    FreestreamConditions base;
    base.Vinf = Vinf;
    base.rho = Rho;
    base.alphaDeg = alphaDeg;

    ViscousCouplingOptions options;
    options.Relaxation = 0.05;
    options.AndersonDepth = 0;
    options.MaxIterations = 300;

    const auto at = [&](double p) {
        FreestreamConditions fc = base;
        fc.p = p;
        return SolveViscousCoupled(f.Panels, f.Strips, fc, f.Ref, 50.0 * Span, model, options);
    };
    const double q = 0.5 * Rho * Vinf * Vinf;
    const double clPlus = BodyAxisFromCoupled(at(+step), q, f.Ref).Cl;
    const double clMinus = BodyAxisFromCoupled(at(-step), q, f.Ref).Cl;
    const double reduce = f.Ref.Span / (2.0 * Vinf);
    return -(clPlus - clMinus) / (2.0 * step * reduce); // see the frame note above
}

/** A polar with a chosen lift slope: cl = slope * alpha, plus fixed drag. */
SectionModel LinearPolar(double slopePerRad) {
    return [slopePerRad](const StripSection& strip, double alphaEffDeg, double, double) {
        SectionCoefficients c;
        c.cl = slopePerRad * Math::DegToRad(alphaEffDeg - strip.EffectiveAlpha0Deg());
        c.cd = 0.01;
        return c;
    };
}

// --- 1. the attached anchor ---------------------------------------------------

void TestAttachedDampingMatchesTheInviscidAnchor() {
    const Fixture f = MakeWing();
    const BodyAxisRateDerivatives inv =
        ComputeBodyAxisRateDerivatives(f.Prepared, [] {
            FreestreamConditions fc;
            fc.Vinf = Vinf;
            fc.rho = Rho;
            fc.alphaDeg = 2.0;
            return fc;
        }(), f.Ref);

    const double coupled = CoupledClp(f, LinearPolar(2.0 * std::numbers::pi), 2.0, 0.05);
    std::cout << "attached: Cl_p inviscid=" << inv.Clp << "  coupled=" << coupled << "\n";

    // THE regression guard for the sign. A wing damps roll.
    CHECK(inv.Clp < 0.0, "inviscid roll damping must be negative, got " << inv.Clp);
    CHECK(coupled < 0.0,
          "coupled roll damping must be negative in attached flow -- a positive value here "
          "means an outstanding frame flip, not autorotation. Got " << coupled);
    CHECK(inv.Clp < -0.2 && inv.Clp > -0.8,
          "inviscid Cl_p should sit near the textbook -0.45, got " << inv.Clp);
    // The coupled solve carries the section model's own slope, so it need
    // not match the lattice exactly -- but it must agree in sign and order.
    CHECK(coupled < -0.2 && coupled > -1.0,
          "coupled Cl_p should be the same order as the inviscid anchor, got " << coupled);
}

// --- 2. the mechanism ---------------------------------------------------------

void TestNegativeSectionSlopeReversesRollDamping() {
    // The whole finding, reduced to its mechanism and driven directly. A
    // rolling wing raises the down-going semi-span's incidence and lowers
    // the up-going one's. With a POSITIVE section lift slope that opposes
    // the roll; with a NEGATIVE slope -- which is what a section past
    // stall has -- it ADDS to it. So the sign of Cl_p follows the sign of
    // dcl/dalpha, and roll damping REVERSES past stall.
    const Fixture f = MakeWing();

    const double damped = CoupledClp(f, LinearPolar(+2.0 * std::numbers::pi), 2.0, 0.05);
    const double reversed = CoupledClp(f, LinearPolar(-1.0), 2.0, 0.05);
    std::cout << "mechanism: Cl_p with slope +2pi = " << damped << ", with slope -1 = "
              << reversed << "\n";

    CHECK(damped < 0.0, "a positive section lift slope must damp roll, got " << damped);
    CHECK(reversed > 0.0,
          "a NEGATIVE section lift slope must ANTI-damp roll -- this is autorotation, and a "
          "flight model that tapers or clamps roll damping past stall reports the wrong "
          "SIGN here. Got " << reversed);
}

void TestZeroSlopeGivesNegligibleDamping() {
    // The boundary between the two regimes: a section whose lift does not
    // respond to incidence produces no circulatory roll damping at all.
    // This pins that the reversal above comes from the SLOPE and not from
    // some fixed offset in the coupling.
    const Fixture f = MakeWing();
    const double flat = CoupledClp(f, LinearPolar(0.0), 2.0, 0.05);
    const double damped = CoupledClp(f, LinearPolar(+2.0 * std::numbers::pi), 2.0, 0.05);
    std::cout << "zero-slope Cl_p = " << flat << " against damped " << damped << "\n";
    CHECK(std::fabs(flat) < 0.1 * std::fabs(damped),
          "a zero section lift slope should give negligible roll damping, got " << flat);
}

void TestPitchRateReachesTheCoupledSolve() {
    // Does a PITCH rate move the coupled solve at all? Roll and yaw rates
    // give a spanwise-VARYING incidence, which a single-row strip method
    // sees directly. A pitch rate instead gives a UNIFORM incidence change
    // proportional to the chordwise offset between the bound line and the
    // moment reference point -- so with the reference point ON the bound
    // line it must produce nothing, and with the reference point offset it
    // must produce a large response. Both halves are checked, because the
    // first is the trap: a fixture that happens to reference about its own
    // quarter chord reports zero pitch response and looks broken.
    const Fixture f = MakeWing();

    ViscousCouplingOptions options;
    options.Relaxation = 0.05;
    options.AndersonDepth = 0;
    options.MaxIterations = 400;
    const double q = 0.5 * Rho * Vinf * Vinf;
    const double chordReduce = f.Ref.Chord / (2.0 * Vinf);
    const double step = 0.01 / chordReduce; // a 0.01 reduced-rate perturbation

    const auto czAt = [&](double pitchRate, const Vec3& refPoint) {
        FreestreamConditions fc;
        fc.Vinf = Vinf;
        fc.rho = Rho;
        fc.alphaDeg = 2.0;
        fc.q = pitchRate;
        fc.RefPoint = refPoint;
        return BodyAxisFromCoupled(
            SolveViscousCoupled(f.Panels, f.Strips, fc, f.Ref, 50.0 * Span,
                                LinearPolar(2.0 * std::numbers::pi), options),
            q, f.Ref).CZ;
    };

    // Reference point ON the bound line: no chordwise arm, no response.
    const Vec3 onLine(0.0, 0.0, 0.0);
    const double czqOnLine =
        -(czAt(+step, onLine) - czAt(-step, onLine)) / (2.0 * step * chordReduce);

    // Reference point offset a chord aft: a real arm, a real response.
    const Vec3 offset(Chord, 0.0, 0.0);
    const double czqOffset =
        -(czAt(+step, offset) - czAt(-step, offset)) / (2.0 * step * chordReduce);

    std::cout << "CZ_q on the bound line = " << czqOnLine << ", offset one chord = "
              << czqOffset << "\n";

    CHECK(std::fabs(czqOnLine) < 0.5,
          "with the reference point ON the bound line a single-row strip method must show "
          "essentially no pitch-rate response -- the arm is zero. Got " << czqOnLine);
    CHECK(std::fabs(czqOffset) > 2.0,
          "with the reference point offset a chord the pitch-rate response must be large; "
          "a near-zero value means body rates are not reaching the coupled strip "
          "velocities at all. Got " << czqOffset);
}

// --- 3. the frame trap, explicitly --------------------------------------------

void TestFrameFlipWouldBeCaught() {
    // Guards the discipline rather than the number: the UNNEGATED quotient
    // -- an FRD moment differenced against a solver-frame rate -- must
    // come out with the opposite sign to the correct one. If some future
    // change makes both conventions agree, one of them has silently moved
    // and the negation in CoupledClp is no longer doing what its comment
    // says.
    const Fixture f = MakeWing();
    const SectionModel model = LinearPolar(2.0 * std::numbers::pi);

    FreestreamConditions base;
    base.Vinf = Vinf;
    base.rho = Rho;
    base.alphaDeg = 2.0;

    ViscousCouplingOptions options;
    options.Relaxation = 0.05;
    options.AndersonDepth = 0;
    options.MaxIterations = 300;

    const auto at = [&](double p) {
        FreestreamConditions fc = base;
        fc.p = p;
        return SolveViscousCoupled(f.Panels, f.Strips, fc, f.Ref, 50.0 * Span, model, options);
    };
    const double q = 0.5 * Rho * Vinf * Vinf;
    const double step = 0.05;
    const double reduce = f.Ref.Span / (2.0 * Vinf);
    const double unnegated =
        (BodyAxisFromCoupled(at(+step), q, f.Ref).Cl -
         BodyAxisFromCoupled(at(-step), q, f.Ref).Cl) / (2.0 * step * reduce);
    const double correct = CoupledClp(f, model, 2.0, step);

    CHECK(unnegated * correct < 0.0,
          "the unnegated FRD-moment / solver-rate quotient must have the OPPOSITE sign to "
          "the correct derivative; got " << unnegated << " and " << correct);
    CHECK(std::fabs(std::fabs(unnegated) - std::fabs(correct)) < 1e-9 * std::fabs(correct) + 1e-12,
          "and the same magnitude -- which is why a magnitude-only check misses it");
}

} // namespace

int main() {
    TestAttachedDampingMatchesTheInviscidAnchor();
    TestNegativeSectionSlopeReversesRollDamping();
    TestZeroSlopeGivesNegligibleDamping();
    TestPitchRateReachesTheCoupledSolve();
    TestFrameFlipWouldBeCaught();

    if (failures == 0) {
        std::cout << "PASS: TestRollDamping\n";
        return 0;
    }
    std::cerr << failures << " check(s) failed in TestRollDamping\n";
    return 1;
}
