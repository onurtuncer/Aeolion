// TestFlapSection.cpp -- the flap a strip carries in its SECTION
// description (Geometry::FlapEffectiveness, StripSection::
// EffectiveAlpha0Deg) and the control response it produces through the
// Level-2 coupling.
//
// This suite guards the fix for a measured blocker. A hinge cannot be
// represented geometrically on a single-row lattice -- PanelBuilder needs
// two chordwise rows to split a strip at the hinge, the coupling requires
// exactly one row per strip -- so a deflection through the coupled path
// produced no control effect at all. The resolution is that to a strip
// method a deflected flap IS a camber change, so it belongs in the
// zero-lift angle the section model is posed against.
//
// The checks below are the ones that would have caught the original
// defect, plus the two exact limits that make the thin-airfoil formula
// verifiable by inspection rather than by trusting a transcription.
#include "Aeolion/Geometry/FlapEffectiveness.h"
#include "Aeolion/Solver/BodyAxes.h"
#include "Aeolion/Solver/Solver.h"
#include "Aeolion/Solver/ViscousCoupling.h"

#include <cmath>
#include <iostream>
#include <numbers>
#include <vector>

using namespace Aeolion::Solver;
namespace Geometry = Aeolion::Geometry;

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

// --- 1. thin-airfoil effectiveness -------------------------------------------

void TestEffectivenessLimitsAreExact() {
    // Hinged at the leading edge, the "flap" is the whole section
    // rotating, so it is worth exactly one degree per degree.
    CHECK(std::fabs(Geometry::FlapEffectiveness(0.0) - 1.0) < 1e-12,
          "a flap hinged at the leading edge must have tau = 1, got "
              << Geometry::FlapEffectiveness(0.0));
    // A flap of zero chord does nothing.
    CHECK(std::fabs(Geometry::FlapEffectiveness(1.0)) < 1e-12,
          "a zero-chord flap must have tau = 0, got " << Geometry::FlapEffectiveness(1.0));
}

void TestEffectivenessAgainstTextbookValues() {
    // tau = 1 - theta/pi + sin(theta)/pi with theta = arccos(1 - 2 x_h).
    // The classical values: a half-chord flap is worth about 0.82, a
    // quarter-chord flap about 0.61.
    const double half = Geometry::FlapEffectiveness(0.5);
    const double quarter = Geometry::FlapEffectiveness(0.75);
    std::cout << "tau: half-chord flap " << half << ", quarter-chord " << quarter
              << ", 12% aileron " << Geometry::FlapEffectiveness(0.88) << "\n";

    CHECK(std::fabs(half - (1.0 - 0.5 + 1.0 / std::numbers::pi)) < 1e-12,
          "half-chord flap effectiveness must be 1 - 1/2 + 1/pi, got " << half);
    CHECK(std::fabs(quarter - 0.6090) < 5e-4,
          "quarter-chord flap effectiveness should be ~0.609, got " << quarter);
}

void TestEffectivenessIsStronglyConcave() {
    // THE number worth pinning for control sizing: a 12%-chord aileron is
    // worth 0.43, not 0.12. A linear-in-chord guess underestimates roll
    // authority by a factor of three and a half, and would look plausible.
    const double tau = Geometry::FlapEffectiveness(0.88);
    CHECK(tau > 0.40 && tau < 0.46,
          "the contract's 12% aileron should have tau ~ 0.43, got " << tau);
    CHECK(tau > 3.0 * 0.12,
          "flap effectiveness is strongly concave in chord fraction: tau must far exceed "
          "the chord fraction itself, got " << tau);

    // Monotone decreasing as the hinge moves aft.
    double previous = 2.0;
    for (double xh = 0.0; xh <= 1.0 + 1e-12; xh += 0.05) {
        const double t = Geometry::FlapEffectiveness(xh);
        CHECK(t <= previous + 1e-12, "tau must decrease as the hinge moves aft, at x_h = " << xh);
        previous = t;
    }
}

void TestZeroLiftShiftSign() {
    // Trailing edge DOWN adds camber, which LOWERS the zero-lift angle.
    const double shift = Geometry::FlapZeroLiftShift(0.88, 10.0);
    CHECK(shift < 0.0, "a positive (TE-down) deflection must lower alpha0, got " << shift);
    CHECK(std::fabs(shift + 10.0 * Geometry::FlapEffectiveness(0.88)) < 1e-12,
          "the shift must be exactly -tau * deflection");
}

// --- 2. the strip carries it --------------------------------------------------

void TestEffectiveAlpha0IsAdditiveAndInert() {
    StripSection strip;
    strip.Alpha0Deg = -4.2;

    // The property that makes this change safe: a strip with no flap is
    // bit-identical to what every existing consumer already had.
    CHECK(strip.EffectiveAlpha0Deg() == strip.Alpha0Deg,
          "with no flap, EffectiveAlpha0Deg must be exactly Alpha0Deg");

    strip.FlapChordFraction = 0.12;
    strip.FlapDeflectionDeg = 0.0;
    CHECK(strip.EffectiveAlpha0Deg() == strip.Alpha0Deg,
          "a flap at zero deflection must change nothing");

    strip.FlapDeflectionDeg = 10.0;
    const double expected = -4.2 - 10.0 * Geometry::FlapEffectiveness(0.88);
    CHECK(std::fabs(strip.EffectiveAlpha0Deg() - expected) < 1e-12,
          "flapped zero-lift angle must be alpha0 - tau*delta, got "
              << strip.EffectiveAlpha0Deg() << " against " << expected);

    // And it is antisymmetric in the deflection, which is what makes an
    // aileron an aileron.
    strip.FlapDeflectionDeg = -10.0;
    CHECK(std::fabs(strip.EffectiveAlpha0Deg() - (-4.2 + 10.0 * Geometry::FlapEffectiveness(0.88)))
              < 1e-12,
          "the shift must be antisymmetric in deflection");
}

// --- 3. the coupled control response ------------------------------------------

struct Fixture {
    std::vector<Panel> Panels;
    std::vector<StripSection> Strips;
    ReferenceGeometry Ref;
};

Fixture MakeWing(double aileronDeg) {
    WingParams w;
    w.Span = Span;
    w.RootChord = Chord;
    w.TipChord = Chord;
    w.NPanelsSemiSpan = 20;

    Fixture f;
    f.Panels = BuildWing(w);
    for (Panel& p : f.Panels) p.Surface = "wing";
    double area = 0.0;
    for (const Panel& p : f.Panels) area += p.PlanformArea;
    f.Ref.Area = area;
    f.Ref.Span = w.Span;
    f.Ref.Chord = area / w.Span;

    for (const Panel& p : f.Panels) {
        const double y = 0.5 * (p.A.y + p.B.y);
        StripSection s;
        s.ChordDir = Vec3(1, 0, 0);
        s.LiftDir = Vec3(0, 0, 1);
        s.Chord = (p.SpanwiseWidth > 0.0) ? p.PlanformArea / p.SpanwiseWidth : Chord;
        s.Width = p.SpanwiseWidth;
        s.Eta = std::fabs(y) / (0.5 * Span);
        s.Alpha0Deg = 0.0;
        if (s.Eta >= 0.7 && aileronDeg != 0.0) { // outboard 30%, antisymmetric
            s.FlapChordFraction = 0.25;
            s.FlapDeflectionDeg = (y >= 0.0) ? aileronDeg : -aileronDeg;
        }
        f.Strips.push_back(s);
    }
    return f;
}

double SolveRollingMoment(const Fixture& f, double alphaDeg) {
    FreestreamConditions fc;
    fc.Vinf = Vinf;
    fc.rho = Rho;
    fc.alphaDeg = alphaDeg;

    ViscousCouplingOptions options;
    options.Relaxation = 0.05;
    options.AndersonDepth = 0;
    options.MaxIterations = 1000;

    const ViscousCoupledResult res = SolveViscousCoupled(
        f.Panels, f.Strips, fc, f.Ref, 50.0 * Span, AnalyticSectionModel{}, options);
    const double q = 0.5 * Rho * Vinf * Vinf;
    // Via the tested extractor: the coupled result's own coefficient
    // members are never populated, and reading them is what made the
    // original measurement report "no control effect at all".
    return BodyAxisFromCoupled(res, q, f.Ref).Cl;
}

void TestDeflectionProducesRollingMoment() {
    // THE regression guard. Before the fix this difference was exactly
    // zero at every attitude, silently, with the solve converging and
    // reporting sensible forces.
    const Fixture neutral = MakeWing(0.0);
    const Fixture rolled = MakeWing(10.0);

    const double clNeutral = SolveRollingMoment(neutral, 4.0);
    const double clRolled = SolveRollingMoment(rolled, 4.0);
    const double delta = clRolled - clNeutral;
    std::cout << "coupled aileron at 4 deg: Cl neutral=" << clNeutral << " rolled=" << clRolled
              << " delta=" << delta << "\n";

    CHECK(std::fabs(clNeutral) < 1e-9,
          "a symmetric wing with no deflection must carry no rolling moment, got " << clNeutral);
    CHECK(std::fabs(delta) > 1e-4,
          "an antisymmetric deflection MUST produce a rolling moment through the coupled "
          "solve -- exactly zero is the defect this suite exists for, got " << delta);

    // Opposite commands must roll opposite ways, by the same amount.
    const Fixture antiRolled = MakeWing(-10.0);
    const double clAnti = SolveRollingMoment(antiRolled, 4.0) - clNeutral;
    CHECK(delta * clAnti < 0.0, "opposite aileron commands must roll opposite ways");
    CHECK(std::fabs(std::fabs(clAnti) - std::fabs(delta)) < 0.05 * std::fabs(delta),
          "and by the same magnitude on a symmetric wing, got " << clAnti << " against "
                                                                << delta);
}

void TestAuthorityDecaysThroughStall() {
    // The question the whole exercise started from: roll authority must
    // fall as the wing stalls. The inviscid lattice cannot show this --
    // it contains no stall -- so a decay here is evidence the section
    // model is genuinely in charge of the control effect.
    const Fixture neutral = MakeWing(0.0);
    const Fixture rolled = MakeWing(10.0);

    const double attached =
        std::fabs(SolveRollingMoment(rolled, 0.0) - SolveRollingMoment(neutral, 0.0));
    const double stalled =
        std::fabs(SolveRollingMoment(rolled, 30.0) - SolveRollingMoment(neutral, 30.0));
    std::cout << "roll authority: |dCl| at alpha 0 = " << attached << ", at 30 = " << stalled
              << " (retained " << 100.0 * stalled / attached << "%)\n";

    CHECK(stalled < attached,
          "roll authority must decay through stall, got " << stalled << " against " << attached);
    CHECK(stalled < 0.6 * attached,
          "and the decay must be substantial, not marginal: a table retaining most of its "
          "attached-flow roll power at 30 degrees models a simulator that recovers from "
          "departures it should not. Retained "
              << 100.0 * stalled / attached << "%");
}

} // namespace

int main() {
    TestEffectivenessLimitsAreExact();
    TestEffectivenessAgainstTextbookValues();
    TestEffectivenessIsStronglyConcave();
    TestZeroLiftShiftSign();
    TestEffectiveAlpha0IsAdditiveAndInert();
    TestDeflectionProducesRollingMoment();
    TestAuthorityDecaysThroughStall();

    if (failures == 0) {
        std::cout << "PASS: TestFlapSection\n";
        return 0;
    }
    std::cerr << failures << " check(s) failed in TestFlapSection\n";
    return 1;
}
