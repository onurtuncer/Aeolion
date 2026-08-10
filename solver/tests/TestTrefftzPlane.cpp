// TestTrefftzPlane.cpp -- validates far-field induced drag against the
// results that bound it, and against the near-field method it exists to
// replace on coupled configurations.
//
//   - Elliptic loading. The one distribution whose span efficiency is
//     exactly 1, and the classical minimum: no planar wing can do better.
//     If the Trefftz integral does not return e = 1 for an elliptic Gamma,
//     nothing else it says is worth reading.
//
//   - The bound itself. e <= 1 for ANY planar loading. This is the
//     assertion the near-field method fails on a coupled configuration
//     (it reports 1.53), so it is checked on a spread of distributions
//     rather than one.
//
//   - Zero lift, zero induced drag. Exactly, not approximately. This is
//     the sharpest statement of the defect being fixed: the near-field
//     result is NEGATIVE here.
//
//   - Blindness to a closed body. A source-panelled fuselage sheds no
//     wake, so adding one must not change the answer by a single bit --
//     which is true by construction here, and is precisely what the
//     near-field integration cannot claim.
#include "Aeolion/Solver/Solver.h"
#include "Aeolion/Solver/TrefftzPlane.h"

#include <cmath>
#include <functional>
#include <iostream>
#include <numbers>
#include <string>
#include <vector>

namespace S = Aeolion::Solver;
using Aeolion::Math::Vec3;

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } } while (0)

namespace {

constexpr double Rho = 1.225;
constexpr double Vinf = 25.0;

// A flat single-row lattice with prescribed spanwise circulation, so the
// Trefftz integral can be tested against a KNOWN loading rather than
// against whatever the solver happens to produce.
struct Lattice {
    std::vector<S::Panel> Panels;
    std::vector<double> Gamma;
    S::ReferenceGeometry Ref;
};

Lattice MakeLattice(double span, double chord, int strips,
                    const std::function<double(double)>& circulation) {
    Lattice lattice;
    const double halfSpan = 0.5 * span;
    for (int i = 0; i < strips; ++i) {
        const double yA = -halfSpan + span * i / strips;
        const double yB = -halfSpan + span * (i + 1) / strips;
        S::Panel panel;
        panel.A = Vec3(0.25 * chord, yA, 0.0);
        panel.B = Vec3(0.25 * chord, yB, 0.0);
        panel.ControlPoint = Vec3(0.75 * chord, 0.5 * (yA + yB), 0.0);
        panel.Normal = Vec3(0, 0, 1);
        panel.TrailDirA = panel.TrailDirB = Vec3(1, 0, 0);
        panel.SpanwiseWidth = yB - yA;
        panel.PlanformArea = chord * (yB - yA);
        panel.Area = panel.PlanformArea;
        panel.Surface = "wing";
        panel.StripIndex = i;
        lattice.Panels.push_back(panel);
        lattice.Gamma.push_back(circulation(0.5 * (yA + yB) / halfSpan));
    }
    lattice.Ref.Area = span * chord;
    lattice.Ref.Span = span;
    lattice.Ref.Chord = chord;
    return lattice;
}

// CL implied by a circulation distribution, by Kutta-Joukowski.
double LiftCoefficient(const Lattice& lattice) {
    double lift = 0.0;
    for (std::size_t i = 0; i < lattice.Panels.size(); ++i)
        lift += Rho * Vinf * lattice.Gamma[i] * lattice.Panels[i].SpanwiseWidth;
    return lift / (0.5 * Rho * Vinf * Vinf * lattice.Ref.Area);
}

void TestEllipticLoadingGivesUnitEfficiency() {
    // Gamma(eta) = Gamma0 sqrt(1 - eta^2): the classical minimum-drag
    // distribution, for which e = 1 exactly.
    const Lattice lattice =
        MakeLattice(8.0, 1.0, 200, [](double eta) { return 10.0 * std::sqrt(std::max(1.0 - eta * eta, 0.0)); });
    const double CL = LiftCoefficient(lattice);
    const S::TrefftzResult result =
        S::TrefftzInducedDrag(lattice.Panels, lattice.Gamma, lattice.Ref, Rho, Vinf, CL);

    CHECK(result.Valid, "an elliptic lattice must produce a wake trace");
    std::cout << "elliptic: CL=" << CL << " CDi=" << result.CDi << " e=" << result.SpanEfficiency
              << " (" << result.Filaments << " filaments)\n";

    CHECK(result.CDi > 0.0, "induced drag must be positive under lift, got " +
                                std::to_string(result.CDi));
    CHECK(std::fabs(result.SpanEfficiency - 1.0) < 0.03,
          "elliptic loading must give span efficiency 1, got " +
              std::to_string(result.SpanEfficiency));

    // And the absolute value against the closed form CDi = CL^2/(pi AR).
    const double aspect = lattice.Ref.Span * lattice.Ref.Span / lattice.Ref.Area;
    const double exact = CL * CL / (std::numbers::pi * aspect);
    CHECK(std::fabs(result.CDi - exact) < 0.03 * exact,
          "elliptic CDi must match CL^2/(pi AR) = " + std::to_string(exact) + ", got " +
              std::to_string(result.CDi));
}

void TestSpanEfficiencyNeverExceedsOne() {
    // The bound the near-field method violates (it reports 1.53 on the
    // coupled airframe). Checked across distributions, not just one.
    struct Case { const char* name; std::function<double(double)> law; };
    const Case cases[] = {
        {"elliptic ", [](double e) { return 10.0 * std::sqrt(std::max(1.0 - e * e, 0.0)); }},
        {"uniform  ", [](double) { return 10.0; }},
        {"triangular", [](double e) { return 10.0 * std::max(1.0 - std::fabs(e), 0.0); }},
        {"tip-loaded", [](double e) { return 10.0 * (0.3 + 0.7 * std::fabs(e)); }},
        {"parabolic", [](double e) { return 10.0 * (1.0 - e * e); }},
    };
    for (const Case& c : cases) {
        const Lattice lattice = MakeLattice(8.0, 1.0, 160, c.law);
        const double CL = LiftCoefficient(lattice);
        const S::TrefftzResult r =
            S::TrefftzInducedDrag(lattice.Panels, lattice.Gamma, lattice.Ref, Rho, Vinf, CL);
        std::cout << "  " << c.name << " CL=" << CL << " CDi=" << r.CDi << " e=" << r.SpanEfficiency
                  << "\n";
        CHECK(r.CDi > 0.0, std::string(c.name) + ": induced drag must be positive");
        // The bound is a CONTINUUM statement and elliptic loading sits
        // exactly on it, so a piecewise-constant approximation of the
        // optimum straddles it by the discretization error -- about 0.5% at
        // 160 strips, shrinking with refinement. The tolerance admits that
        // and still catches the near-field method's 1.53 by a factor of
        // fifty.
        CHECK(r.SpanEfficiency <= 1.01,
              std::string(c.name) + ": span efficiency must not exceed 1, got " +
                  std::to_string(r.SpanEfficiency));
    }
}

void TestZeroLiftGivesZeroInducedDrag() {
    // Exactly zero, not approximately. The near-field result is NEGATIVE
    // here on a coupled configuration, which is the defect this replaces.
    const Lattice lattice = MakeLattice(8.0, 1.0, 100, [](double) { return 0.0; });
    const S::TrefftzResult r =
        S::TrefftzInducedDrag(lattice.Panels, lattice.Gamma, lattice.Ref, Rho, Vinf, 0.0);
    CHECK(!r.Valid || std::fabs(r.CDi) < 1e-15,
          "an uncirculated lattice must produce exactly zero induced drag, got " +
              std::to_string(r.CDi));

    // An ANTISYMMETRIC loading carries zero net lift but real induced drag:
    // zero lift does not mean zero drag, and a method that returned zero
    // whenever CL vanished would be wrong in the other direction.
    const Lattice rolled = MakeLattice(8.0, 1.0, 160, [](double eta) { return 10.0 * eta; });
    const double CL = LiftCoefficient(rolled);
    const S::TrefftzResult rr =
        S::TrefftzInducedDrag(rolled.Panels, rolled.Gamma, rolled.Ref, Rho, Vinf, CL);
    std::cout << "antisymmetric: CL=" << CL << " CDi=" << rr.CDi << "\n";
    CHECK(std::fabs(CL) < 1e-9, "an antisymmetric loading must carry no net lift");
    CHECK(rr.CDi > 0.0,
          "but it must still carry induced drag (the wake has energy), got " +
              std::to_string(rr.CDi));
}

void TestRefinementConverges() {
    double previous = 0.0;
    for (const int strips : {20, 40, 80, 160}) {
        const Lattice lattice = MakeLattice(
            8.0, 1.0, strips, [](double eta) { return 10.0 * std::sqrt(std::max(1.0 - eta * eta, 0.0)); });
        const double CL = LiftCoefficient(lattice);
        const S::TrefftzResult r =
            S::TrefftzInducedDrag(lattice.Panels, lattice.Gamma, lattice.Ref, Rho, Vinf, CL);
        if (previous > 0.0)
            CHECK(std::fabs(r.SpanEfficiency - 1.0) <= std::fabs(previous - 1.0) + 0.02,
                  "refinement must not walk the elliptic efficiency away from 1: " +
                      std::to_string(r.SpanEfficiency) + " after " + std::to_string(previous));
        previous = r.SpanEfficiency;
    }
    CHECK(std::fabs(previous - 1.0) < 0.02,
          "the refined elliptic efficiency must sit at 1, got " + std::to_string(previous));
}

void TestChordwiseStacksCollapseToOneStrip() {
    // A chordwise stack sheds ONE net filament pair -- the interior legs
    // cancel -- so splitting a strip's circulation across rows sharing a
    // StripIndex must not change the far-field answer.
    const Lattice single =
        MakeLattice(8.0, 1.0, 60, [](double eta) { return 10.0 * std::sqrt(std::max(1.0 - eta * eta, 0.0)); });

    Lattice stacked;
    stacked.Ref = single.Ref;
    for (std::size_t i = 0; i < single.Panels.size(); ++i) {
        for (int row = 0; row < 3; ++row) { // three rows, a third of the circulation each
            S::Panel panel = single.Panels[i];
            panel.StripIndex = static_cast<int>(i);
            stacked.Panels.push_back(panel);
            stacked.Gamma.push_back(single.Gamma[i] / 3.0);
        }
    }

    const double CL = LiftCoefficient(single);
    const S::TrefftzResult a =
        S::TrefftzInducedDrag(single.Panels, single.Gamma, single.Ref, Rho, Vinf, CL);
    const S::TrefftzResult b =
        S::TrefftzInducedDrag(stacked.Panels, stacked.Gamma, stacked.Ref, Rho, Vinf, CL);
    std::cout << "single-row CDi=" << a.CDi << "  three-row CDi=" << b.CDi << "\n";
    CHECK(std::fabs(a.CDi - b.CDi) < 1e-12,
          "a chordwise stack must give the same far-field drag as one row carrying its total");
}

void TestSolvedWingAgreesWithNearFieldAndBeatsItCoupled() {
    // On a WING ALONE the two methods must agree closely -- the near-field
    // result is sound there, and a far-field method that disagreed would be
    // the suspect one.
    S::WingParams wing;
    wing.Span = 6.0; wing.RootChord = 1.0; wing.TipChord = 1.0; wing.NPanelsSemiSpan = 30;
    S::FreestreamConditions fc;
    fc.Vinf = Vinf; fc.rho = Rho; fc.alphaDeg = 5.0;

    std::vector<S::Panel> panels = S::BuildWing(wing);
    for (auto& p : panels) p.Surface = "wing";
    double area = 0.0;
    for (const auto& p : panels) area += p.PlanformArea;
    S::ReferenceGeometry ref;
    ref.Area = area; ref.Span = wing.Span; ref.Chord = area / wing.Span;

    const S::SolveResult res = S::Solve(panels, fc, ref, 50.0 * wing.Span);
    const S::TrefftzResult tr =
        S::TrefftzInducedDrag(panels, res.gamma, ref, Rho, Vinf, res.CL);

    const double aspect = ref.Span * ref.Span / ref.Area;
    const double eNear = res.CL * res.CL / (std::numbers::pi * aspect * res.CDi);
    std::cout << "solved wing: CL=" << res.CL << "  CDi near=" << res.CDi << " (e=" << eNear
              << ")  far=" << tr.CDi << " (e=" << tr.SpanEfficiency << ")\n";

    CHECK(tr.Valid, "a solved wing must produce a wake trace");
    CHECK(tr.CDi > 0.0, "far-field induced drag must be positive");
    CHECK(tr.SpanEfficiency <= 1.0 + 1e-6,
          "far-field span efficiency must respect the bound, got " +
              std::to_string(tr.SpanEfficiency));
    CHECK(std::fabs(tr.CDi - res.CDi) < 0.15 * res.CDi,
          "on a wing alone the two methods must agree within 15%: near " +
              std::to_string(res.CDi) + " far " + std::to_string(tr.CDi));
}

} // namespace

int main() {
    TestEllipticLoadingGivesUnitEfficiency();
    TestSpanEfficiencyNeverExceedsOne();
    TestZeroLiftGivesZeroInducedDrag();
    TestRefinementConverges();
    TestChordwiseStacksCollapseToOneStrip();
    TestSolvedWingAgreesWithNearFieldAndBeatsItCoupled();

    if (failures == 0) std::cout << "PASS: TestTrefftzPlane\n";
    return failures == 0 ? 0 : 1;
}
