// TestPostStallSection.cpp -- the anchored post-stall section model
// (Solver/PostStallSection.h) and the section-cm plumbing through the
// Level-2 coupling.
//
// The checks pin the model to its ANCHORS -- the closed-form classical
// results the whole Phase-1 construction exists to be anchored on -- not to
// its own outputs: Viterna's CdMax fit and 90-degree values, Hoerner's 1.98
// plate normal force, Kirchhoff's exact attached and fully separated
// limits, Rayleigh's centre of pressure, continuity at the emergent stall
// junction, and the quarter-chord couple arriving in the coupled solve's
// moments exactly as q c^2 w cm.

#include "Aeolion/Solver/PostStallSection.h"
#include "Aeolion/Solver/Solver.h"
#include "Aeolion/Solver/ViscousCoupling.h"

#include <cmath>
#include <iostream>
#include <numbers>
#include <vector>

using namespace Aeolion;
namespace S = Aeolion::Solver;

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } } while (0)

namespace {

constexpr double Pi = std::numbers::pi;

bool Near(double a, double b, double tol) { return std::fabs(a - b) <= tol; }

// A separation schedule shaped like the computed surveys: attached to
// 6 degrees, then the separation point walks forward and is fully
// separated by 35.
double RampSeparation(double /*eta*/, double alphaFromZeroLiftDeg) {
    const double a = std::fabs(alphaFromZeroLiftDeg);
    if (a <= 6.0) return 1.0;
    if (a >= 35.0) return 0.0;
    return 1.0 - (a - 6.0) / 29.0;
}

void TestClassicalAnchors() {
    // Viterna's fit, its 2-D consistency with Hoerner, and the 90-degree values.
    CHECK(Near(S::ViternaCdMax(6.0), 1.218, 1e-12), "CdMax(AR=6) = 1.11 + 0.018*6 = 1.218");
    CHECK(Near(S::ViternaCdMax(50.0), 2.01, 1e-12), "CdMax caps at the AR=50 edge of the fit");
    CHECK(Near(S::ViternaCdMax(500.0), 2.01, 1e-12), "beyond AR=50 the fit is held, not extrapolated");
    CHECK(Near(S::HoernerPostStallCn(Pi / 2.0), 1.0 / 0.505, 1e-9),
          "Hoerner cn(90) = 1/(0.222+0.283) = 1.980");
    CHECK(Near(S::ViternaCdMax(50.0), S::HoernerPostStallCn(Pi / 2.0), 0.05),
          "the AR->50 Viterna ceiling meets Hoerner's 2-D plate within 0.05");

    // Kirchhoff attenuation limits, and the classical fully separated
    // lift slope pi*alpha/2 (a quarter of the attached 2*pi*alpha).
    CHECK(Near(S::KirchhoffAttenuation(1.0), 1.0, 1e-15), "K(1) = 1: attached");
    CHECK(Near(S::KirchhoffAttenuation(0.0), 0.25, 1e-15), "K(0) = 1/4: the pi*alpha/2 plate slope");

    // Rayleigh centre of pressure: 5/16 at zero incidence, mid-chord at 90,
    // monotone between.
    CHECK(Near(S::RayleighCenterOfPressure(0.0), 5.0 / 16.0, 1e-12), "Rayleigh xcp(0) = 5/16");
    CHECK(Near(S::RayleighCenterOfPressure(Pi / 2.0), 0.5, 1e-12), "Rayleigh xcp(90) = 1/2");
    double prev = 0.0;
    bool monotone = true;
    for (double aDeg = 0.0; aDeg <= 90.0; aDeg += 1.0) {
        const double x = S::RayleighCenterOfPressure(aDeg * Pi / 180.0);
        if (x < prev - 1e-12) monotone = false;
        prev = x;
    }
    CHECK(monotone, "Rayleigh xcp walks monotonically aft");

    // Viterna endpoint behaviour under an arbitrary junction.
    const S::ViternaConstants v = S::MatchViterna(6.0, 20.0 * Pi / 180.0, 1.2, 0.3);
    CHECK(Near(S::ViternaCl(v, 20.0 * Pi / 180.0), 1.2, 1e-9), "Viterna cl continuous at the junction");
    CHECK(Near(S::ViternaCd(v, 20.0 * Pi / 180.0), 0.3, 1e-9), "Viterna cd continuous at the junction");
    CHECK(std::fabs(S::ViternaCl(v, Pi / 2.0)) < 1e-9, "Viterna cl(90) = 0");
    CHECK(Near(S::ViternaCd(v, Pi / 2.0), v.CdMax, 1e-9), "Viterna cd(90) = CdMax");
}

void TestAttachedLimitIsExactThinAirfoil() {
    // With no separation function the branch must reproduce cl = 2 pi sin a
    // and carry NO pressure drag -- only the viscous polar.
    S::PostStallSectionModel model;
    model.SeparationPoint = nullptr;
    S::StripSection strip;
    strip.Eta = 0.5;
    strip.Alpha0Deg = 0.0;

    for (double aDeg = 0.0; aDeg <= 10.0; aDeg += 2.0) {
        const auto c = model(strip, aDeg, model.ReferenceReynolds, 0.0);
        const double a = aDeg * Pi / 180.0;
        CHECK(Near(c.cl, 2.0 * Pi * std::sin(a), 1e-9), "attached cl = 2 pi sin(alpha)");
        const double cdPolar = model.Cd0 + model.KCd * c.cl * c.cl;
        CHECK(Near(c.cd, cdPolar, 1e-9), "attached cd is the viscous polar alone (d'Alembert)");
        CHECK(Near(c.cm, 0.0, 1e-12), "attached cm about c/4 is zero");
    }

    // Odd symmetry.
    const auto plus = model(strip, 8.0, model.ReferenceReynolds, 0.0);
    const auto minus = model(strip, -8.0, model.ReferenceReynolds, 0.0);
    CHECK(Near(plus.cl, -minus.cl, 1e-12), "cl is odd in alpha");
    CHECK(Near(plus.cd, minus.cd, 1e-12), "cd is even in alpha");
    CHECK(Near(plus.cm, -minus.cm, 1e-12), "cm is odd in alpha");

    // The camber convention: alpha is measured from the zero-lift line.
    S::StripSection cambered = strip;
    cambered.Alpha0Deg = -4.0;
    const auto atZeroLift = model(cambered, -4.0, model.ReferenceReynolds, 0.0);
    CHECK(std::fabs(atZeroLift.cl) < 1e-12, "cl vanishes at the zero-lift angle");
}

void TestStallEmergesFromSeparation() {
    S::PostStallSectionModel model;
    model.SeparationPoint = RampSeparation;
    S::StripSection strip;
    strip.Eta = 0.5;

    // The stall angle must be an interior peak of the Kirchhoff branch --
    // an OUTPUT of the separation schedule, not a constant.
    const double aStall = model.StallAngleDeg(strip.Eta);
    CHECK(aStall > 6.0 && aStall < 35.0, "stall angle emerges inside the separation ramp");

    // cl rises to the junction, then decays toward deep stall; cd grows
    // monotonically past the junction toward CdMax at 90.
    double clMax = 0.0, aOfMax = 0.0;
    double prevCd = 0.0;
    bool cdMonotonePastStall = true;
    for (double aDeg = 1.0; aDeg <= 90.0; aDeg += 1.0) {
        const auto c = model(strip, aDeg, model.ReferenceReynolds, 0.0);
        if (c.cl > clMax) { clMax = c.cl; aOfMax = aDeg; }
        if (aDeg > aStall + 1.0 && c.cd < prevCd - 1e-9) cdMonotonePastStall = false;
        prevCd = c.cd;
    }
    CHECK(Near(aOfMax, aStall, 1.0), "the lift peak sits at the reported stall angle");
    CHECK(clMax < 2.0 * Pi * std::sin(aStall * Pi / 180.0),
          "CLmax is attenuated below the attached value");

    const auto at90 = model(strip, 90.0, model.ReferenceReynolds, 0.0);
    CHECK(Near(at90.cd, S::ViternaCdMax(model.AspectRatio), 0.02),
          "cd(90) lands on the AR-aware Viterna ceiling");
    CHECK(std::fabs(at90.cl) < 0.02, "cl(90) ~ 0");
    CHECK(cdMonotonePastStall, "cd grows monotonically through deep stall");

    // Continuity at the junction: the branch hand-off may not jump.
    const auto below = model(strip, aStall - 0.05, model.ReferenceReynolds, 0.0);
    const auto above = model(strip, aStall + 0.05, model.ReferenceReynolds, 0.0);
    CHECK(Near(below.cl, above.cl, 0.02), "cl continuous across the stall junction");
    CHECK(Near(below.cd, above.cd, 0.02), "cd continuous across the stall junction");

    // The centre of pressure walks aft: nose-down cm beyond stall, deeper
    // with incidence, reaching the Rayleigh mid-chord couple at 90.
    const auto deep = model(strip, 60.0, model.ReferenceReynolds, 0.0);
    CHECK(deep.cm < -0.01, "separated flow carries a nose-down quarter-chord cm");
    const double cn90 = at90.cl * 0.0 + (at90.cd - model.Cd0) * 1.0; // cn = cd at 90
    CHECK(Near(at90.cm, -(0.5 - 0.25) * cn90, 0.05), "cm(90) is the mid-chord couple of cn = cd");
}

void TestSectionMomentPlumbing() {
    // Four flat strips; a section model with zero lift and drag but a fixed
    // cm must produce EXACTLY the quarter-chord couple sum q c^2 w cm about
    // +y, and nothing else.
    constexpr double Chord = 0.4, Width = 0.5, V = 20.0, Rho = 1.2, CmFixed = 0.08;
    std::vector<S::Panel> panels;
    std::vector<S::StripSection> strips;
    for (int i = 0; i < 4; ++i) {
        S::Panel p;
        const double y0 = -1.0 + Width * i;
        p.A = S::Vec3(0.0, y0, 0.0);
        p.B = S::Vec3(0.0, y0 + Width, 0.0);
        p.ControlPoint = S::Vec3(0.5 * Chord, y0 + 0.5 * Width, 0.0);
        p.Normal = S::Vec3(0.0, 0.0, 1.0);
        p.TrailDirA = p.TrailDirB = S::Vec3(1.0, 0.0, 0.0);
        p.PlanformArea = Chord * Width;
        p.Area = p.PlanformArea;
        p.SpanwiseWidth = Width;
        p.Surface = "wing";
        panels.push_back(p);

        S::StripSection s;
        s.ChordDir = S::Vec3(1.0, 0.0, 0.0);
        s.LiftDir = S::Vec3(0.0, 0.0, 1.0);
        s.Chord = Chord;
        s.Width = Width;
        s.Eta = 0.5;
        strips.push_back(s);
    }

    const S::SectionModel model = [](const S::StripSection&, double, double, double) {
        return S::SectionCoefficients{0.0, 0.0, CmFixed};
    };

    S::FreestreamConditions fc;
    fc.Vinf = V;
    fc.alphaDeg = 0.0;
    fc.rho = Rho;
    S::ReferenceGeometry ref;
    ref.Area = 4.0 * Chord * Width;
    ref.Span = 2.0;
    ref.Chord = Chord;

    const auto res = S::SolveViscousCoupled(panels, strips, fc, ref, 100.0, model);
    CHECK(res.Converged, "the zero-lift fixed point converges trivially");

    const double q = 0.5 * Rho * V * V;
    const double expected = 4.0 * q * Chord * Chord * Width * CmFixed;
    CHECK(Near(res.SectionMoment.y, expected, 0.01 * expected),
          "SectionMoment.y = sum q c^2 w cm (positive cm pitches nose-up about +y)");
    CHECK(std::fabs(res.SectionMoment.x) < 1e-9 && std::fabs(res.SectionMoment.z) < 1e-9,
          "the couple has no roll/yaw component on a flat wing");
    CHECK(Near(res.Base.My, expected, 0.02 * expected),
          "the couple arrives in the total pitching moment");
    std::cout << "section couple: " << res.SectionMoment.y << " N*m (expected " << expected
              << ")\n";
}

} // namespace

int main() {
    TestClassicalAnchors();
    TestAttachedLimitIsExactThinAirfoil();
    TestStallEmergesFromSeparation();
    TestSectionMomentPlumbing();

    if (failures) {
        std::cerr << failures << " check(s) failed in TestPostStallSection\n";
        return 1;
    }
    std::cout << "TestPostStallSection: all checks passed\n";
    return 0;
}
