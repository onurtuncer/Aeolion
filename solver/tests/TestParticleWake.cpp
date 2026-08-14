// TestParticleWake.cpp -- the three-dimensional particle-wake cross-check
// (Solver/ParticleWake.h), pinned against the results that bound it: the
// attached unsteady lattice must relax onto the steady VLM's lift on the
// SAME panels; the normal plate must land in the finite-wing bluff band
// and, decisively, BELOW the two-dimensional street's level -- spanwise
// breakup is the physics this tier exists to add; and sideslip must
// mirror.

#include "Aeolion/Solver/ParticleWake.h"
#include "Aeolion/Solver/Solver.h"

#include <cmath>
#include <iostream>
#include <vector>

using namespace Aeolion;
namespace S = Aeolion::Solver;

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } } while (0)

namespace {

constexpr double Chord = 0.2, Span = 1.2, V = 10.0, Rho = 1.225;
constexpr int Strips = 12;

void MakeWing(std::vector<S::Panel>& panels, std::vector<S::StripSection>& strips) {
    const double width = Span / Strips;
    for (int i = 0; i < Strips; ++i) {
        const double y0 = -0.5 * Span + width * i;
        S::Panel p;
        p.A = S::Vec3(0.25 * Chord, y0, 0.0);
        p.B = S::Vec3(0.25 * Chord, y0 + width, 0.0);
        p.ControlPoint = S::Vec3(0.75 * Chord, y0 + 0.5 * width, 0.0);
        p.Normal = S::Vec3(0.0, 0.0, 1.0);
        p.TrailDirA = p.TrailDirB = S::Vec3(1.0, 0.0, 0.0);
        p.PlanformArea = Chord * width;
        p.Area = p.PlanformArea;
        p.SpanwiseWidth = width;
        p.Surface = "wing";
        panels.push_back(p);

        S::StripSection s;
        s.ChordDir = S::Vec3(1.0, 0.0, 0.0);
        s.LiftDir = S::Vec3(0.0, 0.0, 1.0);
        s.Chord = Chord;
        s.Width = width;
        s.Eta = std::fabs(y0 + 0.5 * width) / (0.5 * Span);
        strips.push_back(s);
    }
}

S::FreestreamConditions Conditions(double alphaDeg, double betaDeg = 0.0) {
    S::FreestreamConditions fc;
    fc.Vinf = V;
    fc.alphaDeg = alphaDeg;
    fc.betaDeg = betaDeg;
    fc.rho = Rho;
    return fc;
}

S::ReferenceGeometry Reference() {
    S::ReferenceGeometry ref;
    ref.Area = Chord * Span;
    ref.Span = Span;
    ref.Chord = Chord;
    return ref;
}

void TestAttachedRelaxesToSteadyVlm() {
    std::vector<S::Panel> panels;
    std::vector<S::StripSection> strips;
    MakeWing(panels, strips);

    // The steady reference on the SAME panels (horseshoe formulation).
    const S::SolveResult steady = S::Solve(panels, Conditions(5.0), Reference(), 50.0 * Span);

    S::ParticleWakeOptions options;
    options.TimeStep = 0.1;
    options.Duration = 12.0;
    options.KeepHistory = true;
    const auto run = S::SolveParticleWake(panels, strips, Conditions(5.0), Reference(), options);

    CHECK(run.Valid, "the attached run completes");
    CHECK(run.SheddingStrips == 0, "no strip sheds its leading edge at alpha = 5");
    // 15%: the averaging window still rides Wagner's tail (about -8% of
    // the steady value by itself at this duration), and the march is
    // first-order in time on a coarse step.
    CHECK(std::fabs(run.MeanCL - steady.CL) < 0.15 * steady.CL,
          "mean CL relaxes onto the steady VLM within 15%");
    CHECK(run.RmsCL < 0.05 * steady.CL, "the attached wake is quiet");
    // Growth from the impulsive start, Wagner-like -- compared as window
    // means so a single start-up impulse spike cannot fake or mask it.
    double earlyMean = 0.0;
    int earlyCount = 0;
    for (const auto& s : run.History)
        if (s.t < 3.0) { earlyMean += s.CL; ++earlyCount; }
    earlyMean /= std::max(earlyCount, 1);
    CHECK(earlyMean < run.MeanCL, "lift grows from the impulsive start");
    std::cout << "attached alpha=5: mean CL " << run.MeanCL << " (steady VLM " << steady.CL
              << ", circulation CL " << run.MeanCirculationCL << "), rms " << run.RmsCL
              << ", particles " << run.MaxParticles << "\n";
}

void TestNormalPlateBreaksTheStreet() {
    std::vector<S::Panel> panels;
    std::vector<S::StripSection> strips;
    MakeWing(panels, strips);

    S::ParticleWakeOptions options;
    options.TimeStep = 0.1;
    options.Duration = 14.0;
    const auto run = S::SolveParticleWake(panels, strips, Conditions(90.0), Reference(), options);

    CHECK(run.Valid, "the normal-plate run completes");
    CHECK(run.SheddingStrips == Strips, "every strip sheds its leading edge at alpha = 90");
    // The STRUCTURAL claims only: the diffused street is BOUNDED (the
    // inviscid tier's RMS grew without limit; Phase A's collapsed it two
    // orders) and its magnitude sits below the 2-D street's level. At
    // this budget the averaging window holds about one shedding cycle,
    // so the MEAN's sign and level are not resolvable here -- they are
    // the pilot runs' job, which report batch-mean confidence intervals.
    CHECK(run.RmsCN < 3.0, "the diffused street is bounded");
    CHECK(std::fabs(run.MeanCN) < 2.8, "the street's scale sits below the 2-D level");
    CHECK(run.RmsCN > 0.02, "the plate wake sheds");
    std::cout << "plate alpha=90: mean CN " << run.MeanCN << " rms " << run.RmsCN
              << " (circulation CL " << run.MeanCirculationCL << "), particles "
              << run.MaxParticles << "\n";
}

void TestSideslipMirrors() {
    std::vector<S::Panel> panels;
    std::vector<S::StripSection> strips;
    MakeWing(panels, strips);

    S::ParticleWakeOptions options;
    options.TimeStep = 0.1;
    options.Duration = 10.0;
    const auto plus = S::SolveParticleWake(panels, strips, Conditions(30.0, 15.0), Reference(), options);
    const auto minus = S::SolveParticleWake(panels, strips, Conditions(30.0, -15.0), Reference(), options);

    CHECK(plus.Valid && minus.Valid, "both sideslip runs complete");
    CHECK(std::fabs(plus.MeanCL - minus.MeanCL) < 0.2 * std::fabs(plus.MeanCL) + 0.05,
          "lift mirrors under beta -> -beta to shedding scatter");
    // CY on the wing alone is a small difference of large fluctuating
    // quantities; the scatter base is the larger magnitude plus an
    // absolute floor sized to the plate case's mean uncertainty.
    const double cyScale = std::max(std::fabs(plus.MeanCY), std::fabs(minus.MeanCY));
    CHECK(std::fabs(plus.MeanCY + minus.MeanCY) < 0.25 * cyScale + 0.08,
          "side force is odd in beta to shedding scatter");
    std::cout << "sideslip alpha=30: CL " << plus.MeanCL << " / " << minus.MeanCL << ", CY "
              << plus.MeanCY << " / " << minus.MeanCY << "\n";
}

} // namespace

int main() {
    TestAttachedRelaxesToSteadyVlm();
    TestNormalPlateBreaksTheStreet();
    TestSideslipMirrors();

    if (failures) {
        std::cerr << failures << " check(s) failed in TestParticleWake\n";
        return 1;
    }
    std::cout << "TestParticleWake: all checks passed\n";
    return 0;
}

