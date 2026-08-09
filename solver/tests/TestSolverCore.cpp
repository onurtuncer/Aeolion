// TestSolverCore.cpp -- validates Solver.h's core solver against known
// closed-form aerodynamic theory:
//   - CL for a rectangular unswept/untwisted wing should track the
//     lifting-line-corrected thin-wing formula CL = 2*pi*alpha/(1+2/AR)
//     within a few percent (VLM is more accurate than the 2-term formula,
//     so exact match isn't expected -- convergence in the right direction
//     as AR increases is the real check).
//   - Oswald efficiency for that same wing should sit close to 1.0.
//   - CL/CDi should converge monotonically as panel count increases.
#include "Aeolion/Solver/Solver.h"
#include <numbers>
#include <iostream>
#include <cmath>

using namespace Aeolion::Solver;

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } } while (0)

int main() {
    // --- thin-wing lift slope check across a few aspect ratios ---
    for (int AR_target : {6, 8, 12, 20}) {
        WingParams w;
        w.Span = 4.0;
        double c = w.Span / AR_target;
        w.RootChord = c; w.TipChord = c;
        w.SweepQuarterChordDeg = 0; w.DihedralDeg = 0; w.TwistTipDeg = 0;
        w.NPanelsSemiSpan = 30; w.CosineSpacing = true;

        FreestreamConditions fc; fc.Vinf = 30.0; fc.alphaDeg = 4.0; fc.rho = 1.225;
        SolveResult res = Solve(w, fc);
        double AR = w.Span * w.Span / res.ReferenceArea;
        double alphaRad = DegToRad(fc.alphaDeg);
        double CL_theory = Two * std::numbers::pi * alphaRad / (1.0 + Two / AR);
        double relErr = std::fabs(res.CL - CL_theory) / CL_theory;

        std::cout << "AR=" << AR << "  CL=" << res.CL << "  CL_theory=" << CL_theory
                  << "  relErr=" << relErr * 100 << "%\n";
        CHECK(relErr < 0.15, "CL vs thin-wing theory off by more than 15% at AR=" << AR);
        CHECK(res.CDi > 0, "CDi should be positive at AR=" << AR);

        double e = (res.CDi > 1e-9) ? (res.CL * res.CL) / (std::numbers::pi * AR * res.CDi) : 0.0;
        std::cout << "  Oswald e=" << e << "\n";
        CHECK(e > 0.85 && e < 1.05, "Oswald efficiency out of expected range for a plain rectangular wing, AR=" << AR);
    }

    // --- panel convergence: CL/CDi should settle down, not oscillate wildly ---
    double prevCL = 0;
    bool first = true;
    for (int npss : {10, 20, 40, 80}) {
        WingParams w;
        w.Span = 4.0; w.RootChord = 0.5; w.TipChord = 0.5;
        w.NPanelsSemiSpan = npss; w.CosineSpacing = true;
        FreestreamConditions fc; fc.Vinf = 30; fc.alphaDeg = 4.0; fc.rho = 1.225;
        SolveResult res = Solve(w, fc);
        std::cout << "N=" << 2 * npss << "  CL=" << res.CL << "\n";
        if (!first) {
            double delta = std::fabs(res.CL - prevCL);
            CHECK(delta < 0.05, "CL changed too much between panel refinements (N=" << 2 * npss << "), delta=" << delta);
        }
        prevCL = res.CL;
        first = false;
    }

    // --- stability derivatives: the nondimensional rate scaling ------------
    // The per-rad/s derivatives and their reduced-rate counterparts differ
    // by 2V/length, and getting that factor upside down is invisible in a
    // sign check while being wrong by (2V/b)^2 -- a factor of two thousand
    // at these numbers. So the reduced-rate derivatives are pinned against
    // textbook values, which is the only form that HAS textbook values:
    // roll damping of a plain rectangular wing sits near Cl_p = -0.45, and
    // an unswept wing at zero incidence has essentially no yaw-due-to-roll.
    {
        WingParams w;
        w.Span = 6.0; w.RootChord = 1.0; w.TipChord = 1.0;
        w.NPanelsSemiSpan = 20;
        FreestreamConditions fc; fc.Vinf = 25.0; fc.alphaDeg = 2.0; fc.rho = 1.225;

        std::vector<Panel> panels = BuildWing(w);
        for (auto& p : panels) p.Surface = "wing";
        double S = 0.0;
        for (const auto& p : panels) S += p.PlanformArea;
        ReferenceGeometry ref;
        ref.Area = S; ref.Span = w.Span; ref.Chord = S / w.Span;

        StabilityDerivatives d = ComputeDerivatives(panels, fc, ref, 50.0 * w.Span);
        std::cout << "AR=6 rectangular: CL_alpha=" << d.CL_alpha << "  Cl_p(nd)=" << d.Croll_p_nd
                  << "  Cm_q(nd)=" << d.Cm_q_nd << "\n";

        CHECK(d.CL_alpha > 4.0 && d.CL_alpha < 5.5,
              "AR=6 lift slope should sit near 4.5-5 per rad, got " << d.CL_alpha);
        CHECK(d.Croll_p_nd < -0.2 && d.Croll_p_nd > -0.8,
              "AR=6 roll damping should sit near Cl_p = -0.45, got " << d.Croll_p_nd
                  << " (a value near -1e-4 means the reduced-rate factor is inverted)");

        // The reduced-rate derivative is the per-rad/s one times 2V/length.
        // Stating it as an identity is what stops the factor being "fixed"
        // back to its reciprocal by inspection.
        const double expected = d.Croll_p * (2.0 * fc.Vinf / ref.Span);
        CHECK(std::fabs(d.Croll_p_nd - expected) < 1e-9 * std::fabs(expected) + 1e-12,
              "Croll_p_nd must equal Croll_p * 2V/b, got " << d.Croll_p_nd << " against "
                                                            << expected);

        // Pitch damping is EXACTLY zero here, and that is the correct
        // answer rather than a missing term. With one chordwise row every
        // bound vortex lies on the quarter-chord line and every control
        // point a half chord behind it, so a pitch rate about a reference
        // point on that line adds the SAME upwash at every control point:
        // a uniform change in effective incidence, which moves CL but
        // exerts no moment about the line it is applied on. Pitch damping
        // needs either a chordwise-resolved lattice or a reference point
        // off the quarter chord -- checked next.
        CHECK(std::fabs(d.Cm_q_nd) < 1e-9,
              "a single-row wing cannot damp pitch about its own quarter chord, got "
                  << d.Cm_q_nd);

        FreestreamConditions offset = fc;
        offset.RefPoint = Vec3(2.0, 0.0, 0.0); // two chords aft of the quarter-chord line
        StabilityDerivatives shifted = ComputeDerivatives(panels, offset, ref, 50.0 * w.Span);
        std::cout << "  about a point 2c aft: Cm_q(nd)=" << shifted.Cm_q_nd << "\n";
        CHECK(shifted.Cm_q_nd < 0.0,
              "with a real moment arm, pitch damping must be negative, got " << shifted.Cm_q_nd);
    }

    if (failures == 0) { std::cout << "PASS: TestSolverCore\n"; return 0; }
    std::cerr << failures << " check(s) failed in TestSolverCore\n";
    return 1;
}
