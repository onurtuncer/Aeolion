// Solver/TrefftzPlane.h
//
// Induced drag from the FAR field: the kinetic energy the wake carries off,
// measured in a plane far downstream, instead of the streamwise component of
// the near-field forces.
//
// --- why the near-field answer is not enough ---------------------------------
// Solver.h integrates forces by applying Kutta-Joukowski at each bound-vortex
// midpoint, which gives lift and induced drag together and needs no separate
// wake integration. On a WING that is fine. On a coupled configuration --
// lifting surfaces plus a closed source-panelled body -- it is not, and the
// failure is not subtle:
//
//   - CDi comes out NEGATIVE at zero lift, where induced drag must vanish;
//   - fitting CDi = CDi0 + k CL^2 over an alpha sweep of the airframe in
//     tests/Data implies an Oswald efficiency of 1.53, and e <= 1 for any
//     planar wing.
//
// The cause is NOT a violation of d'Alembert, and it is worth saying so
// because that is the tempting explanation and it is measurably false:
// closed bodies in this solver carry zero net force in uniform flow to
// MACHINE PRECISION (|F|/qA ~ 1e-16, TestBodyPanels and TestDuctPanels).
// The body's force in a coupled solve is physical -- it sits in the wing's
// upwash and carries about 8% of the lift on the airframe in tests/Data.
//
// The cause is that near-field induced drag is a small difference of much
// larger quantities. It is the streamwise component of forces dominated by
// lift, so its relative error is amplified by the ratio of lift to drag,
// and it converges slowly for that reason alone. On a wing this is benign
// -- near and far field agree to 0.4% here. Adding a body changes the
// induced velocity at the wing's bound vortices substantially, most of all
// near the root, and adds a pressure integration over the body in a
// strongly non-uniform field. Both errors are negligible beside lift, which
// is why CL and the moments are unaffected, and comparable with induced
// drag, which is a far smaller number.
//
// --- what the far field does about it -----------------------------------------
// The Trefftz plane measures something the body cannot contribute to. A
// closed non-lifting body sheds no wake: its source distribution has no
// trailing vorticity, so it puts nothing through a plane at downstream
// infinity. Only the lifting surfaces' trailing vorticity crosses it. The
// far-field integral therefore never evaluates the body at all, and cannot
// inherit the near-field evaluation's error, by construction rather than by
// a correction that has to be tuned.
//
// The classical result (Munk; Katz & Plotkin) is that the induced drag equals
// the crossflow kinetic energy left in the wake, which for a wake trace of
// potential jump Gamma(y) and downwash w_T(y) reduces to
//
//     D_i = -(rho/2) integral Gamma(y) w_T(y) dy,
//
// with w_T the downwash induced IN THE TREFFTZ PLANE -- twice the downwash at
// the lifting line, since the wake there is doubly infinite rather than
// semi-infinite. That factor of two is the whole reason this is not simply
// the near-field calculation rearranged, and getting it wrong halves the
// answer while leaving every trend intact.
//
// --- the discrete form ---------------------------------------------------------
// A chordwise stack of panels sharing a StripIndex sheds ONE net filament
// pair into the far wake: the interior legs of the stack cancel. So the
// stack's circulations are summed per strip first, and the wake trace is a
// row of 2-D point vortices at the strip edges, each carrying the difference
// between neighbouring strip circulations.
//
// The plane is taken normal to +x, which is where the frozen wake actually
// goes in this solver. Strictly the Trefftz plane is normal to the
// FREESTREAM, and the two differ at incidence -- but the wake model is frozen
// along +x, so a plane normal to the freestream would be inconsistent with
// the vorticity it is supposed to intercept. Matching the plane to the wake
// keeps the two consistent; the error is the wake model's, not this
// integral's, and it is the same O(alpha^2) the frozen wake already carries.
//
// --- what this does NOT fix ------------------------------------------------------
// Profile drag, still absent (the method is inviscid -- see
// Aeolion::DragEstimate). And a body that genuinely produces drag through its
// own trailing vorticity -- a lifting fuselage at incidence shedding vortices
// -- is not represented here either, because a source panelling cannot shed
// vorticity in the first place. This integral reports the induced drag of the
// LIFTING system, which is exactly the quantity whose span efficiency is
// bounded by one.

#pragma once

#include "Aeolion/Lattice/Panel.h"
#include "Aeolion/Math/Constants.h"
#include "Aeolion/Math/Vec3.h"
#include "Aeolion/Solver/ReferenceGeometry.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <numbers>
#include <vector>

namespace Aeolion::Solver {

/** Squared distance below which two wake filaments are the same filament. */
inline constexpr double TrefftzMergeTolerance = 1e-12;
/** Regularizes a vortex's own singularity when it lands on an evaluation point. */
inline constexpr double TrefftzCoreSquared = 1e-14;

/**
 * One trailing filament where it pierces the Trefftz plane: a 2-D point
 * vortex in the (y, z) crossflow plane.
 */
struct WakeTrace {
    double y = 0.0, z = 0.0;
    double Gamma = 0.0; ///< Circulation of the filament [m^2/s].
};

/**
 * The wake trace of a lifting system: strip circulations differenced into
 * shed filaments.
 *
 * Panels sharing a StripIndex are one chordwise stack covering a single
 * spanwise station, and the interior trailing legs of such a stack cancel --
 * only the net strip circulation reaches the far wake. Panels with no strip
 * index are treated as their own strip, which is what a single-row lattice
 * (the wing's usual case) reduces to anyway.
 */
[[nodiscard]] inline std::vector<WakeTrace> BuildWakeTrace(const std::vector<Panel>& panels,
                                                            const std::vector<double>& gamma) {
    std::vector<WakeTrace> trace;
    if (panels.size() != gamma.size() || panels.empty()) return trace;

    // Sum circulation per strip, keeping each strip's own edge positions.
    struct Strip {
        double Gamma = 0.0;
        Vec3 A{0, 0, 0}, B{0, 0, 0};
        bool Seeded = false;
    };
    std::map<int, Strip> strips;
    int synthetic = -1;
    for (std::size_t i = 0; i < panels.size(); ++i) {
        const int key = (panels[i].StripIndex >= 0) ? panels[i].StripIndex : synthetic--;
        Strip& strip = strips[key];
        strip.Gamma += gamma[i];
        if (!strip.Seeded) { // the stack's edges are shared; the first panel states them
            strip.A = panels[i].A;
            strip.B = panels[i].B;
            strip.Seeded = true;
        }
    }

    // Each strip sheds +Gamma at one edge and -Gamma at the other. Filaments
    // landing on the same point are merged, which is what makes adjacent
    // strips' shed vorticity cancel down to the DIFFERENCE in circulation --
    // the physical trailing sheet, not one filament per panel edge.
    const auto shed = [&trace](const Vec3& point, double value) {
        for (WakeTrace& existing : trace) {
            const double dy = existing.y - point.y, dz = existing.z - point.z;
            if (dy * dy + dz * dz < TrefftzMergeTolerance) {
                existing.Gamma += value;
                return;
            }
        }
        trace.push_back({point.y, point.z, value});
    };
    // SIGNS, from the horseshoe circuit rather than from intuition. Solver.h
    // traverses inf -> A -> B -> inf, so the filament at A is traversed in
    // the MINUS x direction and the one at B in plus x. A filament's
    // circulation vector points along its traversal, and the Trefftz plane
    // looks downstream (+x), so A pierces it as -Gamma and B as +Gamma.
    //
    // Worth deriving rather than guessing: getting this backwards leaves
    // every magnitude exactly right -- elliptic loading still returns
    // |e| = 1 -- and flips the drag negative, which looks like the very
    // defect this file exists to fix.
    for (const auto& [key, strip] : strips) {
        shed(strip.A, -strip.Gamma);
        shed(strip.B, +strip.Gamma);
    }

    std::erase_if(trace, [](const WakeTrace& w) { return std::fabs(w.Gamma) < 1e-14; });
    return trace;
}

/**
 * Downwash (the z-component of the crossflow velocity) induced in the
 * Trefftz plane by the whole wake trace, at a point in it.
 *
 * The wake there is DOUBLY infinite -- it runs from minus to plus infinity
 * through the plane -- so each filament acts as a full 2-D point vortex,
 * inducing twice what the corresponding semi-infinite leg induces at the
 * lifting line. That factor of two is the difference between a far-field
 * calculation and a rearranged near-field one.
 */
[[nodiscard]] inline double TrefftzDownwash(const std::vector<WakeTrace>& trace, double y, double z) {
    double w = 0.0;
    for (const WakeTrace& filament : trace) {
        const double dy = y - filament.y, dz = z - filament.z;
        const double r2 = dy * dy + dz * dz;
        if (r2 < TrefftzCoreSquared) continue; // its own singularity
        w += filament.Gamma * dy / (Math::Two * std::numbers::pi * r2);
    }
    return w;
}

/** What a far-field induced-drag evaluation produces. */
struct TrefftzResult {
    bool Valid = false;
    double Di = 0.0;   ///< Induced drag [N].
    double CDi = 0.0;  ///< Referred to ReferenceGeometry::Area.
    /**
     * Span efficiency, CL^2 / (pi AR CDi). Bounded above by ONE for any
     * planar lifting system -- precisely the bound the near-field result
     * violates on a coupled configuration, so it is reported here as the
     * headline diagnostic rather than left for a caller to form.
     *
     * WHICH CL. This drag is the LIFTING SYSTEM's, because only the lifting
     * system sheds a wake, so the CL passed in must be the lifting system's
     * too. Feeding the whole configuration's CL -- wing plus the body's
     * share -- compares a configuration lift against a wing-only drag and
     * inflates the efficiency past one for a perfectly sound integral.
     *
     * On the airframe in tests/Data the fuselage carries about 9% of the
     * lift, and the mistake is worth exactly that: e = 1.17 against the
     * configuration CL, e = 0.97 against the wing's own. Use
     * SolveResult::LiftBySurface to get the latter.
     */
    double SpanEfficiency = 0.0;
    std::size_t Filaments = 0;
};

/**
 * Induced drag of a lifting system by far-field (Trefftz-plane) integration.
 *
 * Takes ONLY the lifting panels and their circulations: source panels are
 * not a parameter because a closed body sheds no wake and therefore cannot
 * contribute. That is the whole point -- the body's near-field force
 * residual is excluded by construction rather than subtracted.
 *
 * `CL` is used solely to report span efficiency; it does not enter the drag.
 */
[[nodiscard]] inline TrefftzResult TrefftzInducedDrag(const std::vector<Panel>& panels,
                                                       const std::vector<double>& gamma,
                                                       const ReferenceGeometry& ref, double rho,
                                                       double Vinf, double CL = 0.0) {
    TrefftzResult result;
    const std::vector<WakeTrace> trace = BuildWakeTrace(panels, gamma);
    if (trace.empty() || !(rho > 0.0) || !(Vinf > 0.0) || !(ref.Area > 0.0)) return result;

    result.Valid = true;
    result.Filaments = trace.size();

    // D_i = -(rho/2) * sum over strips of Gamma * w_T * width. Evaluating the
    // downwash at the strip MIDPOINT rather than at a filament keeps the
    // evaluation off the singularities, which is why the trace and the
    // strips are walked separately here.
    std::map<int, std::pair<Vec3, Vec3>> edges;
    std::map<int, double> stripGamma;
    int synthetic = -1;
    for (std::size_t i = 0; i < panels.size() && i < gamma.size(); ++i) {
        const int key = (panels[i].StripIndex >= 0) ? panels[i].StripIndex : synthetic--;
        stripGamma[key] += gamma[i];
        edges.try_emplace(key, panels[i].A, panels[i].B);
    }

    double drag = 0.0;
    for (const auto& [key, circulation] : stripGamma) {
        const auto& [A, B] = edges.at(key);
        const double yMid = Math::Half * (A.y + B.y), zMid = Math::Half * (A.z + B.z);
        const double width = std::hypot(B.y - A.y, B.z - A.z);
        if (!(width > 0.0)) continue;
        drag += circulation * TrefftzDownwash(trace, yMid, zMid) * width;
    }
    result.Di = -Math::Half * rho * drag;

    const double q = Math::Half * rho * Vinf * Vinf;
    result.CDi = result.Di / (q * ref.Area);
    const double aspect = (ref.Area > 0.0) ? ref.Span * ref.Span / ref.Area : 0.0;
    result.SpanEfficiency = (std::fabs(result.CDi) > 1e-12 && aspect > 0.0)
                                ? CL * CL / (std::numbers::pi * aspect * result.CDi)
                                : 0.0;
    return result;
}

} // namespace Aeolion::Solver
