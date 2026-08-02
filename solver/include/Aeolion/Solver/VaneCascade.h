// Solver/VaneCascade.h
//
// The cascade momentum closure for duct-jet vane loads: annulus-by-annulus
// angular-momentum bookkeeping in place of the vane vortex lattice, for
// the regime where isolated-strip lifting-surface theory has no stable
// answer -- vanes deep in a hover jet's swirl. Full write-up in
// doc/theory.rst ("The cascade momentum closure").
//
// The failure this replaces, measured before it was replaced: the vane
// lattice at hover is bistable once deflected. Each strip's circulation is
// individually capped, but twenty capped strips can mutually inflate one
// another's local velocities, so the FORCE runs away without any gamma
// exceeding its cap -- a single vane at 8 degrees was solved carrying more
// side force than the entire net thrust, on one branch, while the same
// command landed on a sane branch at a neighboring rpm. No stabilizer
// fixes a model whose answer is not unique.
//
// The closure: each vane strip owns one azimuthal sector of its radial
// annulus. The mass flow through the sector is what the jet delivers,
//
//    mdot = rho * u * (2 pi r / N_sector) * width,
//
// and the tangential force the strip may carry is exactly the angular
// momentum it removes from that flow, F_t = mdot * dw. Blade-element
// aerodynamics supplies the same force from the section polar at the
// PASSAGE-MEAN flow (axial u, tangential w - dw/2), so the closure is a
// scalar equation per strip,
//
//    mdot * dw = F_t(dw),
//
// monotone in dw and solved by bisection on a bracket bounded by the
// local dynamic head. Every force is therefore bounded by construction:
// a strip can never carry more momentum than its sector's flow brings it,
// a dead-air strip outside the jet (u ~ 0) carries ~nothing, and the
// undeflected cruciform's recovered torque cannot exceed the jet's
// angular-momentum flux. Deflection enters through the strip's own chord
// frame (the mesher bakes the command into ChordDir/LiftDir), so control
// response comes out of the same bounded bookkeeping.
//
// What the closure gives up, honestly: no mutual induction between strips
// (each sector is independent), no unsteadiness, and the passage-mean
// blade-element picture assumes the vane row fills its annulus the way a
// stator row does. Those are the right trades at hover; in a brisk
// attached-incidence jet the lattice remains available and validated
// (RotorVaneOptions::Closure).

#pragma once

#include "Aeolion/Math/Constants.h"
#include "Aeolion/Solver/Solver.h"
#include "Aeolion/Solver/ViscousCoupling.h"

#include <cmath>
#include <cstddef>
#include <functional>
#include <numbers>
#include <vector>

namespace Aeolion::Solver {

inline constexpr int CascadeBisectionSteps = 48; // ~1e-14 of the bracket

/**
 * Solve the vane row by the cascade momentum closure. `panels` and
 * `strips` are the same aligned pair the lattice path uses (panels carry
 * geometry and orientation; strips carry the section frames with any
 * deflection baked in); `localFlow` is the total onset field at a point
 * (freestream + slipstream). Loads are reported in the same axes and the
 * same ViscousCoupledResult shape as SolveViscousCoupled, and Base.gamma
 * carries the equivalent bound circulations so the azimuthal-mean
 * feedback field (MeanInducedField) works unchanged.
 */
/**
 * `swirlFluxBudget` (optional, N*m): the angular-momentum flux the jet
 * actually carries -- the shaft torque. The azimuthal-mean slipstream
 * reconstruction is not flux-consistent (its far-wake axial doubling is
 * applied without the matching streamtube swirl mapping, and was measured
 * delivering ~3x the shaft torque to the vane plane); when the sampled
 * flux over the vane sectors exceeds the budget, the sampled swirl is
 * scaled ONCE, up front, by the ratio -- a derived normalization, not an
 * iterated feedback -- so conservation holds before any strip is solved.
 * Negative (the default) trusts the field as given.
 */
[[nodiscard]] inline ViscousCoupledResult SolveVaneCascade(
    const std::vector<Panel>& panels, const std::vector<StripSection>& strips,
    const FreestreamConditions& fc, const std::function<Vec3(const Vec3&)>& localFlow,
    const SectionModel& model, double swirlFluxBudget = -1.0) {
    ViscousCoupledResult res;
    const std::size_t n = panels.size();
    res.Strips.resize(n);
    res.Base.gamma.assign(n, 0.0);
    res.Iterations = 1;
    res.Converged = true; // the closure is a direct bracketed solve
    res.MaxResidual = 0.0;
    if (n == 0) return res;

    // Each strip's azimuthal sector: strips sharing a radial band split
    // the annulus evenly (a cruciform's four vanes -> quarter sectors).
    std::vector<double> radius(n), sectorCount(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const Vec3 mid = (panels[i].A + panels[i].B) * Math::Half;
        radius[i] = std::hypot(mid.y, mid.z);
    }
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            if (std::fabs(radius[i] - radius[j]) <
                Math::Half * std::min(strips[i].Width, strips[j].Width) + Math::Tiny)
                sectorCount[i] += 1.0;

    // Wind-axis projections match SolveViscousCoupled's, so consumers
    // read this result like any other solve.
    const double alpha = Math::DegToRad(fc.alphaDeg);
    const double beta = Math::DegToRad(fc.betaDeg);
    const Vec3 Vinf(fc.Vinf * std::cos(alpha) * std::cos(beta), fc.Vinf * std::sin(beta),
                    fc.Vinf * std::sin(alpha) * std::cos(beta));
    const Vec3 dragDir = Vinf.Normalized();
    const Vec3 liftDir = Cross(dragDir, Vec3(0, 1, 0)).Normalized();
    const Vec3 sideDir = Cross(liftDir, dragDir).Normalized();

    // The flux normalization: integrate the SAMPLED angular-momentum flux
    // over the vane sectors first, and if it exceeds the stated budget,
    // scale the sampled swirl once so conservation holds by construction.
    double swirlScale = 1.0;
    if (swirlFluxBudget > 0.0) {
        double sampledFlux = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            const Vec3 mid = (panels[i].A + panels[i].B) * Math::Half;
            const double r = radius[i];
            if (r < Math::Tiny) continue;
            const Vec3 radial = Vec3(0.0, mid.y, mid.z) * (1.0 / r);
            const Vec3 tangent = Cross(Vec3(1.0, 0.0, 0.0), radial);
            const Vec3 flow = localFlow(mid);
            const double u = std::max(Dot(flow, Vec3(1.0, 0.0, 0.0)), 0.0);
            const double w = Dot(flow, tangent);
            sampledFlux += fc.rho * u *
                           (Math::Two * std::numbers::pi * r / sectorCount[i]) * strips[i].Width *
                           w * r;
        }
        if (sampledFlux > swirlFluxBudget)
            swirlScale = swirlFluxBudget / sampledFlux;
    }

    Vec3 totalForce(0, 0, 0);
    for (std::size_t i = 0; i < n; ++i) {
        const StripSection& strip = strips[i];
        StripState& state = res.Strips[i];
        const Vec3 mid = (panels[i].A + panels[i].B) * Math::Half;
        state.Mid = mid;

        const double r = radius[i];
        if (r < Math::Tiny) continue;
        const Vec3 radial = Vec3(0.0, mid.y, mid.z) * (1.0 / r);
        const Vec3 tangent = Cross(Vec3(1.0, 0.0, 0.0), radial); // +x-wise swirl direction

        const Vec3 flow = localFlow(mid);
        const double u = Dot(flow, Vec3(1.0, 0.0, 0.0));
        const double w = Dot(flow, tangent) * swirlScale;

        // The sector's mass flow. A strip in dead air (no axial feed --
        // outside the contracted jet, or reversed) carries no flow and
        // therefore no force: the jet-edge pathology closes itself.
        const double mdot = fc.rho * std::max(u, 0.0) *
                            (Math::Two * std::numbers::pi * r / sectorCount[i]) * strip.Width;

        // Blade-element force at passage-mean flow, as a function of the
        // swirl change dw. Positive F_t (along the swirl) is force ON THE
        // VANE, i.e. angular momentum removed from the jet.
        const Vec3 spanDir = Cross(strip.ChordDir, strip.LiftDir).Normalized();
        const auto evaluate = [&](double dw, Vec3* forceOut) {
            const Vec3 vm = Vec3(1.0, 0.0, 0.0) * u + tangent * (w - Math::Half * dw);
            const double vmag = vm.Norm();
            if (vmag < Math::Tiny) {
                if (forceOut) *forceOut = Vec3(0, 0, 0);
                state.alphaEffDeg = 0.0;
                return 0.0;
            }
            const Vec3 vhat = vm * (1.0 / vmag);
            const double alphaDeg =
                Math::RadToDeg(std::atan2(Dot(vm, strip.LiftDir), Dot(vm, strip.ChordDir)));
            const double Re = vmag * strip.Chord / SeaLevelKinematicViscosity;
            const double Ma = vmag / SeaLevelSpeedOfSound;
            const SectionCoefficients sect = model(strip, alphaDeg, Re, Ma);
            const double q = Math::Half * fc.rho * vmag * vmag;
            const Vec3 liftHat = Cross(spanDir, vhat).Normalized();
            const Vec3 force =
                (liftHat * sect.cl + vhat * sect.cd) * (q * strip.Chord * strip.Width);
            if (forceOut) *forceOut = force;
            state.alphaEffDeg = alphaDeg;
            state.cl = sect.cl;
            state.cd = sect.cd;
            state.Re = Re;
            state.Ma = Ma;
            state.Vrel = vmag;
            return Dot(force, tangent);
        };

        // The closure mdot*dw = F_t(dw), bisected on a bracket bounded by
        // the local dynamic head (the jet cannot be turned harder than
        // its own speed). R(dw) = mdot*dw - F_t is monotone increasing.
        double dw = 0.0;
        if (mdot > Math::Tiny) {
            const double vref = std::hypot(u, w) + Math::Tiny;
            double lo = w - vref, hi = w + vref;
            if (mdot * lo - evaluate(lo, nullptr) > 0.0) {
                dw = lo;
            } else if (mdot * hi - evaluate(hi, nullptr) < 0.0) {
                dw = hi;
            } else {
                for (int k = 0; k < CascadeBisectionSteps; ++k) {
                    dw = Math::Half * (lo + hi);
                    ((mdot * dw - evaluate(dw, nullptr) > 0.0) ? hi : lo) = dw;
                }
            }
        }

        Vec3 force(0, 0, 0);
        evaluate(dw, &force); // final consistent state at the solved dw
        if (mdot <= Math::Tiny) force = Vec3(0, 0, 0);
        state.Force = force;
        state.Residual = 0.0;

        // Equivalent bound circulation for the feedback field, via the
        // same Kutta-Joukowski projection inversion the lattice coupling
        // uses -- no hand-derived sign anywhere.
        const Vec3 vm = Vec3(1.0, 0.0, 0.0) * u + tangent * (w - Math::Half * dw);
        const Vec3 dl = panels[i].B - panels[i].A;
        const double denom = vm.Norm() * vm.Norm() * strip.Chord * strip.Width;
        const double k = (denom > Math::Tiny)
                             ? Math::Two * Dot(Cross(vm, dl), strip.LiftDir) / denom
                             : 0.0;
        res.Base.gamma[i] = (std::fabs(k) > Math::Tiny) ? state.cl / k : 0.0;

        totalForce = totalForce + force;
        const Vec3 arm = mid - fc.RefPoint;
        const Vec3 liftPart = force - vm.Normalized() * (Dot(force, vm.Normalized()));
        res.InducedMoment = res.InducedMoment + Cross(arm, liftPart);
        res.ProfileMoment = res.ProfileMoment + Cross(arm, force - liftPart);
    }

    res.Base.Di = Dot(totalForce, dragDir);
    res.Base.L = Dot(totalForce, liftDir);
    res.Base.Y = Dot(totalForce, sideDir);
    const Vec3 moment = res.InducedMoment + res.ProfileMoment;
    res.Base.Mx = moment.x;
    res.Base.My = moment.y;
    res.Base.Mz = moment.z;
    return res;
}

} // namespace Aeolion::Solver
