// Solver/RotorVaneCoupling.h
//
// Two-way rotor-vane coupling: a block Gauss-Seidel iteration alternating
// the rotating-frame rotor(+duct) solve and the static-frame vane solve,
// with a momentum budget on the swirl. The full method -- the frame
// argument, the azimuthal-mean feedback field, the budget contraction, the
// convergence measures, and what remains outside the model -- is written
// up in doc/theory.rst ("Two-way rotor-vane coupling"); comments here only
// anchor the code to it.
//
// Per outer pass: the rotor solves with the vanes' azimuthal-mean induced
// field in its external-velocity hook (MeanInducedField); the slipstream
// is rebuilt from the fresh rotor loading; and the swirl budget is solved
// as a bracketed 1-D ROOT PROBLEM at frozen rotor state -- find the swirl
// factor s at which the vanes' recovered torque equals the jet's
// angular-momentum flux (the shaft torque), each evaluation remeshing the
// vanes in the rescaled slipstream (their trailing legs follow the local
// mean flow) and re-solving them. Freezing the rotor inside the root find
// is load-bearing: driving s with a feedback law while the rotor re-solves
// was observed both to ratchet past the fixed point onto the floor
// (under-recovering forever after) and, driven the other way, to walk the
// coupled pair onto a runaway branch. Only after the budget is satisfied
// does the pass advance the (under-relaxed) vane feedback into the rotor.
// Convergence is judged on the relative change of the total wrench (net
// thrust, Mx), of the vane circulation, and any remaining budget excess.

#pragma once

#include "Aeolion/Lattice/Panel.h"
#include "Aeolion/Lattice/SourcePanel.h"
#include "Aeolion/Math/Constants.h"
#include "Aeolion/Solver/VaneCascade.h"
#include "Aeolion/Solver/ViscousCoupling.h"

#include <cmath>
#include <cstddef>
#include <functional>
#include <utility>
#include <vector>

namespace Aeolion::Solver {

// --- outer-iteration defaults ------------------------------------------------
inline constexpr int    DefaultOuterIterations = 24;
inline constexpr double DefaultOuterTolerance = 1e-2;   // relative wrench / circulation change
inline constexpr double DefaultOuterRelaxation = 0.35;  // on the vane gammas feeding back
inline constexpr double MinSwirlFactor = 0.05;          // the budget never quite extinguishes the swirl
inline constexpr double SwirlBudgetTolerance = 0.05;    // relative band around recovered = jet torque
inline constexpr int    MaxSwirlBudgetEvaluations = 8;  // vane solves one pass's root find may spend

/**
 * The vane side of the coupling, rebuilt every outer pass: given the
 * current local mean flow (freestream + slipstream), return the vane
 * panels and their aligned strip sections.
 */
using VaneBuilder = std::function<std::pair<std::vector<Panel>, std::vector<StripSection>>(
    const std::function<Vec3(const Vec3&)>& localFlow)>;

/**
 * The rotor side of the self-consistent wake: given the previous pass's
 * solved axial inflow vi(r) (empty on the seed pass), return the blade
 * panels with their trailing-leg pitch posed on it
 * (BuildPropellerLattice's axialInflow parameter).
 */
using RotorBuilder =
    std::function<std::vector<Panel>(const std::function<double(double)>& axialInflow)>;

/**
 * How the vane row is solved each outer pass. Cascade -- the default --
 * is the annulus-by-annulus momentum closure (VaneCascade.h): forces
 * bounded by each sector's mass flow by construction, no swirl-budget
 * factor needed, deterministic and single-valued where the lattice is
 * bistable (a deflected vane deep in the hover swirl). Lattice is the
 * vortex-lattice strip solve with the swirl-budget root find -- the
 * validated choice for brisk attached-incidence jets, where mutual strip
 * induction is real signal rather than a runaway channel.
 */
enum class VaneClosure { Cascade, Lattice };

struct RotorVaneOptions {
    int MaxOuterIterations = DefaultOuterIterations;
    double Tolerance = DefaultOuterTolerance;
    double Relaxation = DefaultOuterRelaxation;
    int AzimuthSamples = DefaultAzimuthSamples;
    VaneClosure Closure = VaneClosure::Cascade;
    bool SwirlBudget = true; ///< Lattice closure only: enforce vane Mx <= jet torque via s.
    ViscousCouplingOptions RotorOptions;
    ViscousCouplingOptions VaneOptions;
};

struct RotorVaneResult {
    ViscousCoupledResult Rotor;  ///< Converged rotor(+duct) state of the final pass.
    ViscousCoupledResult Vanes;  ///< Converged vane state of the final pass.
    std::vector<Panel> BladePanels;           ///< The final pass's blade mesh (self-consistent wake).
    std::vector<Panel> VanePanels;            ///< The final pass's vane mesh...
    std::vector<StripSection> VaneStrips;     ///< ...and its aligned strips.
    double SwirlFactor = 1.0;    ///< The budget's converged swirl scale.
    int OuterIterations = 0;
    bool Converged = false;
    double Residual = 0.0;       ///< Last relative change (wrench + circulation).
};

/**
 * Alternate the rotor(+duct) and vane solves to joint convergence. The
 * rotor arguments mirror SolveViscousCoupled's; the vanes come from
 * `vaneBuilder`, are solved under `vaneFc` (the static frame: its rates
 * should be zero), and their loads complete the wrench. axialSpeed and
 * rho feed the slipstream reconstruction.
 */
[[nodiscard]] inline RotorVaneResult SolveRotorVaneCoupled(
    const std::vector<Panel>& bladePanels, const std::vector<StripSection>& bladeStrips,
    const FreestreamConditions& fc, const ReferenceGeometry& ref, double trail,
    const SectionModel& rotorModel, const std::vector<SourcePanel>& ductSources,
    const VaneBuilder& vaneBuilder, const FreestreamConditions& vaneFc,
    const ReferenceGeometry& vaneRef, double vaneTrail, const SectionModel& vaneModel,
    double axialSpeed, double rho, const RotorVaneOptions& options = {},
    const RotorBuilder& rotorBuilder = nullptr) {
    RotorVaneResult res;
    res.BladePanels = bladePanels;

    std::vector<double> feedbackGamma; // under-relaxed vane circulation driving the rotor
    std::function<double(double)> solvedInflow = {}; // vi(r) of the previous pass; empty = seed
    double previousThrust = 0.0, previousMx = 0.0;

    // The vane warm start (continuation seed across evaluations and outer
    // passes), and the budget's starting point: ASCENDING continuation.
    // The budget seeds at the FLOOR swirl -- the mildest, attached-most
    // vane state -- and lets the root find raise s, so branch selection
    // happens from the attached side (the alpha-continuation philosophy).
    // Seeding at full swirl instead lets a long first inner solve wander
    // into a huge-circulation limit cycle whose state the warm start then
    // faithfully follows forever: recovered stays several-fold over the
    // budget at ANY s and the coupled pair locks onto a runaway branch.
    std::vector<double> vaneWarmStart;
    const bool lattice = options.Closure == VaneClosure::Lattice;
    if (lattice && options.SwirlBudget && vaneBuilder) res.SwirlFactor = MinSwirlFactor;

    for (res.OuterIterations = 1; res.OuterIterations <= options.MaxOuterIterations;
         ++res.OuterIterations) {
        // 0. The self-consistent wake: remesh the blades with the previous
        // pass's solved vi(r) setting each trailing leg's helix pitch (the
        // seed pass uses the builder's fixed-lambda floor).
        if (rotorBuilder) res.BladePanels = rotorBuilder(solvedInflow);

        // 1. Rotor + duct, feeling the vanes' azimuthal-mean field.
        std::function<Vec3(const Vec3&)> vaneField = nullptr;
        if (!res.VanePanels.empty() && !feedbackGamma.empty())
            vaneField = MeanInducedField(res.VanePanels, feedbackGamma, vaneTrail,
                                         options.AzimuthSamples);
        res.Rotor = SolveViscousCoupled(res.BladePanels, bladeStrips, fc, ref, trail, rotorModel,
                                        options.RotorOptions, ductSources, vaneField);

        // 2. Band the fresh loading: the slipstream for the vanes, the
        // inflow for the next pass's wake pitch.
        solvedInflow = AxialInflowFromBands(ComputeSlipstreamBands(res.Rotor, axialSpeed, rho));

        // 3. One vane evaluation at swirl scale s: remesh (the trailing
        // legs follow the local mean flow, which depends on s) and
        // re-solve in the rescaled slipstream, warm-started from the last
        // vane circulation -- continuation, so a vane riding a post-stall
        // branch boundary follows ONE branch across evaluations and
        // passes instead of flipping under tiny input changes. Returns
        // the recovered torque. The rotor state is FROZEN across
        // evaluations.
        const auto solveVanesAt = [&](double s) {
            const auto scaled = SlipstreamField(res.Rotor, axialSpeed, rho, ductSources, s);
            const auto flow = [&](const Vec3& point) {
                return Vec3(axialSpeed, 0.0, 0.0) + scaled(point);
            };
            auto [vanePanels, vaneStrips] = vaneBuilder(flow);
            res.VanePanels = std::move(vanePanels);
            res.VaneStrips = std::move(vaneStrips);
            if (res.VanePanels.empty()) return 0.0;
            if (!lattice) {
                // The jet's true angular-momentum flux is the shaft
                // torque; the cascade normalizes the sampled swirl to it
                // up front (the reconstruction over-carries flux).
                const double flux = -(res.Rotor.InducedMoment.x + res.Rotor.ProfileMoment.x);
                res.Vanes = SolveVaneCascade(res.VanePanels, res.VaneStrips, vaneFc, flow,
                                             vaneModel, flux);
                return res.Vanes.Base.Mx;
            }
            ViscousCouplingOptions warmOptions = options.VaneOptions;
            warmOptions.InitialGamma = vaneWarmStart; // empty = lattice seed
            res.Vanes = SolveViscousCoupled(res.VanePanels, res.VaneStrips, vaneFc, vaneRef,
                                            vaneTrail, vaneModel, warmOptions, {}, scaled);
            vaneWarmStart = res.Vanes.Base.gamma;
            return res.Vanes.Base.Mx;
        };

        double gammaShift = 0.0, gammaScale = 1.0;
        double jetTorque = 0.0, recovered = 0.0;
        if (vaneBuilder) {
            jetTorque = -(res.Rotor.InducedMoment.x + res.Rotor.ProfileMoment.x);
            recovered = solveVanesAt(res.SwirlFactor); // warm start: last pass's factor
        }
        if (!vaneBuilder || res.VanePanels.empty()) {
            if (!rotorBuilder) {
                res.Converged = true; // nothing at all to iterate
                break;
            }
        } else {
            // 4. The swirl budget as a bracketed root find on s in
            // [MinSwirlFactor, 1]: the jet's angular-momentum flux is the
            // shaft torque (blade circulatory + profile moments; the
            // duct's pressure moment is axisymmetrically ~zero, carrying
            // no swirl), and recovered torque grows with the swirl the
            // vanes see, so seek recovered(s) = jetTorque. Untested
            // endpoints are tried before bisecting, so the common
            // settled-warm-start case costs one evaluation, a slack
            // budget accepts s = 1 (recovered(1) <= Q: nothing to scale
            // away, including a swirl-opposing negative recovery), and a
            // row over-recovering even at the floor keeps the floor --
            // where any excess still shows in the convergence residual
            // rather than being declared away.
            // The root find runs only through the FIRST HALF of the outer
            // budget: s is a closure constant, not a state variable, and
            // once its increments fall inside the inner solve's own noise
            // the re-bisection just remeshes jitter into the wrench.
            // Later passes hold s and let the rest of the state settle.
            if (lattice && options.SwirlBudget && jetTorque > Math::Tiny &&
                2 * res.OuterIterations <= options.MaxOuterIterations) {
                const double band = SwirlBudgetTolerance * jetTorque;
                double s = res.SwirlFactor;
                double sLo = MinSwirlFactor, sHi = 1.0;
                bool floorTested = s <= MinSwirlFactor + Math::Tiny;
                bool fullTested = s >= 1.0 - Math::Tiny;
                for (int evals = 1; std::fabs(recovered - jetTorque) > band &&
                                    evals < MaxSwirlBudgetEvaluations;
                     ++evals) {
                    if (recovered > jetTorque) {
                        sHi = s;
                        if (sHi <= MinSwirlFactor + Math::Tiny) break; // floor it is
                        s = floorTested ? Math::Half * (sLo + sHi) : MinSwirlFactor;
                    } else {
                        sLo = s;
                        if (sLo >= 1.0 - Math::Tiny) break; // slack: full swirl stands
                        s = fullTested ? Math::Half * (sLo + sHi) : 1.0;
                    }
                    floorTested = floorTested || s <= MinSwirlFactor + Math::Tiny;
                    fullTested = fullTested || s >= 1.0 - Math::Tiny;
                    recovered = solveVanesAt(s);
                }
                res.SwirlFactor = s;
            }

            // Safety valve on the continuation: a pass still recovering
            // more than twice the budget is on a runaway branch (a valid
            // solve never needs to exceed the flux by that much) -- drop
            // the warm start so the next evaluation reseeds from the
            // lattice instead of faithfully following the bad branch.
            if (lattice && jetTorque > Math::Tiny && recovered > 2.0 * jetTorque)
                vaneWarmStart.clear();

            // 5. Relax the feedback circulation.
            const std::vector<double>& fresh = res.Vanes.Base.gamma;
            gammaScale = Math::Tiny;
            if (feedbackGamma.size() != fresh.size()) {
                feedbackGamma = fresh;
                gammaShift = 1.0; // first sight of this mesh topology: not converged yet
                gammaScale = 1.0;
            } else {
                for (std::size_t i = 0; i < fresh.size(); ++i) {
                    gammaShift = std::max(gammaShift, std::fabs(fresh[i] - feedbackGamma[i]));
                    gammaScale = std::max(gammaScale, std::fabs(fresh[i]));
                    feedbackGamma[i] = (1.0 - options.Relaxation) * feedbackGamma[i] +
                                       options.Relaxation * fresh[i];
                }
            }
        }

        const double thrust = -(res.Rotor.Base.Di + res.Vanes.Base.Di);
        const double mx = res.Rotor.Base.Mx + res.Vanes.Base.Mx;
        const double wrenchScale = std::max({std::fabs(thrust), std::fabs(mx), Math::Tiny});
        // Excess only beyond the root find's own acceptance band -- inside
        // it the budget is satisfied by construction. The cascade closure
        // conserves angular momentum per sector, so it carries no excess
        // term at all.
        const double budgetExcess =
            (lattice && options.SwirlBudget && jetTorque > Math::Tiny)
                ? std::max((recovered - jetTorque) / jetTorque - SwirlBudgetTolerance, 0.0)
                : 0.0;
        res.Residual = std::max({std::fabs(thrust - previousThrust) / wrenchScale,
                                 std::fabs(mx - previousMx) / wrenchScale, gammaShift / gammaScale,
                                 budgetExcess});
        previousThrust = thrust;
        previousMx = mx;

        if (res.OuterIterations > 1 && res.Residual < options.Tolerance) {
            res.Converged = true;
            break;
        }
    }
    return res;
}

} // namespace Aeolion::Solver
