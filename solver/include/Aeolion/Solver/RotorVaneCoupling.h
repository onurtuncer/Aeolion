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
// is rebuilt from the fresh rotor loading with the current swirl budget
// factor; the vanes are REMESHED (their trailing legs follow the local
// mean flow, which just changed) and re-solved in it; and the budget
// factor contracts whenever the vanes extract more torque than the jet's
// angular-momentum flux -- the shaft torque -- delivers. Convergence is
// judged on the relative change of the total wrench (net thrust, Mx) and
// of the vane circulation.

#pragma once

#include "Aeolion/Lattice/Panel.h"
#include "Aeolion/Lattice/SourcePanel.h"
#include "Aeolion/Math/Constants.h"
#include "Aeolion/Solver/ViscousCoupling.h"

#include <cmath>
#include <cstddef>
#include <functional>
#include <utility>
#include <vector>

namespace Aeolion::Solver {

// --- outer-iteration defaults ------------------------------------------------
inline constexpr int    DefaultOuterIterations = 12;
inline constexpr double DefaultOuterTolerance = 1e-2;   // relative wrench / circulation change
inline constexpr double DefaultOuterRelaxation = 0.7;   // on the vane gammas feeding back
inline constexpr double MinSwirlFactor = 0.05;          // the budget never quite extinguishes the swirl

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

struct RotorVaneOptions {
    int MaxOuterIterations = DefaultOuterIterations;
    double Tolerance = DefaultOuterTolerance;
    double Relaxation = DefaultOuterRelaxation;
    int AzimuthSamples = DefaultAzimuthSamples;
    bool SwirlBudget = true; ///< Enforce vane Mx <= jet torque via the swirl factor.
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
        const auto slipstream =
            SlipstreamField(res.Rotor, axialSpeed, rho, ductSources, res.SwirlFactor);
        const auto localFlow = [&](const Vec3& point) {
            return Vec3(axialSpeed, 0.0, 0.0) + slipstream(point);
        };

        // 3. Remesh and re-solve the vanes in it (when there are vanes at
        // all -- the driver also serves the vane-less self-consistent-wake
        // iteration).
        double gammaShift = 0.0, gammaScale = 1.0;
        double jetTorque = 0.0, recovered = 0.0;
        if (vaneBuilder) {
            auto [vanePanels, vaneStrips] = vaneBuilder(localFlow);
            res.VanePanels = std::move(vanePanels);
            res.VaneStrips = std::move(vaneStrips);
        }
        if (!vaneBuilder || res.VanePanels.empty()) {
            if (!rotorBuilder) {
                res.Converged = true; // nothing at all to iterate
                break;
            }
        } else {
            res.Vanes = SolveViscousCoupled(res.VanePanels, res.VaneStrips, vaneFc, vaneRef,
                                            vaneTrail, vaneModel, options.VaneOptions, {},
                                            slipstream);

            // 4. The swirl budget: the jet's angular-momentum flux is the
            // shaft torque (blade circulatory + profile moments; the duct's
            // pressure moment is axisymmetrically ~zero, carrying no swirl).
            jetTorque = -(res.Rotor.InducedMoment.x + res.Rotor.ProfileMoment.x);
            recovered = res.Vanes.Base.Mx;
            if (options.SwirlBudget && jetTorque > Math::Tiny && recovered > jetTorque)
                res.SwirlFactor =
                    std::max(res.SwirlFactor * jetTorque / recovered, MinSwirlFactor);

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
        const double budgetExcess =
            (options.SwirlBudget && jetTorque > Math::Tiny)
                ? std::max(recovered - jetTorque, 0.0) / jetTorque
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
