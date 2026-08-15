// PropulsionMapExport.cpp -- the propulsion side of the DAVE-ML flight
// model (models/README.md): the CONTRACT's ducted propulsor -- rotor from
// the propulsion_bemt block via Geometry::ToPropeller, duct shroud from
// the contract's duct block, the four DuctJet vanes -- swept over advance
// ratio J and the three vane command modes, each condition solved by the
// two-way rotor-vane coupling (cascade closure, the production default;
// the same method the SciTech paper's reference case validates).
//
// AXIAL INFLOW ONLY. The rotor-vane machinery is built on axisymmetry
// end to end (slipstream bands, azimuthal-mean vane feedback,
// AxialInflowFromBands), so disk incidence is not representable here:
// alphaDisk carries no axis in this map, and the DAVE-ML spec records the
// nonaxial model as blocked future work. Fabricating an alphaDisk sweep
// from this solver would be an artifact, not a result.
//
// ROTATION SENSE. The validated solve path poses the blades and spins
// them with fc.p = +Omega about solver +x. Under the contract->solver
// frame flip that is a rotation VECTOR of -Omega * x_frd (clockwise seen
// from the nose). The contract states the rotation axis but not the
// handedness of reference_rpm about it, so the map is generated in the
// solver's validated sense and states it in meta ("rotationVectorFrd");
// if the aircraft's rotor spins opposite, regenerate mirrored rather
// than silently flipping signs downstream (models/README.md open item).
//
// FRAMES. Contract ingest and wrench export are the only two conversion
// sites (the 180-degree rotation about y; ADR-0016, and the moment-arm
// bug note in AttachmentSweepExport.cpp). Everything between them lives
// in the solver frame; everything in the JSON is FRD. Moments are about
// the contract's moment_reference_point, converted once at ingest
// relative to the duct center, which is where BuildPropellerDuct poses
// the rotor plane (origin, mid-chord).
//
// VANE COMMAND MODES (models/README.md): with all four hinge axes
// pointing radially OUTWARD, roll is the common mode, pitch the y-pair
// differential, yaw the z-pair differential. Deflections are right-hand
// rule about the contract's own FRD axes -- BuildDuctVanes converts the
// axis internally, and a proper rotation preserves the rule. The numeric
// sign meaning is pinned by the DAVE-ML checkData shots, not by comments.
//
// Both coupling sides run the analytic section polar -- the configuration
// every rotor-vane suite validates (the helical jet's swirl angles sit
// beyond the BL section model's convergent envelope; see the note in
// TestRotorVaneCoupling.cpp). A per-section BL rotor polar is a recorded
// refinement, not this map.
//
// Rows are flushed to the output JSON as they finish (the
// ParticleCrosscheckExport precedent), each carrying its own Converged
// flag and residual -- an unconverged row is exported as such, never
// laundered.
//
// Usage:
//   aeolion_propulsion_map <handoff.json> <out.json> [rpm] [modes] [jMax]
//                          [jSteps] [maxOuter]
//     rpm      default: the contract's reference_rpm
//     modes    "all" (default) | "baseline" | "pitch" | "yaw" | "roll"
//     jMax     default 1.0
//     jSteps   default 11 (J = 0 .. jMax inclusive)
//     maxOuter outer rotor-vane passes; default 24 (the solver's own
//              DefaultOuterIterations).
//     relax    vane-feedback under-relaxation; default 0.35 (the solver's
//              DefaultOuterRelaxation).
//
// MEASURED, so it is not re-litigated. At POSITIVE common-mode (roll)
// deflection and J = 0.6, three conditions do not converge, and NEITHER
// standard cure works:
//
//   - Not a budget problem. Residuals at 24 and at 48 passes agree to
//     three digits (+5 deg: 0.0656 both; +10: 0.0459 -> 0.0465; +15:
//     0.0268 -> 0.0278 -- two marginally WORSE) against a 0.01
//     tolerance. A stationary limit cycle, not slow convergence.
//   - Not a feedback-gain problem. Damping the vane-feedback relaxation
//     0.35 -> 0.15 still fails at 40 passes.
//
// The failure is SIGN-ASYMMETRIC: negative roll converges in 2-3 passes
// at the same J. That asymmetry is the physics -- a common-mode
// deflection either adds to or opposes the residual swirl the vanes sit
// in, and only one sign drives them toward the incidence where the vane
// row is multivalued (the bistability RotorVaneCoupling.h's VaneClosure
// comment describes).
//
// What the two independent iteration paths DO give is an uncertainty on
// the exported value, which is the useful deliverable: comparing
// relaxation 0.35 against 0.15, thrust agrees to 0.004% at +5 deg and
// 0.06% at +15 deg, but differs by 3.5% at +10 deg. So two of the three
// cells are reproducible despite missing tolerance, and one is genuinely
// path-dependent. The assembler must carry these as DAVE-ML uncertainty
// bounds rather than as ordinary table entries.

#include "Aeolion/Geometry/HandoffContract.h"
#include "Aeolion/PanelBuilder/PanelBuilder.h"
#include "Aeolion/Solver/RotorVaneCoupling.h"
#include "Aeolion/Solver/Solver.h"
#include "Aeolion/Solver/ViscousCoupling.h"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <numbers>
#include <string>
#include <utility>
#include <vector>

using namespace Aeolion;
namespace PB = Aeolion::PanelBuilder;
namespace S = Aeolion::Solver;

namespace {

constexpr double Rho = 1.225;
constexpr double HoverFloorSpeed = 0.01; // the house floor for fc.Vinf at J = 0

// The vane-mode deltas, within the contract's +-15 deg soft limit
// (models/README.md vaneBp). Zero lives in the baseline rows.
const double ModeDeltasDeg[] = {-15.0, -10.0, -5.0, 5.0, 10.0, 15.0};

// --- the two frame-conversion sites -----------------------------------------
// Contract FRD (x-fwd/z-down) <-> solver (x-aft/z-up): 180 degrees about y,
// a proper rotation, applied componentwise to points, forces and moments
// alike. These two helpers are the ONLY frame logic in this driver.
S::Vec3 IngestFrd(const Math::Vec3& v) { return S::Vec3(-v.x, v.y, -v.z); }

struct WrenchFrd {
    double Fx = 0.0, Fy = 0.0, Fz = 0.0; // [N]
    double Mx = 0.0, My = 0.0, Mz = 0.0; // [N*m], about the contract's moment reference point
};

// Solver-frame loads (Di streamwise +x_vlm, Y +y, L +z_vlm; moments about
// RefPoint in the same axes) exported to FRD.
WrenchFrd ExportWrench(const S::SolveResult& base) {
    WrenchFrd w;
    w.Fx = -base.Di;
    w.Fy = base.Y;
    w.Fz = -base.L;
    w.Mx = -base.Mx;
    w.My = base.My;
    w.Mz = -base.Mz;
    return w;
}

WrenchFrd operator+(const WrenchFrd& a, const WrenchFrd& b) {
    return {a.Fx + b.Fx, a.Fy + b.Fy, a.Fz + b.Fz, a.Mx + b.Mx, a.My + b.My, a.Mz + b.Mz};
}

// --- vane classification -----------------------------------------------------
// The cruciform's four DuctJet surfaces, identified by their OUTWARD radial
// hinge axes in the contract's FRD frame (z down): [0,0,1] extends down,
// [0,0,-1] up, [0,1,0] right, [0,-1,0] left. Indices are into the full
// ControlSurfaces list, which is also how BuildDuctVanes' deflection
// vector is addressed (viewer precedent). Classified by axis, never by
// name or position -- contract vane names collide by design.
struct VaneIndices {
    int Bottom = -1, Left = -1, Top = -1, Right = -1;
};

bool ClassifyVanes(const std::vector<Geometry::ControlSurface>& surfaces, VaneIndices& vanes) {
    for (std::size_t i = 0; i < surfaces.size(); ++i) {
        const Geometry::ControlSurface& s = surfaces[i];
        if (s.Binding != Geometry::ControlSurfaceBinding::DuctJet) continue;
        const Math::Vec3& a = s.HingeAxis;
        int* slot = nullptr;
        if (a.z > 0.9)
            slot = &vanes.Bottom;
        else if (a.z < -0.9)
            slot = &vanes.Top;
        else if (a.y > 0.9)
            slot = &vanes.Right;
        else if (a.y < -0.9)
            slot = &vanes.Left;
        if (!slot || *slot >= 0) return false; // oblique axis, or a duplicate
        *slot = static_cast<int>(i);
    }
    return vanes.Bottom >= 0 && vanes.Left >= 0 && vanes.Top >= 0 && vanes.Right >= 0;
}

// Mode mixing (models/README.md): roll = common mode on all four outward
// axes; pitch = y-pair differential (right +, left -); yaw = z-pair
// differential (bottom +, top -). This is the mixing matrix of the
// report's Eq. (mixing), applied as one linear map so a combined command
// is expressed exactly as the allocator would issue it -- which is what
// makes the superposition measurement meaningful rather than circular.
std::vector<double> MixDeflections(double pitchDeg, double yawDeg, double rollDeg,
                                   const VaneIndices& vanes, std::size_t surfaceCount) {
    std::vector<double> deflections(surfaceCount, 0.0);
    deflections[vanes.Bottom] = rollDeg + yawDeg;
    deflections[vanes.Left] = rollDeg - pitchDeg;
    deflections[vanes.Top] = rollDeg - yawDeg;
    deflections[vanes.Right] = rollDeg + pitchDeg;
    return deflections;
}

/** One swept condition: a labelled point in the three-mode command space. */
struct Command {
    std::string Mode; ///< "baseline", "pitch", "yaw", "roll", a pair label, or "vane<n>".
    double DeltaDeg = 0.0; ///< The commanded magnitude (both members, for a pair).
    double Pitch = 0.0, Yaw = 0.0, Roll = 0.0;
    /**
     * Per-vane angles applied directly, bypassing the mode mixing. Used by
     * the "single" sweep, which deflects ONE vane at a time to test
     * whether the wrench is a sum of independent per-vane responses --
     * the structure the mode-superposition measurement points to.
     * Empty means the (Pitch, Yaw, Roll) mixing applies.
     */
    std::vector<double> Explicit;
};

void WriteWrench(std::ofstream& out, const char* key, const WrenchFrd& w) {
    out << '"' << key << R"(":{"fx":)" << w.Fx << ",\"fy\":" << w.Fy << ",\"fz\":" << w.Fz
        << ",\"mx\":" << w.Mx << ",\"my\":" << w.My << ",\"mz\":" << w.Mz << '}';
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "usage: aeolion_propulsion_map <handoff.json> <out.json>"
                     " [rpm] [modes] [jMax] [jSteps]\n";
        return 1;
    }
    const std::string handoffPath = argv[1];
    const std::string outPath = argv[2];

    Geometry::HandoffContract contract;
    try {
        contract = Geometry::LoadHandoff(handoffPath);
    } catch (const std::exception& error) {
        std::cerr << "cannot load " << handoffPath << ": " << error.what() << "\n";
        return 1;
    }
    if (!contract.Duct.IsStated) {
        std::cerr << "contract states no duct block; the propulsor map needs one\n";
        return 1;
    }
    if (!contract.MomentReferencePointStated) {
        std::cerr << "contract states no moment_reference_point; refusing to invent one\n";
        return 1;
    }

    const double rpm =
        (argc > 3) ? std::atof(argv[3]) : contract.Propulsion.ReferenceRpm;
    const std::string modeFilter = (argc > 4) ? argv[4] : "all";
    const double jMax = (argc > 5) ? std::atof(argv[5]) : 1.0;
    const int jSteps = (argc > 6) ? std::atoi(argv[6]) : 11;
    const int maxOuter = (argc > 7) ? std::atoi(argv[7]) : S::DefaultOuterIterations;
    const double relax = (argc > 8) ? std::atof(argv[8]) : S::DefaultOuterRelaxation;
    if (rpm <= 0.0 || jSteps < 1 || maxOuter < 2 || !(relax > 0.0 && relax <= 1.0)) {
        std::cerr << "rpm must be positive, jSteps >= 1, maxOuter >= 2, 0 < relax <= 1\n";
        return 1;
    }

    const auto prop = Geometry::ToPropeller(contract.Propulsion);
    const double omega = rpm * 2.0 * std::numbers::pi / 60.0;
    const double revsPerSec = rpm / 60.0;
    const double diameter = 2.0 * prop.Radius;

    const double shroudInner = 0.5 * contract.Duct.InnerDiameter;
    const double shroudOuter = 0.5 * contract.Duct.OuterDiameter;
    const double shroudChord = contract.Duct.Chord;

    VaneIndices vanes;
    if (!ClassifyVanes(contract.ControlSurfaces, vanes)) {
        std::cerr << "contract's DuctJet surfaces are not the expected cruciform"
                     " (four vanes, outward radial axes along +-y/+-z)\n";
        return 1;
    }

    // Ingest site: the moment reference point, taken relative to the duct
    // center because BuildPropellerDuct poses the rotor plane at the
    // origin. The flip happens HERE and at ExportWrench, nowhere else.
    const Math::Vec3 refFrd(contract.MomentReferencePoint.x - contract.Duct.Center.x,
                            contract.MomentReferencePoint.y - contract.Duct.Center.y,
                            contract.MomentReferencePoint.z - contract.Duct.Center.z);
    const S::Vec3 refPoint = IngestFrd(refFrd);

    const auto strips = PB::BuildPropellerStrips(prop);
    const auto duct = PB::BuildPropellerDuct(shroudInner, shroudOuter, shroudChord);
    const double trail = S::DefaultTrailSpanFactor * diameter;
    const double vaneTrail = S::DefaultTrailSpanFactor * 2.0 * shroudInner;

    S::ReferenceGeometry vaneRef;
    vaneRef.Area = 1.0; // dimensional readouts only
    vaneRef.Span = 2.0 * shroudInner;
    vaneRef.Chord = shroudChord;

    std::ofstream out(outPath);
    if (!out) {
        std::cerr << "cannot open " << outPath << " for writing\n";
        return 1;
    }
    out.precision(9);

    out << "{\n\"meta\":{\"designId\":\"" << contract.DesignId << "\",\"schema\":\""
        << contract.SchemaVersion << "\",\"rpm\":" << rpm << ",\"revsPerSec\":" << revsPerSec
        << ",\"omegaRadps\":" << omega << ",\"diameterM\":" << diameter
        << ",\"radiusM\":" << prop.Radius << ",\"hubRadiusM\":" << prop.HubRadius
        << ",\"bladeCount\":" << prop.BladeCount << ",\"shroudInnerM\":" << shroudInner
        << ",\"shroudOuterM\":" << shroudOuter << ",\"shroudChordM\":" << shroudChord
        << ",\"rho\":" << Rho << ",\"hoverFloorMps\":" << HoverFloorSpeed
        << ",\"closure\":\"cascade\",\"sectionModel\":\"analytic\""
        << ",\"maxOuterIterations\":" << maxOuter << ",\"relaxation\":" << relax
        << ",\"axialInflowOnly\":true,\"rotationVectorFrd\":[" << -omega << ",0,0]"
        << ",\"momentRefFrd\":[" << contract.MomentReferencePoint.x << ','
        << contract.MomentReferencePoint.y << ',' << contract.MomentReferencePoint.z
        << "],\"ductCenterFrd\":[" << contract.Duct.Center.x << ',' << contract.Duct.Center.y
        << ',' << contract.Duct.Center.z << "],\"tipMach\":" << omega * prop.Radius / 340.0
        << ",\"vaneSoftLimitDeg\":15,\"frames\":\"wrench components FRD; moments about"
           " moment_reference_point\"},\n\"rows\":[\n";
    out.flush();

    // The solve for one (J, deflection) condition. Each row is
    // independent -- the coupled solve carries its own internal
    // continuation (ascending swirl-budget seed, vane warm starts).
    const auto solveRow = [&](double j, const std::vector<double>& deflections) {
        const double axialSpeed = j * revsPerSec * diameter;

        S::FreestreamConditions fc;
        fc.Vinf = std::max(axialSpeed, HoverFloorSpeed);
        fc.alphaDeg = 0.0;
        fc.betaDeg = 0.0;
        fc.rho = Rho;
        fc.p = omega;
        fc.RefPoint = refPoint;

        const auto panels = PB::BuildPropellerLattice(prop, axialSpeed, omega);
        S::ReferenceGeometry ref;
        for (const auto& panel : panels) ref.Area += panel.PlanformArea;
        ref.Span = diameter;
        ref.Chord = ref.Area / ref.Span;

        S::FreestreamConditions vaneFc = fc;
        vaneFc.p = 0.0; // static frame

        const S::VaneBuilder vaneBuilder =
            [&](const std::function<S::Vec3(const S::Vec3&)>& flow)
            -> std::pair<std::vector<S::Panel>, std::vector<S::StripSection>> {
            return {PB::BuildDuctVanes(contract.ControlSurfaces, shroudInner, 0.5 * shroudChord,
                                       shroudChord, deflections, flow),
                    PB::BuildDuctVaneStrips(contract.ControlSurfaces, shroudInner, shroudChord,
                                            deflections)};
        };
        const S::RotorBuilder rotorBuilder = [&](const std::function<double(double)>& inflow) {
            return PB::BuildPropellerLattice(prop, axialSpeed, omega, inflow);
        };

        S::RotorVaneOptions options;
        options.MaxOuterIterations = maxOuter;
        options.Relaxation = relax;

        return S::SolveRotorVaneCoupled(panels, strips, fc, ref, trail,
                                        S::AnalyticSectionModel{}, duct, vaneBuilder, vaneFc,
                                        vaneRef, vaneTrail, S::AnalyticSectionModel{}, axialSpeed,
                                        Rho, options, rotorBuilder);
    };

    const bool wantAll = modeFilter == "all";
    std::vector<Command> conditions;
    if (wantAll || modeFilter == "baseline") conditions.push_back({"baseline", 0.0, 0, 0, 0});
    for (const double d : ModeDeltasDeg) {
        if (wantAll || modeFilter == "pitch") conditions.push_back({"pitch", d, d, 0, 0});
        if (wantAll || modeFilter == "yaw") conditions.push_back({"yaw", d, 0, d, 0});
        if (wantAll || modeFilter == "roll") conditions.push_back({"roll", d, 0, 0, d});
    }
    // Simultaneous commands, for the superposition measurement: each pair
    // at HALF the soft limit so the per-vane sum -- which is what the
    // mechanism actually sees, and what an allocator must saturate --
    // still lands inside +-15 deg. Comparing these against the sum of the
    // matching single-mode increments is the only way to learn whether
    // the mode tables may be added, which the buildup assumes.
    if (wantAll || modeFilter == "combo") {
        constexpr double H = 7.5;
        for (const double s : {-1.0, 1.0}) {
            // The single-mode references at the same magnitude. Without
            // these the pair rows have nothing to be compared against --
            // +-7.5 is not on the ModeDeltasDeg grid.
            conditions.push_back({"pitch", s * H, s * H, 0, 0});
            conditions.push_back({"yaw", s * H, 0, s * H, 0});
            conditions.push_back({"roll", s * H, 0, 0, s * H});
            conditions.push_back({"pitch+yaw", s * H, s * H, s * H, 0});
            conditions.push_back({"pitch+roll", s * H, s * H, 0, s * H});
            conditions.push_back({"yaw+roll", s * H, 0, s * H, s * H});
        }
    }
    // Single-vane sweep: the per-vane response, from which a summation
    // model would predict every command. Deltas cover the magnitudes the
    // combo rows put on individual vanes (0, +-7.5, +-15).
    if (modeFilter == "single") {
        const int order[4] = {vanes.Bottom, vanes.Left, vanes.Top, vanes.Right};
        const char* names[4] = {"vaneBottom", "vaneLeft", "vaneTop", "vaneRight"};
        for (int v = 0; v < 4; ++v) {
            for (const double d : {-15.0, -7.5, 7.5, 15.0}) {
                Command cmd;
                cmd.Mode = names[v];
                cmd.DeltaDeg = d;
                cmd.Explicit.assign(contract.ControlSurfaces.size(), 0.0);
                cmd.Explicit[order[v]] = d;
                conditions.push_back(std::move(cmd));
            }
        }
    }
    if (conditions.empty()) {
        std::cerr << "unknown mode filter '" << modeFilter
                  << "' (all | baseline | pitch | yaw | roll | combo | single)\n";
        return 1;
    }

    const double normForce = Rho * revsPerSec * revsPerSec * std::pow(diameter, 4);
    const double normMoment = normForce * diameter;

    int rowsWritten = 0, unconverged = 0;
    for (const Command& cmd : conditions) {
        const auto deflections =
            cmd.Explicit.empty() ? MixDeflections(cmd.Pitch, cmd.Yaw, cmd.Roll, vanes,
                                                  contract.ControlSurfaces.size())
                                 : cmd.Explicit;
        for (int step = 0; step < jSteps; ++step) {
            const double j = (jSteps == 1) ? 0.0 : jMax * step / (jSteps - 1);
            const auto start = std::chrono::steady_clock::now();
            const auto result = solveRow(j, deflections);
            const double seconds =
                std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

            const WrenchFrd rotor = ExportWrench(result.Rotor.Base);
            const WrenchFrd vane = ExportWrench(result.Vanes.Base);
            const WrenchFrd total = rotor + vane;
            // The jet's angular-momentum flux = shaft torque (solver-frame
            // diagnostic, the budget's own quantity) -- NOT the FRD Mx,
            // which already carries the vane recovery and the arm to the
            // reference point.
            const double jetTorque =
                -(result.Rotor.InducedMoment.x + result.Rotor.ProfileMoment.x);

            if (rowsWritten++) out << ",\n";
            out << R"( {"mode":")" << cmd.Mode << R"(","deltaDeg":)" << cmd.DeltaDeg
                << ",\"cmdPitchDeg\":" << cmd.Pitch << ",\"cmdYawDeg\":" << cmd.Yaw
                << ",\"cmdRollDeg\":" << cmd.Roll << ",\"J\":" << j
                << ",\"VinfMps\":" << j * revsPerSec * diameter
                << ",\"converged\":" << (result.Converged ? "true" : "false")
                << ",\"outer\":" << result.OuterIterations << ",\"residual\":" << result.Residual
                << ",\"swirlFactor\":" << result.SwirlFactor << ',';
            WriteWrench(out, "rotorDuctFrd", rotor);
            out << ',';
            WriteWrench(out, "vanesFrd", vane);
            out << ',';
            WriteWrench(out, "totalFrd", total);
            out << ",\"thrustN\":" << total.Fx << ",\"jetTorqueNm\":" << jetTorque
                << ",\"ct\":" << total.Fx / normForce << ",\"cq\":" << jetTorque / normMoment
                << ",\"seconds\":" << seconds << '}';
            out.flush(); // a killed run keeps every finished row

            if (!result.Converged) ++unconverged;
            std::cout << cmd.Mode << " delta=" << cmd.DeltaDeg << " J=" << j << ": T=" << total.Fx
                      << " N  Q=" << jetTorque << " N*m  s=" << result.SwirlFactor
                      << (result.Converged ? "" : "  UNCONVERGED")
                      << " (outer=" << result.OuterIterations << ", " << seconds << " s)"
                      << std::endl;
        }
    }

    out << "\n]}\n";
    std::cout << "wrote " << outPath << " (" << rowsWritten << " rows, " << unconverged
              << " unconverged)" << std::endl;
    return 0;
}
