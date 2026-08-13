// AttachmentSweepExport.cpp -- run the Journal of Aircraft paper's
// application case: an alpha/beta MATRIX of the coupled wing--body solve on
// the real geometry handoff, exporting, per attitude,
//
//   - the six force/moment coefficients and the stability derivatives,
//     computed on the COUPLED system (Solver::ComputeDerivatives cannot:
//     it takes a wing-only panel list, and the fuselage is exactly what
//     Cn_beta and Cm_alpha need to see);
//   - the per-strip attachment line;
//   - the separation survey -- each strip marched from its real attachment
//     point to turbulent separation (Solver/AttachmentBoundaryLayer.h);
//   - the fuselage skin-flow topology,
//
// as JSON for the figure and table renderer in
// papers/journal-of-aircraft/figures/.
//
// Two coupled systems are solved per flight condition, sharing the same
// fuselage+duct source panels:
//
//   carry -- TrimWingAtBody + CarryThroughLift, the physical solve. Its
//            forces and its fuselage surface field are what the paper
//            reports. But its two innermost strips extend their bound
//            segments to the centreline (that extension IS the
//            carry-through vortex), which puts their segment midpoints --
//            and therefore their reconstructed leading-edge points --
//            inside the body, where an attachment station is meaningless
//            and where the leading-edge direction differencing of the
//            NEIGHBOURING stations is corrupted too.
//
//   clean -- TrimWingAtBody only. Every exposed strip's geometry is true,
//            so the attachment line is computed on this one. The cost is
//            the root circulation the carry-through restores; exporting
//            BOTH lines lets the renderer (and the paper) show that the
//            outboard stations agree and only the near-root stations feel
//            the difference.
//
// The mesh is the contract's own, except chordwise_panels is forced to 1:
// ComputeAttachmentLine reads one local velocity per strip off the bound
// midpoints, the same one-Weissinger-row-per-strip contract
// SolveViscousCoupled states (see TODO item on multi-row lattices).

#include "Aeolion/Geometry/HandoffContract.h"
#include "Aeolion/PanelBuilder/PanelBuilder.h"
#include "Aeolion/Solver/AttachmentBoundaryLayer.h"
#include "Aeolion/Solver/AttachmentLine.h"
#include "Aeolion/Solver/DiskInduction.h"
#include "Aeolion/Solver/Solver.h"
#include "Aeolion/Solver/SurfaceFlow.h"
#include "Aeolion/Solver/TrefftzPlane.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <numbers>
#include <span>
#include <string>
#include <vector>

using namespace Aeolion;
namespace PB = Aeolion::PanelBuilder;
namespace S = Aeolion::Solver;

namespace {

constexpr double FlightSpeed = 25.0; // [m/s] -- VBAT-class cruise, same as TestAirframe
constexpr double Rho = 1.225;
constexpr double TrailSpans = 50.0;

// The condition whose fuselage streamlines and surface Cp become the
// skin-flow figure: combined incidence, so the windward attachment node
// and the leeward flow are both off the symmetry plane.
constexpr double FigureAlphaDeg = 8.0;
constexpr double FigureBetaDeg = 10.0;

constexpr int StreamlinePointBudget = 400; // per exported line, decimated

// Consumer-side body discretization. The azimuthal count is the same knob
// LatticeOptions already exposes; the nose refinement below is its axial
// counterpart.
//
// MEASURED WARNING: do not raise the sector count without checking the
// clean system's pivot ratio. The trimmed-without-carry-through wing's
// inner trailing legs run along the body flank at exactly the trim
// radius, and 24 sectors moves the flank control points to within ~6 mm
// of that line: the clean solve's MinPivotRatio drops 0.0077 -> 0.0015
// and its CL inflates by up to 80% (1.75 vs the carry solve's 0.96 at
// alpha=8). The carry solve is immune -- its root segment runs at the
// centreline, far from any source panel. 16 sectors keeps the flank
// control points 9.6 mm off the leg line and both solves consistent.
constexpr int BodySectors = 16;
constexpr int NoseRefineStations = 20;
constexpr double NoseRefineFraction = 0.10; // of body length

// The attitude matrix. Wider and finer than the figures need, because the
// coefficient/derivative table is the deliverable and the separation
// boundary has to be bracketed rather than assumed -- so alpha runs well
// past where the wing is expected to let go.
const double Alphas[] = {-4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0};
const double Betas[] = {-10.0, -5.0, 0.0, 5.0, 10.0};

// With the fan running the study is an alpha x THRUST matrix rather than
// alpha x beta: the unpowered sweep already established that the separation
// boundary is insensitive to sideslip, so spending solves on beta again
// would buy nothing. Sideslip collapses to the symmetry plane.
const double PoweredBetas[] = {0.0};

// Conditions the figures are drawn at (a subset of the matrix above).
const double FigureAlphas[] = {0.0, 4.0, 8.0};
const double FigureBetas[] = {-10.0, 0.0, 10.0};

// Resample the body's station list with cosine clustering toward the nose,
// keeping every original breakpoint so added stations lie exactly on the
// contract's own piecewise-linear radius law -- the geometry is unchanged,
// only its discretization is refined. This is the fix SurfaceFlow.h's
// nose-resolution caveat calls for: at moderate incidence the stagnation
// point sits a few millimetres from the apex, inside the first panel ring
// of the contract's 25-station list, where the honest answer is
// AttachesUpstream rather than a located node. Station spacing is a
// consumer choice exactly like the azimuthal sector count: the handoff
// states the shape, not the mesh.
void RefineNoseStations(Geometry::BodyGeometry& body) {
    if (body.Stations.size() < 2) return;
    const double xNose = body.Stations.front().x;
    const double xTail = body.Stations.back().x;
    const double depth = NoseRefineFraction * (xNose - xTail);

    std::vector<double> xs;
    xs.reserve(body.Stations.size() + NoseRefineStations);
    for (const Geometry::BodyStation& station : body.Stations) xs.push_back(station.x);
    for (int k = 1; k <= NoseRefineStations; ++k) {
        const double s =
            1.0 - std::cos(0.5 * std::numbers::pi * static_cast<double>(k) / NoseRefineStations);
        xs.push_back(xNose - depth * s);
    }
    std::ranges::sort(xs, std::greater<>());
    const auto duplicates =
        std::ranges::unique(xs, [](double a, double b) { return std::fabs(a - b) < 1e-9; });
    xs.erase(duplicates.begin(), duplicates.end());

    std::vector<Geometry::BodyStation> refined;
    refined.reserve(xs.size());
    for (const double x : xs) refined.push_back({x, Geometry::RadiusAt(body, x)});
    body.Stations = std::move(refined);
}

// One StripSection per single-row panel. The frame is the TRUE CHORD
// frame -- StripSection's own contract ("ChordDir: leading edge ->
// trailing edge") -- NOT the panel's: the quarter-to-three-quarter-chord
// direction follows the camber line and the panel Normal carries the
// camber slope at the control point, so a frame built from them is
// pitched by roughly that slope (~4 deg on this section). The incidence
// ComputeAttachmentLine resolves in that frame feeds SolveSectionContour,
// whose contour ALREADY carries the camber -- the tilt double-counts it,
// biasing alpha_n high by ~4 deg and with it the separation table and the
// attachment-line Reynolds margins. (Found via the post-stall driver,
// where the same tilted reconstruction shifted the coupled zero-lift to
// the sum of the lattice's and the thin-airfoil camber angles.) The
// handoff wing is rectangular, unswept and untwisted, so its true chord
// frame is the solver frame itself; a swept or twisted wing needs a
// PanelBuilder-side strip builder that carries the section plane.
std::vector<S::StripSection> StripsFromPanels(const std::vector<S::Panel>& panels, double halfSpan) {
    std::vector<S::StripSection> strips;
    strips.reserve(panels.size());
    for (const S::Panel& panel : panels) {
        const S::Vec3 mid = (panel.A + panel.B) * 0.5;
        S::StripSection strip;
        strip.ChordDir = S::Vec3(1.0, 0.0, 0.0);
        strip.LiftDir = S::Vec3(0.0, 0.0, 1.0);
        strip.Chord = (panel.SpanwiseWidth > 0.0) ? panel.PlanformArea / panel.SpanwiseWidth : 0.0;
        strip.Width = panel.SpanwiseWidth;
        strip.Eta = std::fabs(mid.y) / halfSpan;
        strips.push_back(strip);
    }
    return strips;
}

// A carry-through strip's bound segment starts AT the centreline.
bool IsCarryThroughStrip(const S::Panel& panel) {
    return std::min(std::fabs(panel.A.y), std::fabs(panel.B.y)) < 1e-9;
}

// The fan as an actuator disk, positioned from the contract.
//
// The blades sweep to propulsion_bemt.disk_radius inside the duct bore, and
// the tail boom occupies the middle, so the disk is an ANNULUS between the
// body radius at that station and the blade tip. Contract frame is x-forward
// and the solver x-aft, so the placement flips and the downstream axis is
// +x in solver axes -- this is a PUSHER installation, and the sign is the
// whole point: it puts the wing in the disk's upstream induction rather than
// in its wake.
S::ActuatorDisk MakeFanDisk(const Geometry::HandoffContract& contract, double thrust, double rho,
                            double axialSpeed) {
    S::ActuatorDisk disk;
    if (!(thrust > 0.0) || !contract.Duct.IsStated || !(contract.Propulsion.DiskRadius > 0.0))
        return disk;

    const double diskContractX = contract.Duct.Center.x;
    disk.Center = S::Vec3(-diskContractX, contract.Duct.Center.y, -contract.Duct.Center.z);
    disk.Axis = S::Vec3(1.0, 0.0, 0.0); // downstream, solver aft
    disk.Radius = contract.Propulsion.DiskRadius;
    disk.HubRadius = Geometry::RadiusAt(contract.Body, diskContractX);
    if (!(disk.HubRadius < disk.Radius)) disk.HubRadius = 0.0;
    disk.InducedVelocity =
        S::InducedVelocityFromThrust(thrust, rho, disk.Area(), axialSpeed);
    return disk;
}

// Central-difference derivatives about `base`, on the COUPLED system.
//
// Solver::ComputeDerivatives does exactly this but takes a wing-only panel
// list and prepares its own system, so it cannot see the fuselage. Since
// the whole point of the table is that the body contributes -- Cn_beta and
// Cm_alpha are where a fuselage shows up most -- the loop is written here
// against the already-factorized coupled system instead. The factorization
// depends only on geometry, so all ten perturbed solves reuse it.
S::StabilityDerivatives CoupledDerivatives(const S::PreparedSystem& prepared,
                                           const S::FreestreamConditions& base,
                                           const S::ReferenceGeometry& ref) {
    S::StabilityDerivatives d;
    const S::SolveResult s0 = S::SolveWithSystem(prepared, base, ref);
    d.CL0 = s0.CL; d.CDi0 = s0.CDi; d.CY0 = s0.CY;
    d.Cm0 = s0.Cm; d.Croll0 = s0.Croll; d.Cn0 = s0.Cn;

    const auto central = [&](S::FreestreamConditions plus, S::FreestreamConditions minus,
                             double step, auto&& assign) {
        const S::SolveResult rp = S::SolveWithSystem(prepared, plus, ref);
        const S::SolveResult rm = S::SolveWithSystem(prepared, minus, ref);
        assign(rp, rm, 2.0 * step);
    };

    const double dAlphaDeg = S::DefaultAlphaStepDeg, dBetaDeg = S::DefaultBetaStepDeg;
    {
        S::FreestreamConditions p = base, m = base;
        p.alphaDeg += dAlphaDeg;
        m.alphaDeg -= dAlphaDeg;
        central(p, m, Math::DegToRad(dAlphaDeg), [&](const auto& rp, const auto& rm, double h) {
            d.CL_alpha = (rp.CL - rm.CL) / h;
            d.CDi_alpha = (rp.CDi - rm.CDi) / h;
            d.Cm_alpha = (rp.Cm - rm.Cm) / h;
        });
    }
    {
        S::FreestreamConditions p = base, m = base;
        p.betaDeg += dBetaDeg;
        m.betaDeg -= dBetaDeg;
        central(p, m, Math::DegToRad(dBetaDeg), [&](const auto& rp, const auto& rm, double h) {
            d.CY_beta = (rp.CY - rm.CY) / h;
            d.Croll_beta = (rp.Croll - rm.Croll) / h;
            d.Cn_beta = (rp.Cn - rm.Cn) / h;
        });
    }

    const double cbar = ref.Chord, span = ref.Span;

    // Per-rad/s -> per reduced rate. The reduced rate is the rate TIMES
    // length/(2V), so converting divides by that, i.e. multiplies by
    // 2V/length. Written once, here, because having a second copy of this
    // expression is precisely how it came to be inverted in Solver.h.
    const auto reduce = [&](double perRadPerSecond, double length) {
        return perRadPerSecond * (2.0 * base.Vinf / length);
    };

    const double dq = S::DefaultRateStepFraction * base.Vinf / cbar;
    {
        S::FreestreamConditions p = base, m = base;
        p.q += dq;
        m.q -= dq;
        central(p, m, dq, [&](const auto& rp, const auto& rm, double h) {
            d.CL_q = (rp.CL - rm.CL) / h;
            d.Cm_q = (rp.Cm - rm.Cm) / h;
        });
        d.CL_q_nd = reduce(d.CL_q, cbar);
        d.Cm_q_nd = reduce(d.Cm_q, cbar);
    }
    const double dp = S::DefaultRateStepFraction * base.Vinf / span;
    {
        S::FreestreamConditions p = base, m = base;
        p.p += dp;
        m.p -= dp;
        central(p, m, dp, [&](const auto& rp, const auto& rm, double h) {
            d.Croll_p = (rp.Croll - rm.Croll) / h;
            d.Cn_p = (rp.Cn - rm.Cn) / h;
        });
        d.Croll_p_nd = reduce(d.Croll_p, span);
        d.Cn_p_nd = reduce(d.Cn_p, span);
    }
    const double dr = S::DefaultRateStepFraction * base.Vinf / span;
    {
        S::FreestreamConditions p = base, m = base;
        p.r += dr;
        m.r -= dr;
        central(p, m, dr, [&](const auto& rp, const auto& rm, double h) {
            d.Croll_r = (rp.Croll - rm.Croll) / h;
            d.Cn_r = (rp.Cn - rm.Cn) / h;
        });
        d.Croll_r_nd = reduce(d.Croll_r, span);
        d.Cn_r_nd = reduce(d.Cn_r, span);
    }
    return d;
}

const char* SeparationModeName(S::SeparationMode mode) {
    switch (mode) {
        case S::SeparationMode::Attached: return "attached";
        case S::SeparationMode::TurbulentSeparation: return "turbulent";
    }
    return "unknown";
}

const char* StateName(S::AttachmentLineState state) {
    switch (state) {
        case S::AttachmentLineState::Laminar: return "laminar";
        case S::AttachmentLineState::Contaminated: return "contaminated";
        case S::AttachmentLineState::Turbulent: return "turbulent";
    }
    return "unknown";
}

const char* TypeName(S::CriticalPointType type) {
    switch (type) {
        case S::CriticalPointType::AttachmentNode: return "attachment_node";
        case S::CriticalPointType::SeparationNode: return "separation_node";
        case S::CriticalPointType::Saddle: return "saddle";
        case S::CriticalPointType::AttachmentFocus: return "attachment_focus";
        case S::CriticalPointType::SeparationFocus: return "separation_focus";
        case S::CriticalPointType::Degenerate: return "degenerate";
    }
    return "unknown";
}

const char* ExitName(S::StreamlineExit exit) {
    switch (exit) {
        case S::StreamlineExit::StepLimit: return "step_limit";
        case S::StreamlineExit::LeftPatch: return "left_patch";
        case S::StreamlineExit::Stalled: return "stalled";
        case S::StreamlineExit::InvalidGrid: return "invalid_grid";
    }
    return "unknown";
}

void WriteVec3(std::ofstream& out, const Math::Vec3& v) {
    out << '[' << v.x << ',' << v.y << ',' << v.z << ']';
}

void WriteStation(std::ofstream& out, const S::AttachmentStation& station, double signedEta,
                  bool carryStrip, bool first) {
    if (!first) out << ",\n";
    out << R"(    {"eta":)" << signedEta << R"(,"found":)" << (station.Found ? "true" : "false")
        << R"(,"atKink":)" << (station.AtKink ? "true" : "false") << R"(,"carryStrip":)"
        << (carryStrip ? "true" : "false") << R"(,"sweepDeg":)" << Math::RadToDeg(station.SweepRad)
        << R"(,"effectiveSweepDeg":)" << Math::RadToDeg(station.EffectiveSweepRad)
        << R"(,"alphaNormalDeg":)" << station.AlphaNormalDeg << R"(,"stagnationOffset":)"
        << station.StagnationOffset << R"(,"strainRate":)" << station.StrainRate
        << R"(,"theta0":)" << station.MomentumThickness << R"(,"Rbar":)"
        << station.AttachmentLineReynolds << R"(,"state":")" << StateName(station.State)
        << R"(","point":)";
    WriteVec3(out, station.AttachmentPoint);
    out << '}';
}

void WriteAttachmentLine(std::ofstream& out, const S::AttachmentLine& line,
                         const std::vector<S::Panel>& panels, double halfSpan) {
    out << "[\n";
    for (std::size_t i = 0; i < line.Stations.size(); ++i) {
        const S::Vec3 mid = (panels[i].A + panels[i].B) * 0.5;
        WriteStation(out, line.Stations[i], mid.y / halfSpan, IsCarryThroughStrip(panels[i]), i == 0);
    }
    out << "\n  ]";
}

// The viewer's quad reconstruction (LatticeRenderer.cpp, and the same
// helper in LatticeFigureExport.cpp): the drawn quad is the panel's
// planform footprint around its bound quarter-chord segment.
void WriteWingQuad(std::ofstream& out, const S::Panel& panel, double gamma, bool first) {
    const double chord =
        (panel.SpanwiseWidth > 0.0) ? panel.PlanformArea / panel.SpanwiseWidth : 0.0;
    S::Vec3 chordDir = panel.ControlPoint - (panel.A + panel.B) * 0.5;
    const double len = chordDir.Norm();
    chordDir = (len > 1e-12) ? chordDir * (1.0 / len) : S::Vec3(1.0, 0.0, 0.0);

    const S::Vec3 corners[4] = {panel.A - chordDir * (0.25 * chord),
                                panel.B - chordDir * (0.25 * chord),
                                panel.B + chordDir * (0.75 * chord),
                                panel.A + chordDir * (0.75 * chord)};
    if (!first) out << ",\n";
    out << R"(    {"gamma":)" << gamma << R"(,"corners":[)";
    for (int i = 0; i < 4; ++i) {
        if (i) out << ',';
        WriteVec3(out, corners[i]);
    }
    out << "]}";
}

void WriteStreamline(std::ofstream& out, const S::SurfaceStreamline& line, const char* origin,
                     bool first) {
    if (!first) out << ",\n";
    const std::size_t count = line.Points.size();
    const std::size_t stride =
        (count > StreamlinePointBudget) ? (count + StreamlinePointBudget - 1) / StreamlinePointBudget
                                        : 1;
    out << R"(    {"origin":")" << origin << R"(","exit":")" << ExitName(line.Exit)
        << R"(","upstream":)" << (line.Upstream ? "true" : "false") << R"(,"length":)" << line.Length
        << R"(,"points":[)";
    bool firstPoint = true;
    for (std::size_t k = 0; k < count; k += stride) {
        const std::size_t index = (k + stride < count) ? k : count - 1; // always keep the endpoint
        const S::StreamlinePoint& point = line.Points[index];
        if (!firstPoint) out << ',';
        firstPoint = false;
        out << '[' << point.Point.x << ',' << point.Point.y << ',' << point.Point.z << ','
            << point.Cp << ',' << point.Speed << ']';
    }
    out << "]}";
}

} // namespace

int main(int argc, char** argv) {
    const std::string handoffPath =
        (argc > 1) ? argv[1] : "tests/Data/AeolionGeometryHandoff-1.8.0.json";
    const std::string outPath = (argc > 2) ? argv[2] : "attachment-sweep.json";
    // Fan thrust [N]. Zero (the default) is the unpowered airframe of the
    // second paper; positive turns the aft fan on and switches the sweep
    // from alpha x beta to alpha at the symmetry plane.
    const double thrust = (argc > 3) ? std::atof(argv[3]) : 0.0;
    // Airspeed [m/s]. A powered sweep MUST be able to vary this: transition
    // is low speed at high thrust, and running the fan at the cruise speed
    // instead understates its induction several-fold (the induced velocity
    // is comparable to flight speed only when the flight speed is low).
    const double flightSpeed = (argc > 4) ? std::atof(argv[4]) : FlightSpeed;

    // Optional alpha range: start, end, step. The default matrix steps 2
    // degrees, which is fine for a coefficient table and far too coarse to
    // locate a separation BOUNDARY -- a boundary is a threshold crossing,
    // and a shift smaller than the grid moves the crossing within an
    // interval without ever moving it across one. Resolving how much
    // incidence the fan buys needs the curve, not the crossing.
    std::vector<double> alphaList(std::begin(Alphas), std::end(Alphas));
    if (argc > 7) {
        const double start = std::atof(argv[5]), end = std::atof(argv[6]);
        const double step = std::atof(argv[7]);
        if (step > 0.0 && end >= start) {
            alphaList.clear();
            for (double a = start; a <= end + 1e-9; a += step) alphaList.push_back(a);
        }
    }

    Geometry::HandoffContract contract;
    try {
        contract = Geometry::LoadHandoff(handoffPath);
    } catch (const std::exception& error) {
        std::cerr << "cannot load " << handoffPath << ": " << error.what() << "\n";
        return 1;
    }
    contract.Mesh.ChordwisePanels = 1; // one Weissinger row per strip -- the attachment-line contract
    RefineNoseStations(contract.Body);

    PB::LatticeOptions carryOptions; // trim + carry-through: the physical solve
    carryOptions.BodyCircumferentialPanels = BodySectors;
    PB::LatticeBuilder carryBuilder(contract, carryOptions);
    PB::LatticeOptions cleanOptions = carryOptions;
    cleanOptions.CarryThroughLift = false; // true strip geometry for the attachment line
    PB::LatticeBuilder cleanBuilder(contract, cleanOptions);

    const auto wingCarry = carryBuilder.Build();
    const auto wingClean = cleanBuilder.Build();
    const auto body = carryBuilder.BuildBody();
    const auto duct = carryBuilder.BuildDuct();
    std::vector<Lattice::SourcePanel> sources = body;
    sources.insert(sources.end(), duct.begin(), duct.end());

    const double halfSpan = 0.5 * contract.Span;
    const auto stripsCarry = StripsFromPanels(wingCarry, halfSpan);
    const auto stripsClean = StripsFromPanels(wingClean, halfSpan);
    const auto& sections = contract.AirfoilSections;

    const double trail = TrailSpans * contract.Span;
    const auto preparedCarry = S::Prepare(S::PanelSystem{wingCarry, sources}, trail);
    const auto preparedClean = S::Prepare(S::PanelSystem{wingClean, sources}, trail);
    std::cout << "pivot ratios: carry " << preparedCarry.Factorization.MinPivotRatio << ", clean "
              << preparedClean.Factorization.MinPivotRatio << "\n";

    S::ReferenceGeometry ref;
    ref.Area = carryBuilder.GrossPlanformArea();
    ref.Span = contract.Span;
    ref.Chord = ref.Area / contract.Span;

    S::FreestreamConditions fc;
    fc.Vinf = flightSpeed;
    fc.rho = Rho;
    // Contract frame is x-forward / z-down, the solver x-aft / z-up: the
    // same 180-degree rotation about y that LatticeBuilder applies to the
    // placement anchor, the body stations and the hinge axes. Skipping it
    // puts the moment reference point an equal distance the WRONG side of
    // the origin, which shows up as a moment arm of twice its x -- roughly
    // a body length here, and a Cm_alpha an order of magnitude too large.
    if (contract.MomentReferencePointStated)
        fc.RefPoint = S::Vec3(-contract.MomentReferencePoint.x, contract.MomentReferencePoint.y,
                              -contract.MomentReferencePoint.z);

    std::ofstream out(outPath);
    if (!out) {
        std::cerr << "cannot open " << outPath << " for writing\n";
        return 1;
    }
    out.precision(9);

    out << "{\n\"meta\":{\"handoff\":\"" << contract.SchemaVersion << "\",\"Vinf\":" << flightSpeed
        << ",\"rho\":" << Rho << ",\"nu\":" << S::SeaLevelKinematicViscosity
        << ",\"span\":" << contract.Span << ",\"area\":" << ref.Area << ",\"chord\":" << ref.Chord
        << ",\"trimEta\":" << carryBuilder.TrimEta() << ",\"wingStrips\":" << wingCarry.size()
        << ",\"bodyPanels\":" << body.size() << ",\"ductPanels\":" << duct.size()
        << ",\"RbarContamination\":" << S::AttachmentLineContaminationReynolds
        << ",\"RbarTransition\":" << S::AttachmentLineTransitionReynolds
        << ",\"figureAlphaDeg\":" << FigureAlphaDeg << ",\"figureBetaDeg\":" << FigureBetaDeg
        << ",\"fanThrust\":" << thrust;
    {
        const S::ActuatorDisk probe = MakeFanDisk(contract, thrust, Rho, flightSpeed);
        out << ",\"fanDisk\":{\"active\":" << (probe.Valid() ? "true" : "false")
            << ",\"radius\":" << probe.Radius << ",\"hubRadius\":" << probe.HubRadius
            << ",\"area\":" << (probe.Valid() ? probe.Area() : 0.0)
            << ",\"x\":" << probe.Center.x << ",\"viAtCruise\":" << probe.InducedVelocity << "}";
    }
    out << ",\"figureAlphas\":[";
    for (std::size_t i = 0; i < std::size(FigureAlphas); ++i)
        out << (i ? "," : "") << FigureAlphas[i];
    out << "],\"figureBetas\":[";
    for (std::size_t i = 0; i < std::size(FigureBetas); ++i)
        out << (i ? "," : "") << FigureBetas[i];
    out << "]},\n\"conditions\":[\n";

    bool firstCondition = true;
    for (const double alphaDeg : alphaList) {
        for (const double betaDeg : (thrust > 0.0 ? std::span<const double>(PoweredBetas)
                                                  : std::span<const double>(Betas))) {
            fc.alphaDeg = alphaDeg;
            fc.betaDeg = betaDeg;

            // The disk is rebuilt per condition: momentum theory's induced
            // velocity depends on the axial speed through the disk, so a
            // fixed-thrust sweep is NOT a fixed-induction one.
            const S::ActuatorDisk disk = MakeFanDisk(contract, thrust, Rho, fc.Vinf);
            const auto fanField =
                disk.Valid() ? S::DiskInductionField(disk) : std::function<S::Vec3(const S::Vec3&)>{};

            const S::SolveResult carry = S::SolveWithSystem(preparedCarry, fc, ref, fanField);
            const S::SolveResult clean = S::SolveWithSystem(preparedClean, fc, ref, fanField);
            const S::FlowField fieldCarry =
                S::MakeFlowField(preparedCarry, fc, carry.gamma, carry.sigma, fanField);
            const S::FlowField fieldClean =
                S::MakeFlowField(preparedClean, fc, clean.gamma, clean.sigma, fanField);

            const S::AttachmentLine lineClean =
                S::ComputeAttachmentLine(fieldClean, wingClean, stripsClean, sections);
            const S::AttachmentLine lineCarry =
                S::ComputeAttachmentLine(fieldCarry, wingCarry, stripsCarry, sections);

            const S::SurfaceGrid grid = S::BuildSurfaceGrid(fieldCarry, PB::BodySurfaceName);
            const S::SurfaceFlowTopology topology = S::AnalyzeSurfaceFlow(grid);

            // Separation is surveyed on the CLEAN line for the same reason
            // the attachment line is: a carry-through strip's section
            // geometry is not its own.
            const S::SeparationSurvey separation = S::SurveySeparation(lineClean);
            const S::StabilityDerivatives derivatives =
                CoupledDerivatives(preparedCarry, fc, ref);

            if (disk.Valid() && alphaDeg == alphaList.front() && betaDeg == 0.0)
                std::cout << "fan: T=" << thrust << " N, disk r=" << disk.Radius << " hub="
                          << disk.HubRadius << " at x=" << disk.Center.x
                          << ", vi=" << disk.InducedVelocity << " m/s\n";
            std::cout << "alpha=" << alphaDeg << " beta=" << betaDeg << ": CL=" << carry.CL
                      << " Cm=" << carry.Cm << " CLa=" << derivatives.CL_alpha
                      << " Cnb=" << derivatives.Cn_beta << " sep=" << separation.SeparatedStations
                      << "/" << separation.ResolvedStations;
            if (separation.AnySeparated())
                std::cout << " fwd x/c=" << separation.ForwardmostSeparationPsi;
            std::cout << (topology.AttachesUpstream ? "  body-attaches-upstream" : "") << "\n";

            if (!firstCondition) out << ",\n";
            firstCondition = false;
            // Far-field induced drag alongside the near-field one. On a
            // coupled configuration the near-field number collects a
            // spurious thrust from the discretized closed bodies; the
            // Trefftz integral cannot, because a closed body sheds no wake.
            const S::TrefftzResult trefftz = S::TrefftzInducedDrag(
                wingCarry, carry.gamma, ref, Rho, fc.Vinf, carry.CL);

            out << R"( {"alphaDeg":)" << alphaDeg << R"(,"betaDeg":)" << betaDeg << R"(,"CL":)"
                << carry.CL << R"(,"CDi":)" << carry.CDi << R"(,"CDiTrefftz":)" << trefftz.CDi
                << R"(,"spanEfficiency":)" << trefftz.SpanEfficiency
                << R"(,"liftWing":)" << (carry.LiftBySurface.count("wing") ? carry.LiftBySurface.at("wing") : 0.0)
                << R"(,"liftTotal":)" << carry.L << R"(,"CY":)" << carry.CY << R"(,"Cm":)"
                << carry.Cm << R"(,"Croll":)" << carry.Croll << R"(,"Cn":)" << carry.Cn
                << R"(,"CLClean":)" << clean.CL << ",\n  \"derivatives\":{"
                << R"("CL_alpha":)" << derivatives.CL_alpha << R"(,"CDi_alpha":)"
                << derivatives.CDi_alpha << R"(,"Cm_alpha":)" << derivatives.Cm_alpha
                << R"(,"CY_beta":)" << derivatives.CY_beta << R"(,"Croll_beta":)"
                << derivatives.Croll_beta << R"(,"Cn_beta":)" << derivatives.Cn_beta
                << R"(,"CL_q_nd":)" << derivatives.CL_q_nd << R"(,"Cm_q_nd":)"
                << derivatives.Cm_q_nd << R"(,"Croll_p_nd":)" << derivatives.Croll_p_nd
                << R"(,"Cn_p_nd":)" << derivatives.Cn_p_nd << R"(,"Croll_r_nd":)"
                << derivatives.Croll_r_nd << R"(,"Cn_r_nd":)" << derivatives.Cn_r_nd << "},\n"
                << R"(  "separation":{"resolved":)" << separation.ResolvedStations
                << R"(,"separated":)" << separation.SeparatedStations << R"(,"fraction":)"
                << separation.SeparatedFraction() << R"(,"forwardmostPsi":)"
                << separation.ForwardmostSeparationPsi << R"(,"forwardmostEta":)"
                << separation.ForwardmostSeparationEta << R"(,"stations":[)";
            {
                // StripSection::Eta is a semi-span FRACTION and therefore
                // unsigned, which would collapse the two wings onto each
                // other -- fatal for reading a sideslip asymmetry. The
                // signed station comes from the panel the strip was built
                // from, the same way the attachment line reports it.
                bool firstStation = true;
                for (const S::StationSeparation& entry : separation.Stations) {
                    if (!entry.Resolved) continue;
                    const std::size_t k = static_cast<std::size_t>(entry.Strip);
                    const double signedEta =
                        (k < wingClean.size())
                            ? ((wingClean[k].A + wingClean[k].B) * 0.5).y / halfSpan
                            : entry.Eta;
                    if (!firstStation) out << ',';
                    firstStation = false;
                    out << R"({"eta":)" << signedEta << R"(,"mode":")"
                        << SeparationModeName(entry.Upper.Mode) << R"(","psi":)"
                        << entry.Upper.SeparationPsi << R"(,"transitionPsi":)"
                        << (entry.Upper.Transitioned ? entry.Upper.TransitionPsi : -1.0)
                        << R"(,"bubble":)" << (entry.Upper.BubbleFormed ? "true" : "false") << '}';
                }
            }
            out << "]},\n  \"attachmentLine\":";
            WriteAttachmentLine(out, lineClean, wingClean, halfSpan);
            out << ",\n  \"attachmentLineCarry\":";
            WriteAttachmentLine(out, lineCarry, wingCarry, halfSpan);

            out << ",\n  \"body\":{\"gridValid\":" << (grid.Valid() ? "true" : "false")
                << ",\"attachesUpstream\":" << (topology.AttachesUpstream ? "true" : "false")
                << ",\"noseAttachment\":";
            WriteVec3(out, topology.NoseAttachment);
            out << ",\"criticalPoints\":[\n";
            for (std::size_t i = 0; i < topology.CriticalPoints.size(); ++i) {
                const S::CriticalPoint& point = topology.CriticalPoints[i];
                if (i) out << ",\n";
                out << R"(    {"type":")" << TypeName(point.Type) << R"(","point":)";
                WriteVec3(out, point.Point);
                out << R"(,"station":)" << point.Station << R"(,"sector":)" << point.Sector
                    << R"(,"trace":)" << point.Trace << R"(,"determinant":)" << point.Determinant
                    << '}';
            }
            out << "\n  ]}";

            // The portrait condition additionally carries the fuselage
            // surface mesh with Cp, the wing lattice, and the streamlines.
            const bool isFigure = (alphaDeg == FigureAlphaDeg) && (betaDeg == FigureBetaDeg);
            if (isFigure && grid.Valid()) {
                out << ",\n  \"surface\":{\"stations\":" << grid.Stations
                    << ",\"sectors\":" << grid.Sectors << ",\"points\":[";
                for (std::size_t k = 0; k < grid.Samples.size(); ++k) {
                    const S::SurfaceSample& sample = grid.Samples[k];
                    if (k) out << ',';
                    out << '[' << sample.Point.x << ',' << sample.Point.y << ',' << sample.Point.z
                        << ',' << sample.Cp << ']';
                }
                out << "]},\n  \"wingPanels\":[\n";
                // The clean solve's wing: its quads stop at the body surface
                // (no carry-through segment crossing the fuselage) and it is
                // the solve the attachment line is computed on.
                for (std::size_t i = 0; i < wingClean.size(); ++i)
                    WriteWingQuad(out, wingClean[i],
                                  (i < clean.gamma.size()) ? clean.gamma[i] : 0.0, i == 0);
                out << "\n  ],\n  \"streamlines\":[\n";

                bool firstLine = true;
                // The distinguished lines of every critical point: the
                // attachment node's principal lines ARE the body's
                // attachment lines, a separation node's are the separation
                // lines (traced upstream by TraceFromCriticalPoint itself).
                for (const S::CriticalPoint& point : topology.CriticalPoints) {
                    for (const S::SurfaceStreamline& line : S::TraceFromCriticalPoint(grid, point)) {
                        WriteStreamline(out, line, TypeName(point.Type), firstLine);
                        firstLine = false;
                    }
                }
                // A ring of ordinary streamlines seeded just aft of the nose,
                // one per sector: the skin-flow portrait itself.
                for (int j = 0; j < grid.Sectors; ++j) {
                    const S::SurfaceStreamline line = S::TraceSurfaceStreamline(grid, 1.5, j);
                    WriteStreamline(out, line, "ring", firstLine);
                    firstLine = false;
                }
                out << "\n  ]";
            }
            out << '}';
        }
    }
    out << "\n]\n}\n";

    std::cout << "wrote " << outPath << "\n";
    return 0;
}
