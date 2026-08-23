// Solver/SolveResult.h
//
// Full result of a VLM solve: circulation distribution, per-station
// breakdown, integrated force/moment coefficients and their dimensional
// counterparts, and per-surface bookkeeping.
#pragma once
#include <vector>
#include <string>
#include <map>
#include "Aeolion/Solver/StationResult.h"

namespace Aeolion::Solver {

/** Full result of a VLM solve. */
struct SolveResult {
    std::vector<double> gamma;
    /**
     * Converged SOURCE strengths, aligned with the system's source panels
     * (empty for a wing-only solve). Kept for the same reason gamma is:
     * together they are the full state of the solved field, and anything
     * that wants to re-evaluate the velocity after the fact -- a surface
     * streamline, a stagnation point, a wake survey -- needs both to
     * rebuild a Solver::FlowField.
     */
    std::vector<double> sigma;
    std::vector<StationResult> Stations;
    double CL = 0.0;
    /**
     * INDUCED drag only -- VLM is a potential-flow method and cannot predict
     * viscous/profile drag. Total CD = CDi + your own CD0 estimate.
     *
     * NEAR-FIELD, and on a COUPLED CONFIGURATION that is not the number you
     * want. This is the streamwise component of the Kutta-Joukowski forces,
     * so it is a small difference of much larger lift-dominated quantities,
     * and adding a source-panelled body to a lifting surface degrades it
     * visibly: it goes NEGATIVE at zero lift, where induced drag must vanish,
     * and an alpha sweep of the airframe in tests/Data fits to an Oswald
     * efficiency of 1.53, which is impossible for a planar wing. CL and the
     * moments are unaffected -- the same absolute error is negligible beside
     * lift and comparable with induced drag. For a bare wing near and far
     * field agree to 0.4% and this field is fine.
     *
     * Use Solver::TrefftzInducedDrag (Solver/TrefftzPlane.h) when the induced
     * drag itself matters on a configuration with a body. The far-field
     * integral cannot inherit this error by construction: a closed body sheds
     * no trailing vorticity, so it puts nothing through a plane at downstream
     * infinity and is never evaluated. That header carries the full argument,
     * including why d'Alembert is NOT the explanation (closed bodies here
     * carry zero net force in uniform flow to machine precision).
     *
     * It is deliberately not computed here, and deliberately does not replace
     * this field. Three reasons, each sufficient on its own. Propeller thrust
     * is -Di (PanelBuilder.h) and CDi = Di/(q S), so changing one without the
     * other breaks the identity and changing both breaks the rotor. A
     * rotating-frame rotor sheds a helical wake, which is not what a plane at
     * downstream infinity models. And the coupled solver calls Solve on the
     * order of a thousand times per condition without ever reading CDi, so
     * an O(strips^2) wake integral on every call would be paid entirely in
     * sweeps that discard it.
     */
    double CDi = 0.0;
    double CY = 0.0;
    double Croll = 0.0; ///< Rolling moment coefficient about RefPoint (positive: right wingtip up).
    double Cm = 0.0;    ///< Pitching moment coefficient about RefPoint (positive: nose up).
    double Cn = 0.0;    ///< Yawing moment coefficient about RefPoint (positive: nose right, about +z).
    double L = 0.0, Di = 0.0, Y = 0.0;
    double Mx = 0.0, My = 0.0, Mz = 0.0; ///< Moments about RefPoint [N*m], same axis convention as Croll/Cm/Cn.
    double ReferenceArea = 0.0, ReferenceChord = 0.0, ReferenceSpan = 0.0;
    std::map<std::string, double> LiftBySurface;  ///< [N], per source surface tag.
    std::map<std::string, double> DragBySurface;  ///< [N]
    std::map<std::string, double> AreaBySurface;  ///< [m^2]
};

} // namespace Aeolion::Solver
