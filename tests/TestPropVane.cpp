// TestPropVane.cpp -- integration test tying Bemt.h's slipstream field
// into Vlm.h's externalField hook: a deflected vane sitting behind a
// running propeller, at essentially zero airspeed (hover), should produce
// real side force / yaw moment purely from the propwash -- while the same
// vane with no propwash present produces ~nothing. This is the literal
// physical mechanism that gives a thrust-vectoring-vane aircraft hover
// control authority, so it's worth a standing regression check.
#include "Aeolion/Bemt.h"
#include "Aeolion/Vlm.h"
#include <iostream>
#include <cmath>

using namespace Aeolion;
using namespace Aeolion::Vlm;

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } } while (0)

int main() {
    Bemt::PropGeometry geom;
    geom.NBlades = 2; geom.Radius = 0.15; geom.HubRadius = 0.02;
    geom.RotationSign = 1;
    int N = 12;
    for (int i = 0; i < N; ++i) {
        double r = geom.HubRadius + (geom.Radius - geom.HubRadius) * (i + 0.5) / (N + 1);
        double eta = r / geom.Radius;
        geom.Stations.push_back({r, 0.025 * (1 - 0.3 * eta), 30.0 - 20.0 * eta});
    }
    Bemt::Polar polar;
    double rho = 1.225, rpm = 6000.0;

    WingParams vaneParams;
    vaneParams.Span = 0.10; vaneParams.RootChord = 0.04; vaneParams.TipChord = 0.04;
    vaneParams.NPanelsSemiSpan = 8; vaneParams.CosineSpacing = true;
    std::vector<Panel> vanePanelsFlat = BuildWing(vaneParams);

    std::vector<Panel> vanePanels;
    for (auto p : vanePanelsFlat) {
        auto rot = [](Vec3 v) { return Vec3(v.x, -v.z, v.y); }; // 90deg about x: span along z instead of y
        p.A = rot(p.A); p.B = rot(p.B); p.ControlPoint = rot(p.ControlPoint); p.Normal = rot(p.Normal);
        p.TrailDirA = Vec3(1, 0, 0); p.TrailDirB = Vec3(1, 0, 0);
        p.Surface = "vane";
        vanePanels.push_back(p);
    }
    double defRad = 15.0 * Pi / 180.0;
    for (auto& p : vanePanels) {
        Vec3 n = p.Normal;
        double c = std::cos(defRad), s = std::sin(defRad);
        p.Normal = Vec3(n.x * c - n.y * s, n.x * s + n.y * c, n.z).Normalized();
        p.A.x += 0.08; p.B.x += 0.08; p.ControlPoint.x += 0.08;
    }

    ReferenceGeometry ref;
    ref.Area = 0.0; for (auto& p : vanePanels) ref.Area += p.Area;
    ref.Span = vaneParams.Span; ref.Chord = ref.Area / ref.Span;

    FreestreamConditions fc;
    fc.Vinf = 0.001; fc.alphaDeg = 0; fc.rho = rho; // near-zero airspeed: pure hover

    SolveResult noProp = Solve(vanePanels, fc, ref, 5.0);
    std::cout << "no propwash: Y=" << noProp.Y << " N\n";
    CHECK(std::fabs(noProp.Y) < 1e-3, "vane with no propwash and no airspeed should produce ~zero force");

    auto bemtRes = Bemt::Solve(geom, polar, rpm, 0.0, rho);
    CHECK(bemtRes.Converged, "propeller solve should converge");
    CHECK(bemtRes.Thrust > 0, "propeller should produce positive thrust in hover");

    Bemt::SlipstreamField field;
    field.BemtResult = bemtRes;
    field.HubCenter = Vec3(0, 0, 0);
    field.AxisDir = Vec3(1, 0, 0);
    field.DevelopmentLength = 0.15;

    SolveResult withProp = Solve(vanePanels, fc, ref, 5.0, field);
    std::cout << "with propwash: Y=" << withProp.Y << " N,  Mz=" << withProp.Mz << " N*m\n";
    CHECK(std::fabs(withProp.Y) > 0.01, "deflected vane in the propwash MUST produce meaningful side force at zero airspeed -- this is the whole point of a thrust-vectoring vane");
    CHECK(std::fabs(withProp.Mz) > 0.001, "deflected vane in the propwash should produce a real yaw moment");

    if (failures == 0) { std::cout << "PASS: TestPropVane\n"; return 0; }
    std::cerr << failures << " check(s) failed in TestPropVane\n";
    return 1;
}
