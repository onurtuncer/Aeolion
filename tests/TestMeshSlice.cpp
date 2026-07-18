// TestMeshSlice.cpp -- validates the mesh-slicing pipeline
// (MeshSlice.h: slice a triangle-mesh wing, extract camber line, build
// VLM panels) by comparing against the parametric VLM::BuildWing() path
// for the SAME planform, which was independently validated in
// TestVLMCore.cpp. Agreement here confirms the slicing/camber-extraction
// geometry is correct, not just that the solver itself is correct.
#include "Aeolion/VLM.h"
#include "Aeolion/MeshSlice.h"
#include "fixtures/SyntheticWing.h"
#include <iostream>
#include <cmath>

using namespace Aeolion;
using namespace Aeolion::VLM;

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } } while (0)

static void RunCase(double sweep, double dihedral, double twist, double tolCL, double tolCDi) {
    double span = 4.0, rootC = 0.5, tipC = 0.5;
    FreestreamConditions fc; fc.Vinf = 30.0; fc.alphaDeg = 4.0; fc.rho = 1.225;

    WingParams wp;
    wp.Span = span; wp.RootChord = rootC; wp.TipChord = tipC;
    wp.SweepQuarterChordDeg = sweep; wp.DihedralDeg = dihedral; wp.TwistTipDeg = twist;
    wp.NPanelsSemiSpan = 30; wp.CosineSpacing = true;
    SolveResult refRes = Solve(wp, fc);

    MeshIO::PartMesh mesh = BuildSyntheticWingMesh(span, rootC, tipC, sweep, dihedral, twist, 80, 24, 0.03);

    MeshSlice::SurfaceMeshParams sp;
    sp.SurfaceName = "wing";
    sp.SpanAxis = 1; sp.ChordAxis = 0; sp.ThicknessAxis = 2;
    sp.NPanelsSpan = 60; sp.CosineSpacing = true;
    sp.ReferenceUp = Vec3(0, 0, 1);

    std::vector<Panel> meshPanels = MeshSlice::MeshToPanels(mesh, sp);
    double S = 0.0; for (auto& p : meshPanels) S += p.Area;
    ReferenceGeometry ref; ref.Area = S; ref.Span = span; ref.Chord = S / span;
    SolveResult meshRes = Solve(meshPanels, fc, ref, 50.0 * span);

    double relErrCL = std::fabs(meshRes.CL - refRes.CL) / refRes.CL;
    double relErrCDi = std::fabs(meshRes.CDi - refRes.CDi) / refRes.CDi;
    std::cout << "sweep=" << sweep << " dihedral=" << dihedral << " twist=" << twist
              << "  CL: param=" << refRes.CL << " mesh=" << meshRes.CL << " (" << relErrCL * 100 << "%)"
              << "  CDi: param=" << refRes.CDi << " mesh=" << meshRes.CDi << " (" << relErrCDi * 100 << "%)\n";
    CHECK(relErrCL < tolCL, "mesh-derived CL diverged too far from the parametric result for sweep=" << sweep << " dihedral=" << dihedral << " twist=" << twist);
    CHECK(relErrCDi < tolCDi, "mesh-derived CDi diverged too far from the parametric result for sweep=" << sweep << " dihedral=" << dihedral << " twist=" << twist);
}

int main() {
    RunCase(0.0, 0.0, 0.0, 0.01, 0.05);       // plain rectangular wing
    RunCase(15.0, 5.0, -4.0, 0.02, 0.05);     // combined sweep+dihedral+twist

    if (failures == 0) { std::cout << "PASS: TestMeshSlice\n"; return 0; }
    std::cerr << failures << " check(s) failed in TestMeshSlice\n";
    return 1;
}
