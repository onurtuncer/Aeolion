// TestDiskInduction.cpp -- validates the actuator disk's vortex-cylinder
// field against closed-form answers, not against itself.
//
//   - Momentum theory. The three values a loaded disk must produce on its
//     own axis are v_i at the disk plane, 2 v_i far downstream and zero far
//     upstream. The discretized sheet has to reproduce them without being
//     told: they are consequences of the vortex system, not inputs to it.
//
//   - The discretized field against the exact axis solution, at a spread of
//     stations INCLUDING upstream ones -- which is the half of the field
//     this module exists for and the half SlipstreamField does not carry.
//
//   - Convergence. Refining the rings and the polygon must move the answer
//     toward the closed form, not merely change it.
//
//   - The annulus. A uniformly loaded annular disk sheds at both edges, so
//     its axis field is the difference of two cylinders. That difference
//     vanishes AT the disk plane and asymptotically, and in between it does
//     not -- downstream inside the bore it goes NEGATIVE, which is the
//     centerbody wake and a real feature of annular jets rather than a
//     defect. Checking the discretized annulus against its own closed form
//     is what pins the inner sheet's sign.
//
//   - Direction. The induced velocity must point DOWNSTREAM inside the
//     cylinder for a thrusting disk. Getting the ring circulation sense
//     backwards would turn the upstream acceleration this whole paper rests
//     on into a deceleration, and every magnitude check above would still
//     pass.
#include "Aeolion/Solver/DiskInduction.h"

#include <cmath>
#include <iostream>
#include <string>

namespace S = Aeolion::Solver;
using Aeolion::Math::Vec3;

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } } while (0)

namespace {

constexpr double Radius = 0.1045; // the JETTAIL duct bore
constexpr double Vi = 12.0;       // a representative disk-plane induction [m/s]

S::ActuatorDisk MakeDisk(double hubRadius = 0.0) {
    S::ActuatorDisk disk;
    disk.Center = Vec3(0, 0, 0);
    disk.Axis = Vec3(1, 0, 0); // downstream = +x, the solver's aft direction
    disk.Radius = Radius;
    disk.HubRadius = hubRadius;
    disk.InducedVelocity = Vi;
    return disk;
}

void TestMomentumTheoryOnTheAxis() {
    const S::ActuatorDisk disk = MakeDisk();

    CHECK(std::fabs(S::DiskAxisInducedVelocity(disk, 0.0) - Vi) < 1e-12,
          "the closed form must give v_i at the disk plane, got " +
              std::to_string(S::DiskAxisInducedVelocity(disk, 0.0)));

    const double far = S::DiskAxisInducedVelocity(disk, 400.0 * Radius);
    CHECK(std::fabs(far - 2.0 * Vi) < 0.01 * Vi,
          "the fully developed wake must reach 2 v_i, got " + std::to_string(far));

    const double ahead = S::DiskAxisInducedVelocity(disk, -400.0 * Radius);
    CHECK(std::fabs(ahead) < 0.01 * Vi,
          "the induction must vanish far upstream, got " + std::to_string(ahead));

    // One radius upstream -- roughly where a wing sits -- the exact value is
    // v_i (1 - 1/sqrt(2)). Worth pinning explicitly: it is the number the
    // whole aft-fan argument turns on.
    const double atOneRadius = S::DiskAxisInducedVelocity(disk, -Radius);
    const double expected = Vi * (1.0 - 1.0 / std::sqrt(2.0));
    CHECK(std::fabs(atOneRadius - expected) < 1e-12,
          "one radius upstream the induction must be v_i(1 - 1/sqrt2) = " +
              std::to_string(expected) + ", got " + std::to_string(atOneRadius));
    CHECK(atOneRadius > 0.25 * Vi,
          "and it must be a substantial fraction of v_i, not a rounding error");
}

void TestDiscretizedSheetMatchesTheAxisSolution() {
    const S::ActuatorDisk disk = MakeDisk();
    const S::VortexCylinderMesh mesh = S::BuildVortexCylinder(disk);
    CHECK(!mesh.Empty(), "a loaded disk must build a sheet");

    std::cout << "disk axis field (" << mesh.Count() << " filaments):\n"
              << "      s/R      exact     discrete\n";
    double worst = 0.0;
    for (const double sOverR : {-4.0, -2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 4.0, 10.0}) {
        const double s = sOverR * Radius;
        const double exact = S::DiskAxisInducedVelocity(disk, s);
        const Vec3 v = S::VortexCylinderVelocity(mesh, Vec3(s, 0.0, 0.0));
        std::printf("  %8.2f  %9.4f  %9.4f\n", sOverR, exact, v.x);

        // Off-axis components must vanish by symmetry on the axis itself.
        CHECK(std::fabs(v.y) < 1e-6 * Vi && std::fabs(v.z) < 1e-6 * Vi,
              "the on-axis field must be purely axial at s/R=" + std::to_string(sOverR));
        worst = std::max(worst, std::fabs(v.x - exact));
    }
    std::cout << "  worst axial error " << worst << " m/s (" << 100.0 * worst / Vi << "% of v_i)\n";
    CHECK(worst < 0.02 * Vi,
          "the discretized sheet must track the closed-form axis solution to 2% of v_i, worst " +
              std::to_string(worst));
}

void TestRefinementConvergesToTheClosedForm() {
    const S::ActuatorDisk disk = MakeDisk();
    const double station = -Radius; // one radius upstream
    const double exact = S::DiskAxisInducedVelocity(disk, station);

    double previous = 1e30;
    for (const int rings : {4, 8, 16}) {
        S::VortexCylinderOptions options;
        options.RingsPerRadius = rings;
        options.SegmentsPerRing = 6 * rings;
        const S::VortexCylinderMesh mesh = S::BuildVortexCylinder(disk, options);
        const double error =
            std::fabs(S::VortexCylinderVelocity(mesh, Vec3(station, 0, 0)).x - exact);
        std::cout << "  rings/R=" << rings << " -> error " << error << "\n";
        CHECK(error < previous,
              "refinement must reduce the error, not merely change it: " + std::to_string(error) +
                  " after " + std::to_string(previous));
        previous = error;
    }
}

void TestAnnulusAxisFieldAndItsWake() {
    const S::ActuatorDisk annulus = MakeDisk(0.0406); // the JETTAIL tail boom
    CHECK(annulus.Valid(), "an annular disk must be valid");
    CHECK(annulus.Area() < MakeDisk().Area(), "the bore must remove area");

    // The two sheets cancel exactly where their own terms vanish -- at the
    // disk plane -- and asymptotically far from it.
    CHECK(std::fabs(S::DiskAxisInducedVelocity(annulus, 0.0)) < 1e-12,
          "at the disk plane the two sheets must cancel on the axis");
    CHECK(std::fabs(S::DiskAxisInducedVelocity(annulus, 500.0 * Radius)) < 0.01 * Vi,
          "far downstream the annulus's axis induction must decay");
    CHECK(std::fabs(S::DiskAxisInducedVelocity(annulus, -500.0 * Radius)) < 0.01 * Vi,
          "far upstream likewise");

    // Ahead of the disk the axis flow is accelerated, but far more weakly
    // than the jet itself -- the hub shields the axis.
    const double ahead = S::DiskAxisInducedVelocity(annulus, -2.0 * Radius);
    CHECK(ahead > 0.0, "upstream on the axis the annulus must still accelerate the flow");
    CHECK(ahead < S::DiskAxisInducedVelocity(MakeDisk(), -2.0 * Radius),
          "but less than the equivalent full disk, which has no bore to shield the axis");

    // Downstream INSIDE the bore the inner sheet dominates and the axial
    // induction reverses: the centerbody wake.
    const double behind = S::DiskAxisInducedVelocity(annulus, Radius);
    std::cout << "  annulus axis: " << ahead << " m/s at s=-2R, " << behind << " m/s at s=+R\n";
    CHECK(behind < 0.0,
          "inside the bore behind an annular disk the flow must be retarded (the centerbody "
          "wake), got " + std::to_string(behind));

    // The discretized annulus must reproduce that whole structure, which is
    // the real check on the inner sheet's sign and strength.
    const S::VortexCylinderMesh mesh = S::BuildVortexCylinder(annulus);
    double worst = 0.0;
    for (const double sOverR : {-4.0, -2.0, -0.5, 0.0, 1.0, 4.0}) {
        const double exact = S::DiskAxisInducedVelocity(annulus, sOverR * Radius);
        const Vec3 v = S::VortexCylinderVelocity(mesh, Vec3(sOverR * Radius, 0, 0));
        worst = std::max(worst, std::fabs(v.x - exact));
    }
    std::cout << "  annulus worst axial error " << worst << " m/s\n";
    CHECK(worst < 0.03 * Vi,
          "the discretized annulus must track its own closed form, worst " + std::to_string(worst));

    // And between the radii, where the jet actually is, it must accelerate.
    const double rMid = 0.5 * (annulus.HubRadius + annulus.Radius);
    const Vec3 inJet = S::VortexCylinderVelocity(mesh, Vec3(0.0, rMid, 0.0));
    std::cout << "  annulus at mid-bore radius, disk plane: u_x=" << inJet.x << "\n";
    CHECK(inJet.x > 0.3 * Vi,
          "the annular jet itself must carry a substantial axial induction, got " +
              std::to_string(inJet.x));
}

void TestInductionPointsDownstream() {
    // The sign test the magnitude checks cannot catch. A thrusting disk
    // accelerates the flow ALONG its axis, both in the jet and ahead of it.
    const S::ActuatorDisk disk = MakeDisk();
    const S::VortexCylinderMesh mesh = S::BuildVortexCylinder(disk);

    for (const double sOverR : {-2.0, -1.0, 0.0, 2.0}) {
        const Vec3 v = S::VortexCylinderVelocity(mesh, Vec3(sOverR * Radius, 0, 0));
        CHECK(v.x > 0.0, "the induction must point downstream (+x) at s/R=" +
                             std::to_string(sOverR) + ", got " + std::to_string(v.x));
    }

    // And it must reverse with the axis, not with the coordinate: a disk
    // facing the other way pushes the other way.
    S::ActuatorDisk reversed = disk;
    reversed.Axis = Vec3(-1, 0, 0);
    const S::VortexCylinderMesh flipped = S::BuildVortexCylinder(reversed);
    const Vec3 v = S::VortexCylinderVelocity(flipped, Vec3(Radius, 0, 0)); // now UPSTREAM of it
    CHECK(v.x < 0.0, "reversing the disk axis must reverse the induced flow, got " +
                         std::to_string(v.x));
}

void TestMomentumTheoryInducedVelocity() {
    const double rho = 1.225, area = MakeDisk().Area();

    // Hover: v_i = sqrt(T / (2 rho A)).
    const double hover = S::InducedVelocityFromThrust(25.0, rho, area, 0.0);
    const double expected = std::sqrt(25.0 / (2.0 * rho * area));
    CHECK(std::fabs(hover - expected) < 1e-9,
          "at zero airspeed v_i must be the hover value " + std::to_string(expected) + ", got " +
              std::to_string(hover));

    // Forward flight relieves the disk: same thrust needs less induction.
    const double cruise = S::InducedVelocityFromThrust(25.0, rho, area, 25.0);
    CHECK(cruise < hover, "forward speed must reduce the induced velocity at fixed thrust");
    CHECK(cruise > 0.0, "but it must stay positive for a thrusting disk");

    // And the momentum balance it was inverted from must hold.
    const double thrust = 2.0 * rho * area * cruise * (25.0 + cruise);
    CHECK(std::fabs(thrust - 25.0) < 1e-9,
          "v_i must satisfy T = 2 rho A v_i (V + v_i), recovered " + std::to_string(thrust));

    CHECK(S::InducedVelocityFromThrust(0.0, rho, area, 10.0) == 0.0,
          "an unloaded disk induces nothing");
    CHECK(S::InducedVelocityFromThrust(-5.0, rho, area, 10.0) == 0.0,
          "a windmilling disk is a different momentum balance and must be declined");
}

void TestDegenerateDisksAreDeclined() {
    S::ActuatorDisk empty;
    CHECK(!empty.Valid(), "a zero-radius disk is not valid");
    CHECK(S::BuildVortexCylinder(empty).Empty(), "and must build no sheet");

    S::ActuatorDisk unloaded = MakeDisk();
    unloaded.InducedVelocity = 0.0;
    CHECK(S::BuildVortexCylinder(unloaded).Empty(), "an unloaded disk must build no sheet");

    S::ActuatorDisk swallowed = MakeDisk();
    swallowed.HubRadius = swallowed.Radius; // bore equals the tip
    CHECK(!swallowed.Valid(), "a hub at the tip radius leaves no annulus");
}

} // namespace

int main() {
    TestMomentumTheoryOnTheAxis();
    TestDiscretizedSheetMatchesTheAxisSolution();
    TestRefinementConvergesToTheClosedForm();
    TestAnnulusAxisFieldAndItsWake();
    TestInductionPointsDownstream();
    TestMomentumTheoryInducedVelocity();
    TestDegenerateDisksAreDeclined();

    if (failures == 0) std::cout << "PASS: TestDiskInduction\n";
    return failures == 0 ? 0 : 1;
}
