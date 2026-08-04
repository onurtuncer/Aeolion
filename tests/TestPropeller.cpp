// TestPropeller.cpp -- the handoff-to-metric-propeller conversion
// (Geometry::ToPropeller).
//
// The conversion is unit arithmetic, which is exactly the sort of code that
// looks obviously right and silently isn't: the contract states radii as
// fractions of the disk radius and consumers want metres. A factor-of-R
// error here would shift every radius and still produce a plausible-looking
// propeller.
#include "Aeolion/Geometry/HandoffContract.h"

#include <cmath>
#include <iostream>
#include <string>

using Aeolion::Geometry::ContractError;

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } } while (0)

namespace {

Aeolion::Geometry::HandoffContract Fixture() {
    return Aeolion::Geometry::LoadHandoff(std::string(AEOLION_TEST_DATA_DIR) +
                                          "/AeolionGeometryHandoff-1.5.0.json");
}

void TestConversionConvertsUnits() {
    const auto contract = Fixture();
    const auto& spec = contract.Propulsion;
    const auto prop = Aeolion::Geometry::ToPropeller(spec);

    CHECK(prop.BladeCount == 3, "blade count should come from the contract");
    CHECK(prop.Stations.size() == spec.BladeStations.size(), "every blade station should survive the conversion");
    CHECK(std::fabs(prop.Radius - spec.DiskRadius) < 1e-15, "tip radius is the disk radius");

    // The contract defines the blade only outboard of its first station, so
    // that station IS the hub cutout.
    CHECK(std::fabs(prop.HubRadius - spec.BladeStations.front().RadiusFraction * spec.DiskRadius) < 1e-15,
          "hub radius is the innermost station, in metres");

    // Fractions became metres, and chord/twist passed through untouched.
    for (std::size_t i = 0; i < prop.Stations.size(); ++i) {
        const double expected = spec.BladeStations[i].RadiusFraction * spec.DiskRadius;
        CHECK(std::fabs(prop.Stations[i].r - expected) < 1e-15,
              "station " + std::to_string(i) + " radius should be r/R times R");
        CHECK(prop.Stations[i].Chord == spec.BladeStations[i].Chord, "chord is already metric");
        CHECK(prop.Stations[i].TwistDeg == spec.BladeStations[i].TwistDeg, "twist passes through");
        if (i > 0)
            CHECK(prop.Stations[i].r > prop.Stations[i - 1].r, "radii must increase outboard");
    }
    CHECK(prop.Stations.back().r <= prop.Radius + 1e-15, "no station may sit outboard of the tip");
}

void TestConversionRefusesWhatItCannotKnow() {
    auto contract = Fixture();

    // Blade count arrived at schema 1.4.0. An older contract states none,
    // and guessing would change everything a consumer derives from the
    // propeller while looking fine.
    contract.Propulsion.BladeCount = 0;
    bool rejected = false;
    try {
        (void)Aeolion::Geometry::ToPropeller(contract.Propulsion);
    } catch (const ContractError&) { rejected = true; }
    CHECK(rejected, "a contract with no blade count should be refused, not defaulted");
    CHECK(Aeolion::Geometry::ToPropeller(contract.Propulsion, 3).BladeCount == 3,
          "the caller should be able to supply the missing blade count");

    // A blade needs at least two stations to have a shape at all.
    auto empty = Fixture().Propulsion;
    empty.BladeStations.clear();
    rejected = false;
    try {
        (void)Aeolion::Geometry::ToPropeller(empty);
    } catch (const ContractError&) { rejected = true; }
    CHECK(rejected, "a propulsion block with no blade stations should be refused");
}

void TestBladeSectionsRideAlong() {
    // Blade CST sections arrived at schema 1.8.0: that fixture must carry
    // them through to the metric propeller (as AirfoilSections keyed by
    // Eta = r/R), and the 1.5.0 fixture, which predates them, must yield
    // flat blades rather than inventing shape data.
    const auto old = Aeolion::Geometry::ToPropeller(Fixture().Propulsion);
    CHECK(old.Sections.empty(), "a pre-1.8.0 contract states no blade sections");

    const auto contract = Aeolion::Geometry::LoadHandoff(std::string(AEOLION_TEST_DATA_DIR) +
                                                         "/AeolionGeometryHandoff-1.8.0.json");
    const auto prop = Aeolion::Geometry::ToPropeller(contract.Propulsion);
    CHECK(prop.Sections.size() == contract.Propulsion.BladeAirfoilSections.size(),
          "every stated blade section should survive the conversion");
    for (std::size_t i = 0; i < prop.Sections.size(); ++i) {
        CHECK(prop.Sections[i].Eta == contract.Propulsion.BladeAirfoilSections[i].RadiusFraction,
              "section Eta is the stated radius fraction");
        CHECK(!prop.Sections[i].CoefficientsUpper.empty() && !prop.Sections[i].CoefficientsLower.empty(),
              "stated sections should carry their CST coefficient sets");
        if (i > 0)
            CHECK(prop.Sections[i].Eta > prop.Sections[i - 1].Eta, "sections must stay ordered by radius");
    }
}

} // namespace

int main() {
    TestConversionConvertsUnits();
    TestConversionRefusesWhatItCannotKnow();
    TestBladeSectionsRideAlong();

    if (failures == 0) { std::cout << "PASS: TestPropeller\n"; return 0; }
    std::cerr << failures << " check(s) failed in TestPropeller\n";
    return 1;
}
