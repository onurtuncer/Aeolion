// The supplier that completes C1: stagnation strain tabulated per (eta,
// alpha) so the camber-line boundary-layer march can be seeded with the
// momentum thickness its unresolved stagnation region would have delivered.
//
// The properties checked are the ones a wrong implementation would still
// satisfy a panel count on: that the strain is POSITIVE and finite where the
// march uses it, that a sharper nose strains harder than a blunt one (which
// is what says the number is a nose property and not an artifact of the
// panelling), that camber makes it ASYMMETRIC in incidence (the reason this
// table is signed where the separation tables are not), and that the edge
// rule clamps rather than extrapolating into negative strain, for which the
// Hiemenz seed has no real square root.
#include "Aeolion/Solver/StagnationStrainTables.h"

#include <cmath>
#include <iostream>
#include <vector>

namespace S = Aeolion::Solver;
namespace G = Aeolion::Geometry;

static int failures = 0;
#define CHECK(cond, msg) \
    do { if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } } while (0)

namespace {

// CST: first coefficient sets the leading-edge radius as r_LE = A0^2 / 2.
G::AirfoilSection Symmetric(double eta, double a0) {
    G::AirfoilSection s;
    s.Eta = eta;
    s.CoefficientsUpper = {a0, a0, a0};
    s.CoefficientsLower = {-a0, -a0, -a0};
    return s;
}

G::AirfoilSection Cambered(double eta, double a0, double camber) {
    G::AirfoilSection s;
    s.Eta = eta;
    s.CoefficientsUpper = {a0 + camber, a0 + camber, a0 + camber};
    s.CoefficientsLower = {-a0 + camber, -a0 + camber, -a0 + camber};
    return s;
}

} // namespace

int main() {
    // A blunt section and a sharp one, so the nose-radius trend is testable.
    const std::vector<G::AirfoilSection> sections{Symmetric(0.0, 0.30), Symmetric(1.0, 0.15)};
    const auto tables = S::BuildStagnationStrainTables(sections);
    CHECK(tables.size() == 2, "one table per section");
    CHECK(!tables[0].AlphaDeg.empty(), "the blunt section produced no table");

    const auto strain = S::MakeStagnationStrainFunction(tables, {0.0, 1.0});

    // 1. Positive and finite across the envelope the march actually uses.
    //    A zero here silently disables the seed; a negative one has no
    //    Hiemenz solution at all.
    for (double a = -13.0; a <= 13.0; a += 1.0) {
        const double v = strain(0.0, a);
        CHECK(v > 0.0 && std::isfinite(v), "strain must be positive and finite at alpha " << a
                                                                                          << ": " << v);
    }

    // 2. A SHARPER nose strains harder. r_LE = A0^2/2, so the a0 = 0.15
    //    section has a quarter the leading-edge radius of the a0 = 0.30 one,
    //    and the flow must turn correspondingly faster around it. This is the
    //    check that says the number is a property of the nose rather than of
    //    the discretization.
    const double blunt = strain(0.0, 0.0);
    const double sharp = strain(1.0, 0.0);
    CHECK(sharp > blunt, "the sharper nose must strain harder: blunt " << blunt << " vs sharp "
                                                                       << sharp);

    // 3. A SYMMETRIC section must be symmetric in incidence -- the control
    //    that stops check 4 from passing on a table that is simply noisy.
    const double plus = strain(0.0, 6.0), minus = strain(0.0, -6.0);
    CHECK(std::fabs(plus - minus) < 0.02 * std::max(plus, minus),
          "a symmetric section must strain equally either way: " << plus << " vs " << minus);

    // 4. A CAMBERED section must NOT be. This is the whole reason the table is
    //    signed where the separation tables key on |alpha|: at +6 and -6 the
    //    stagnation point sits on opposite surfaces at different curvatures.
    const std::vector<G::AirfoilSection> camberedSections{Cambered(0.5, 0.25, 0.12)};
    const auto camberedTables = S::BuildStagnationStrainTables(camberedSections);
    const auto camberedStrain = S::MakeStagnationStrainFunction(camberedTables, {0.5});
    const double cPlus = camberedStrain(0.5, 6.0), cMinus = camberedStrain(0.5, -6.0);
    CHECK(std::fabs(cPlus - cMinus) > 0.02 * std::max(cPlus, cMinus),
          "camber must break the incidence symmetry, got " << cPlus << " vs " << cMinus);

    // 5. The edge rule CLAMPS. Linear extrapolation of a falling strain
    //    eventually goes negative, and sqrt(0.075 nu / a) then has no real
    //    value. Separation tables extrapolate deliberately; this must not.
    const double far = strain(0.0, 400.0);
    CHECK(far >= 0.0, "strain must never go negative, got " << far << " at alpha 400");
    CHECK(std::fabs(far - strain(0.0, S::StrainTableMaxAlphaDeg)) < 1e-12,
          "beyond the table the strain must clamp to its last value");

    // 6. An empty table means no seed, which reproduces the unseeded march.
    const auto none = S::MakeStagnationStrainFunction({}, {});
    CHECK(none(0.5, 4.0) == 0.0, "an empty table must report zero strain, i.e. no seed");

    if (failures == 0) {
        std::cout << "TestStagnationStrainTables: blunt " << blunt << ", sharp " << sharp
                  << ", cambered " << cPlus << "/" << cMinus << " -- OK" << std::endl;
        return 0;
    }
    std::cerr << failures << " check(s) failed\n";
    return 1;
}
