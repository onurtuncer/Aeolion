// TestParticleTree.cpp -- the Barnes-Hut particle tree
// (Solver/ParticleTree.h), pinned by EQUIVALENCE against direct summation:
// the tree is an evaluation strategy, not a model, so its one obligation
// is to reproduce the direct sum -- exactly when every cell is opened,
// and within the opening criterion's documented error band otherwise, for
// the velocity AND the gradient, with self-exclusion intact.

#include "Aeolion/Solver/ParticleTree.h"

#include <chrono>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

using namespace Aeolion;
namespace S = Aeolion::Solver;

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } } while (0)

namespace {

struct Cloud {
    S::Vec3 X, Alpha;
    double Core2;
};

std::vector<Cloud> MakeCloud(int count, unsigned seed) {
    // A street-like cloud: elongated slab, mixed strengths and signs,
    // cores spanning a spreading wake's range.
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> ux(0.0, 8.0), uy(-1.0, 1.0), uz(-0.6, 0.6);
    std::uniform_real_distribution<double> ua(-1.0, 1.0), uc(0.01, 0.06);
    std::vector<Cloud> cloud(static_cast<std::size_t>(count));
    for (Cloud& p : cloud) {
        p.X = S::Vec3(ux(rng), uy(rng), uz(rng));
        p.Alpha = S::Vec3(ua(rng), ua(rng), ua(rng)) * 0.02;
        const double c = uc(rng);
        p.Core2 = c * c;
    }
    return cloud;
}

void DirectSum(const std::vector<Cloud>& cloud, const S::Vec3& at, double targetCore2,
               int exclude, S::Vec3& u, S::Vec3 grad[3]) {
    for (std::size_t k = 0; k < cloud.size(); ++k) {
        if (static_cast<int>(k) == exclude) continue;
        S::Detail::TreeKernel(at - cloud[k].X, cloud[k].Alpha,
                              0.5 * (targetCore2 + cloud[k].Core2), u, grad, true);
    }
}

void TestEquivalence() {
    const auto cloud = MakeCloud(3000, 7);
    S::ParticleTree tree;
    tree.Build(cloud);

    // Targets: a mix of free-space points and points ON particles (the
    // self-exclusion path), the latter with their own index excluded.
    std::mt19937 rng(11);
    std::uniform_int_distribution<int> pick(0, 2999);

    for (double theta : {0.0, 0.25, 0.5}) {
        tree.Theta = theta;
        double worstU = 0.0, worstG = 0.0, refU = 0.0, refG = 0.0;
        for (int t = 0; t < 60; ++t) {
            const int self = (t % 2 == 0) ? pick(rng) : -1;
            const S::Vec3 at = (self >= 0)
                                   ? cloud[static_cast<std::size_t>(self)].X
                                   : S::Vec3(8.0 * (t / 60.0), 0.3, 0.1);
            const double tc2 = (self >= 0) ? cloud[static_cast<std::size_t>(self)].Core2 : 0.0;

            S::Vec3 ud(0, 0, 0), gd[3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
            DirectSum(cloud, at, tc2, self, ud, gd);
            S::Vec3 ut(0, 0, 0), gt[3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
            tree.Evaluate(at, tc2, self, ut, gt);

            worstU = std::max(worstU, (ut - ud).Norm());
            refU = std::max(refU, ud.Norm());
            for (int l = 0; l < 3; ++l) {
                worstG = std::max(worstG, (gt[l] - gd[l]).Norm());
                refG = std::max(refG, gd[l].Norm());
            }
        }
        const double relU = worstU / std::max(refU, 1e-12);
        const double relG = worstG / std::max(refG, 1e-12);
        std::cout << "theta " << theta << ": rel err u " << relU << ", grad " << relG << "\n";
        if (theta == 0.0) {
            // Fully opened, the tree IS the direct sum.
            CHECK(relU < 1e-12, "theta = 0 reproduces direct summation exactly (velocity)");
            CHECK(relG < 1e-12, "theta = 0 reproduces direct summation exactly (gradient)");
        } else if (theta == 0.25) {
            CHECK(relU < 5e-3, "theta = 0.25 velocity inside the tight band");
            CHECK(relG < 2e-2, "theta = 0.25 gradient inside the tight band");
        } else {
            CHECK(relU < 2e-2, "theta = 0.5 velocity inside the working band");
            CHECK(relG < 8e-2, "theta = 0.5 gradient inside the working band");
        }
    }
}

void TestScaling() {
    // Not a benchmark, a sanity ratio: at N = 6000 the tree evaluation
    // over all particles must beat direct by a clear factor.
    const auto cloud = MakeCloud(6000, 23);
    S::ParticleTree tree;
    tree.Build(cloud);
    tree.Theta = 0.5;

    const auto t0 = std::chrono::steady_clock::now();
    S::Vec3 sink(0, 0, 0);
    for (std::size_t i = 0; i < cloud.size(); ++i) {
        S::Vec3 u(0, 0, 0), g[3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        tree.Evaluate(cloud[i].X, cloud[i].Core2, static_cast<int>(i), u, g);
        sink = sink + u;
    }
    const auto t1 = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < 600; ++i) { // a tenth of the targets, directly
        S::Vec3 u(0, 0, 0), g[3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        DirectSum(cloud, cloud[i].X, cloud[i].Core2, static_cast<int>(i), u, g);
        sink = sink + u;
    }
    const auto t2 = std::chrono::steady_clock::now();
    const double treeMs = std::chrono::duration<double, std::milli>(t1 - t0).count();
    const double directMs = 10.0 * std::chrono::duration<double, std::milli>(t2 - t1).count();
    std::cout << "N=6000 all-target sweep: tree " << treeMs << " ms vs direct ~" << directMs
              << " ms (sink " << sink.Norm() << ")\n";
    CHECK(treeMs < 0.5 * directMs, "the tree beats direct summation clearly at N = 6000");
}

} // namespace

int main() {
    TestEquivalence();
    TestScaling();

    if (failures) {
        std::cerr << failures << " check(s) failed in TestParticleTree\n";
        return 1;
    }
    std::cout << "TestParticleTree: all checks passed\n";
    return 0;
}
