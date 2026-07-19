#include "Aeolion/GeometryContract.h"
#include "Aeolion/VLM.h"
#include <print>

int main(int argc, char** argv) {
    if (argc != 2) {
        std::println(stderr, "usage: aeolion_geometry <aeolion_geometry.json>");
        return 2;
    }
    try {
        const auto geometry = Aeolion::Geometry::Load(argv[1]);
        Aeolion::VLM::FreestreamConditions conditions;
        const auto result = Aeolion::VLM::Solve(geometry.Wing, conditions);
        std::print("surface={}\n"
                   "airfoil={}\n"
                   "panels={}\n"
                   "reference_area_m2={:.10g}\n"
                   "CL={:.10g}\n"
                   "CDi={:.10g}\n",
                   geometry.SurfaceName, geometry.Airfoil, result.gamma.size(),
                   result.ReferenceArea, result.CL, result.CDi);
    } catch (const std::exception& error) {
        std::println(stderr, "geometry error: {}", error.what());
        return 1;
    }
}
