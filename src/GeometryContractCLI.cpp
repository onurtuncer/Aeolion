#include "Aeolion/GeometryContract.h"
#include "Aeolion/VLM.h"
#include <iomanip>
#include <iostream>

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: aeolion_geometry <aeolion_geometry.json>\n";
        return 2;
    }
    try {
        const auto geometry = Aeolion::Geometry::Load(argv[1]);
        Aeolion::VLM::FreestreamConditions conditions;
        const auto result = Aeolion::VLM::Solve(geometry.Wing, conditions);
        std::cout << std::setprecision(10)
                  << "surface=" << geometry.SurfaceName << "\n"
                  << "airfoil=" << geometry.Airfoil << "\n"
                  << "panels=" << result.gamma.size() << "\n"
                  << "reference_area_m2=" << result.ReferenceArea << "\n"
                  << "CL=" << result.CL << "\n"
                  << "CDi=" << result.CDi << "\n";
    } catch (const std::exception& error) {
        std::cerr << "geometry error: " << error.what() << "\n";
        return 1;
    }
}
