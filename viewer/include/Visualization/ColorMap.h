// Visualization/ColorMap.h
//
// Viridis colormap for scalar fields on the lattice: perceptually uniform,
// colorblind-safe, and the de-facto standard in aero/CFD post-processing.
// Polynomial fit by Matt Zucker (public domain, shadertoy 'viridis quintic'
// family) -- accurate to well under 1% over [0,1].

#pragma once

#include <glm/glm.hpp>

#include <algorithm>

namespace Aeolion::Viewer {

[[nodiscard]] inline glm::vec3 Viridis(float t) {
    t = std::clamp(t, 0.0f, 1.0f);
    const glm::vec3 c0(0.2777273272234177f, 0.005407344544966578f, 0.3340998053353061f);
    const glm::vec3 c1(0.1050930431085774f, 1.404613529898575f, 1.384590162594685f);
    const glm::vec3 c2(-0.3308618287255563f, 0.214847559468213f, 0.09509516302823659f);
    const glm::vec3 c3(-4.634230498983486f, -5.799100973351585f, -19.33244095627987f);
    const glm::vec3 c4(6.228269936347081f, 14.17993336680509f, 56.69055260068105f);
    const glm::vec3 c5(4.776384997670288f, -13.74514537774601f, -65.35303263337234f);
    const glm::vec3 c6(-5.435455855934631f, 4.645852612178535f, 26.3124352495832f);
    return c0 + t * (c1 + t * (c2 + t * (c3 + t * (c4 + t * (c5 + t * c6)))));
}

} // namespace Aeolion::Viewer
