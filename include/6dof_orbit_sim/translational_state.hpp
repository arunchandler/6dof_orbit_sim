#pragma once
#include "types.hpp"

namespace orb {

/// Translational state in Earth-Centered Inertial frame
struct ECIState {
    Vec3 position;   // r [m], from Earth center
    Vec3 velocity;   // v [m/s], inertial frame
};

struct ECIStateDot {
    Vec3 velocity;      // dr/dt = v
    Vec3 acceleration;  // dv/dt = f/m
};

} // namespace orb