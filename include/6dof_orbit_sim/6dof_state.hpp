#pragma once
#include <array>
#include <cmath>
#include "types.hpp"
#include "attitude_state.hpp"
#include "translational_state.hpp"

namespace orb {

struct SixDoFState {
    ECIState eci;           // Position and velocity in ECI frame
    AttitudeState attitude; // Orientation and angular velocity
};

struct SixDoFStateDot {
    ECIStateDot eci_dot;           // Time derivative of ECI state
    AttitudeStateDot attitude_dot; // Time derivative of attitude state
};

} // namespace orb