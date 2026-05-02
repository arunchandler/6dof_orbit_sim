#pragma once
#include "types.hpp"
#include "attitude_state.hpp"

namespace orb {

/// @brief Compute the time derivative of the attitude state from gravity gradient
/// @param state Current attitude state (quaternion and angular velocity)
/// @param mu Gravitational parameter (default: Earth's mu)
/// @param R Reference radius (default: Earth's mean radius)
/// @return  Time derivative of attitude state (q_dot and omega_dot)
AttitudeStateDot ggradAttitudeDynamics(const AttitudeState& state, Real mu = constants::MU_EARTH, Real R = constants::R_EARTH);

/// @brief Compute the time derivative of the attitude state from aerodynamic torque
/// @param state Current attitude state (quaternion and angular velocity)
/// @param Cd Drag coefficient (default: 2.2)
/// @param A Reference area (default: 1.0)
/// @param m Mass (default: 100.0)
/// @return  Time derivative of attitude state (q_dot and omega_dot)
AttitudeStateDot aeroAttitudeDynamics(const AttitudeState& state, Real Cd = 2.2, Real A = 1.0, Real m = 100.0);

/// @brief Compute the time derivative of the attitude state from solar radiation pressure torque
/// @param state Current attitude state (quaternion and angular velocity)
/// @param Cr Reflectivity coefficient (default: 1.5)
/// @param A Reference area (default: 1.0)
/// @param m Mass (default: 100.0)
/// @return  Time derivative of attitude state (q_dot and omega_dot)
AttitudeStateDot srpAttitudeDynamics(const AttitudeState& state, Real Cr = 1.5, Real A = 1.0, Real m = 100.0);

} // namespace orb