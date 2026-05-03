#pragma once
#include "types.hpp"
#include "translational_state.hpp"
#include "attitude_state.hpp"

namespace orb {

/// @brief Compute the time derivative of the attitude state from gravity gradient
/// @param attitude_state Current attitude state (quaternion and angular velocity)
/// @param eci_state Current state in ECI frame (position and velocity)
/// @param inertia_tensor Inertia tensor in body frame [kg·m²]
/// @param inertia_tensor_inv Inverse of inertia tensor for dynamics calculations
/// @param mu Gravitational parameter (default: Earth's mu)
/// @return  Time derivative of attitude state (q_dot and omega_dot)
AttitudeStateDot ggradAttitudeDynamics(const AttitudeState& attitude_state, const ECIState& eci_state, Mat3& inertia_tensor, Mat3& inertia_tensor_inv, Real mu = constants::MU_EARTH);

/// @brief Compute the time derivative of the attitude state from aerodynamic torque
/// @param attitude_state Current attitude state (quaternion and angular velocity)
/// @param Cd Drag coefficient (default: 2.2)
/// @param A Reference area (default: 1.0)
/// @param m Mass (default: 100.0)
/// @return  Time derivative of attitude state (q_dot and omega_dot)
AttitudeStateDot aeroAttitudeDynamics(const AttitudeState& attitude_state, Real Cd = 2.2, Real A = 1.0, Real m = 100.0);

/// @brief Compute the time derivative of the attitude state from solar radiation pressure torque
/// @param attitude_state Current attitude state (quaternion and angular velocity)
/// @param Cr Reflectivity coefficient (default: 1.5)
/// @param A Reference area (default: 1.0)
/// @param m Mass (default: 100.0)
/// @return  Time derivative of attitude state (q_dot and omega_dot)
AttitudeStateDot srpAttitudeDynamics(const AttitudeState& attitude_state, Real Cr = 1.5, Real A = 1.0, Real m = 100.0);

} // namespace orb