#pragma once
#include "types.hpp"
#include "translational_state.hpp"
#include "attitude_state.hpp"
#include "6dof_state.hpp"
#include "translational_dynamics.hpp"
#include "attitude_dynamics.hpp"

namespace orb {

/// @brief Compute the time derivative of the full 6DoF state (translational + attitude) based on selected torque models.
/// @param state Current 6DoF state (position, velocity, attitude, angular velocity)
/// @param inertia_tensor Inertia tensor in body frame [kg·m²]
/// @param inertia_tensor_inv Inverse of inertia tensor for dynamics calculations
/// @param force_config Configuration for which forces to include in translational dynamics
/// @param torque_config Configuration for which torques to include in attitude dynamics
/// @return Time derivative of the full 6DoF state (position_dot, velocity_dot, q_dot, omega_dot)
SixDoFStateDot sixDoFDynamics(const SixDoFState& state,
                             Mat3& inertia_tensor,
                             Mat3& inertia_tensor_inv,
                             const ForceModelConfig& force_config,
                             const TorqueModelConfig& torque_config,
                             ActuatorSuite<4,3>& suite,
                             ActuatorCommands<4,3>& cmds,
                             const Quaternion& q_desired,
                             const QuaternionPDController& ctrl,
                             Real t_sec,
                             Real dt);

DerivFunc makeSixDoFDeriv(Mat3& inertia, Mat3& inertia_inv,
                           const ForceModelConfig& force_config,
                           const TorqueModelConfig& torque_config,
                           ActuatorSuite<4,3>& suite,
                           ActuatorCommands<4,3>& cmds,
                           const Quaternion& q_desired,
                           const QuaternionPDController& ctrl,
                           Real dt);

} // namespace orb