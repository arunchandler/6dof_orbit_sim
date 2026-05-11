#pragma once
#include "types.hpp"
#include "translational_state.hpp"
#include "attitude_state.hpp"
#include "state_conversions.hpp"

namespace orb {

/// @brief Compute the time derivative of the attitude state from gravity gradient
/// @param attitude_state Current attitude state (quaternion and angular velocity)
/// @param eci_state Current state in ECI frame (position and velocity)
/// @param inertia_tensor Inertia tensor in body frame [kg·m²]
/// @param inertia_tensor_inv Inverse of inertia tensor for dynamics calculations
/// @param mu Gravitational parameter (default: Earth's mu)
/// @return  Time derivative of attitude state (q_dot and omega_dot)
AttitudeStateDot computeGravGradAttitudeDynamics(const AttitudeState& attitude_state, const ECIState& eci_state, Mat3& inertia_tensor, Mat3& inertia_tensor_inv, Real mu = constants::MU_EARTH);

/// @brief Compute the time derivative of the attitude state from aerodynamic torque
/// @param attitude_state Current attitude state (quaternion and angular velocity)
/// @param Cd Drag coefficient (default: 2.2)
/// @param A Reference area (default: 1.0)
/// @param m Mass (default: 100.0)
/// @return  Time derivative of attitude state (q_dot and omega_dot)
AttitudeStateDot computeDragAttitudeDynamics(const AttitudeState& attitude_state,
                                      const ECIState& eci_state,
                                      const Vec3& center_of_mass,
                                      Mat3& inertia_tensor,
                                      Mat3& inertia_tensor_inv,
                                      Real Cd, Real area);

/// @brief Compute the time derivative of the attitude state from solar radiation pressure torque
/// @param attitude_state Current attitude state (quaternion and angular velocity)
/// @param Cr Reflectivity coefficient (default: 1.5)
/// @param A Reference area (default: 1.0)
/// @param m Mass (default: 100.0)
/// @return  Time derivative of attitude state (q_dot and omega_dot)
AttitudeStateDot computeSRPAttitudeDynamics(const AttitudeState& attitude_state,
                                     const ECIState& eci_state,
                                     const Vec3& sun_dir_eci,
                                     const Vec3& center_of_mass,
                                     Mat3& inertia_tensor,
                                     Mat3& inertia_tensor_inv,
                                     Real Cr, Real area, Real P_srp);

/// @brief Configuration for selecting which torques to include in total attitude dynamics.
struct TorqueModelConfig {
    bool useGravityGradient = true;
    bool useDrag            = false;
    bool useSRP             = false;

    // Drag parameters
    Real Cd   = 2.2;
    Real area = 1.0;   // m^2

    // SRP parameters
    Real Cr    = 1.5;
    Real P_srp = 4.56e-6; // N/m^2 at 1 AU

    // Shared geometry
    Vec3 center_of_mass = Vec3::Zero();

    // Required for SRP if enabled
    Vec3 sun_dir_eci = Vec3::UnitX();
};

/// @brief Compute total attitude state derivative from selected torque models.
/// @param attitude_state Current attitude state
/// @param eci_state Current translational state in ECI
/// @param inertia_tensor Inertia tensor in body frame
/// @param inertia_tensor_inv Inverse inertia tensor
/// @param config Torque model configuration
/// @param mu Gravitational parameter
/// @return Total attitude state derivative (q_dot, omega_dot)
AttitudeStateDot computeTotalAttitudeDynamics(const AttitudeState& attitude_state,
                             const ECIState& eci_state,
                             Mat3& inertia_tensor,
                             Mat3& inertia_tensor_inv,
                             const TorqueModelConfig& config,
                             Real mu = constants::MU_EARTH);

DerivFunc makeAttDeriv(const ECIState& eci_state,
                        Mat3& inertia,
                        Mat3& inertia_inv,
                        const TorqueModelConfig& config);

} // namespace orb