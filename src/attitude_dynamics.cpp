#include "6dof_orbit_sim/attitude_dynamics.hpp"
#include <cmath>
#include <stdexcept>

namespace orb {

AttitudeStateDot computeGravGradAttitudeDynamics(const AttitudeState& attitude_state, const ECIState& eci_state, Mat3& inertia_tensor, Mat3& inertia_tensor_inv, Real mu) {
    Real r = eci_state.position.norm();
    Vec3 r_hat_eci = eci_state.position / r;

    // Rotate into body frame
    Vec3 r_hat_body = attitude_state.q.rotate(r_hat_eci);

    // τ_gg = (3μ/r³) * r̂_body × (I · r̂_body)
    Real coeff  = 3.0 * mu / (r * r * r);
    Vec3 torque = coeff * r_hat_body.cross(inertia_tensor * r_hat_body);

    // --- Euler's rotational equation: ω̇ = I⁻¹(τ - ω × Iω) ---
    Vec3 omega = attitude_state.omega;
    Vec3 alpha = inertia_tensor_inv * (torque - omega.cross(inertia_tensor * omega));

    // --- Quaternion kinematics: q̇ = 0.5 * q ⊗ [0, ω] ---
    QuaternionDot q_dot = quaternion_kinematics(attitude_state.q, omega);

    return AttitudeStateDot{q_dot, alpha};
}

AttitudeStateDot computeDragAttitudeDynamics(const AttitudeState& attitude_state,
                                      const ECIState& eci_state,
                                      const Vec3& center_of_mass,
                                      Mat3& inertia_tensor,
                                      Mat3& inertia_tensor_inv,
                                      Real Cd, Real area) {

    Real rho = 1e-10; // Placeholder - replace with actual atmosphere model
    
    Real v_rel = eci_state.velocity.norm();

    Vec3 force_eci = -0.5 * rho * Cd * area * v_rel * eci_state.velocity;
    Vec3 force_body = attitude_state.q.rotate(force_eci);
    Vec3 torque = (-center_of_mass).cross(force_body); //Assumes center of pressure is at origin of body frame

    Vec3 omega = attitude_state.omega;
    Vec3 alpha = inertia_tensor_inv * (torque - omega.cross(inertia_tensor * omega));

    QuaternionDot q_dot = quaternion_kinematics(attitude_state.q, omega);

    return AttitudeStateDot{q_dot, alpha};
}

AttitudeStateDot computeSRPAttitudeDynamics(const AttitudeState& attitude_state,
                                     const ECIState& eci_state,
                                     const Vec3& sun_dir_eci,
                                     const Vec3& center_of_mass,
                                     Mat3& inertia_tensor,
                                     Mat3& inertia_tensor_inv,
                                     Real Cr, Real area, Real P_srp) {

    Vec3 force_eci = -P_srp * Cr * area * sun_dir_eci;
    Vec3 force_body = attitude_state.q.rotate(force_eci);
    Vec3 torque = (-center_of_mass).cross(force_body);

    Vec3 omega = attitude_state.omega;
    Vec3 alpha = inertia_tensor_inv * (torque - omega.cross(inertia_tensor * omega));

    QuaternionDot q_dot = quaternion_kinematics(attitude_state.q, omega);

    return AttitudeStateDot{q_dot, alpha};
}

AttitudeStateDot computeTotalAttitudeDynamics(const AttitudeState& attitude_state,
                             const ECIState& eci_state,
                             Mat3& inertia_tensor,
                             Mat3& inertia_tensor_inv,
                             const TorqueModelConfig& config,
                             Real mu) {

    Vec3 omega_dot = Vec3::Zero();

    if (config.useGravityGradient) {
        AttitudeStateDot x = computeGravGradAttitudeDynamics(
            attitude_state,
            eci_state,
            inertia_tensor,
            inertia_tensor_inv,
            mu
        );
        omega_dot += x.omega_dot;
    }

    if (config.useDrag) {
        AttitudeStateDot x = computeDragAttitudeDynamics(
            attitude_state,
            eci_state,
            config.center_of_mass,
            inertia_tensor,
            inertia_tensor_inv,
            config.Cd,
            config.area
        );
        omega_dot += x.omega_dot;
    }

    if (config.useSRP) {
        AttitudeStateDot x = computeSRPAttitudeDynamics(
            attitude_state,
            eci_state,
            config.sun_dir_eci,
            config.center_of_mass,
            inertia_tensor,
            inertia_tensor_inv,
            config.Cr,
            config.area,
            config.P_srp
        );
        omega_dot += x.omega_dot;
    }

    QuaternionDot q_dot = quaternion_kinematics(attitude_state.q, attitude_state.omega);

    return AttitudeStateDot{q_dot, omega_dot};
}

DerivFunc makeAttDeriv(const ECIState& eci_state,
                        Mat3& inertia,
                        Mat3& inertia_inv,
                        const TorqueModelConfig& config) {
    return [&inertia, &inertia_inv, config, eci_state](Real t, const VecX& state_vec) -> VecX {
        AttitudeState att = unpackAttitude(state_vec.head<ATTITUDE_STATE_DIM>());
        AttitudeStateDot dot = computeTotalAttitudeDynamics(att, eci_state,
                                                             inertia, inertia_inv, config);
        return packAttitudeDot(dot);
    };
}

} // namespace orb