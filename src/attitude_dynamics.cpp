#include "6dof_orbit_sim/attitude_dynamics.hpp"
#include <cmath>
#include <stdexcept>

namespace orb {

AttitudeStateDot ggradAttitudeDynamics(const AttitudeState& attitude_state, const ECIState& eci_state, Mat3& inertia_tensor, Mat3& inertia_tensor_inv, Real mu) {
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
    Quaternion q_dot = 0.5 * attitude_state.q * Quaternion{0.0, omega.x(), omega.y(), omega.z()};

    return AttitudeStateDot{q_dot, alpha};
}

AttitudeStateDot dragAttitudeDynamics(const AttitudeState& attitude_state,
                                      const ECIState& eci_state,
                                      const Vec3& center_of_mass,
                                      Mat3& inertia_tensor,
                                      Mat3& inertia_tensor_inv,
                                      Real Cd, Real area) {

    Real rho = 1e-10; // Placeholder - replace with actual atmosphere model
    
    Real v_rel = eci_state.velocity.norm();

    Vec3 force_eci = -0.5 * rho * Cd * area * v_rel * eci_state.velocity.normalized();
    Vec3 force_body = attitude_state.q.rotate(force_eci);
    Vec3 torque = (-center_of_mass).cross(force_body);

    Vec3 omega = attitude_state.omega;
    Vec3 alpha = inertia_tensor_inv * (torque - omega.cross(inertia_tensor * omega));

    Quaternion q_dot = 0.5 * attitude_state.q * Quaternion{0.0, omega.x(), omega.y(), omega.z()};

    return AttitudeStateDot{q_dot, alpha};
}

AttitudeStateDot srpAttitudeDynamics(const AttitudeState& attitude_state,
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

    Quaternion q_dot = 0.5 * attitude_state.q * Quaternion{0.0, omega.x(), omega.y(), omega.z()};

    return AttitudeStateDot{q_dot, alpha};
}

} // namespace orb