#include "6dof_orbit_sim/actuator.hpp"
 
#include <stdexcept>
 
namespace orb {
 
// ─────────────────────────────────────────────────────────────────────────────
// Magnetic field model
// ─────────────────────────────────────────────────────────────────────────────
Vec3 dipole_field_eci(const Vec3& r_eci, Real t_sec)
{
    const Real theta    = constants::EARTH_ROT_RATE * t_sec
                        + constants::MAG_DIPOLE_LON_J2000;
    const Real sin_tilt = std::sin(constants::MAG_DIPOLE_TILT);
    const Real cos_tilt = std::cos(constants::MAG_DIPOLE_TILT);
 
    const Vec3 m_hat(sin_tilt * std::cos(theta),
                     sin_tilt * std::sin(theta),
                     cos_tilt);
 
    const Real r = r_eci.norm();
    if (r < Real(1))
        throw std::domain_error("dipole_field_eci: |r_eci| is near zero");
 
    const Vec3 r_hat = r_eci / r;
    const Real scale = constants::MU0_OVER_4PI
                     * constants::EARTH_MAG_MOM / (r * r * r);
 
    return scale * (Real(3) * m_hat.dot(r_hat) * r_hat - m_hat);
}
 
// ─────────────────────────────────────────────────────────────────────────────
// ReactionWheel
// ─────────────────────────────────────────────────────────────────────────────
ActuatorOutput ReactionWheel::update(Real tau_cmd, Real dt)
{
    const Real tau     = std::min(std::max(tau_cmd, -max_torque), max_torque);
    const Real tau_net = tau - friction_coef * omega_w;
 
    omega_w += (tau_net / inertia) * dt;
    omega_w = std::min(std::max(omega_w, -max_speed), max_speed);
 
    ActuatorOutput out;
    out.torque_b = -tau * spin_axis_b;
    return out;
}
 
Vec3 ReactionWheel::momentum() const
{
    return inertia * omega_w * spin_axis_b;
}
 
bool ReactionWheel::is_saturated() const
{
    return std::abs(omega_w) >= max_speed * Real(0.999);
}
 
// ─────────────────────────────────────────────────────────────────────────────
// Thruster
// ─────────────────────────────────────────────────────────────────────────────
ActuatorOutput Thruster::fire(Real duty_cycle) const
{
    const Vec3 F = (max_thrust * std::min(std::max(duty_cycle, Real(0)), Real(1)))
                 * direction_b;
    ActuatorOutput out;
    out.force_b  = F;
    out.torque_b = position_b.cross(F);
    return out;
}
 
// ─────────────────────────────────────────────────────────────────────────────
// ThrusterPair
// ─────────────────────────────────────────────────────────────────────────────
ActuatorOutput ThrusterPair::fire(Real torque_cmd) const
{
    if (torque_cmd > Real(0)) return pos_thruster.fire(Real(1));
    if (torque_cmd < Real(0)) return neg_thruster.fire(Real(1));
    return {};
}
 
// ─────────────────────────────────────────────────────────────────────────────
// Magnetorquer
// ─────────────────────────────────────────────────────────────────────────────
ActuatorOutput Magnetorquer::apply(Real dipole_cmd, const Vec3& B_body) const
{
    const Real m_mag = max_dipole * std::min(std::max(dipole_cmd, Real(-1)), Real(1));
    ActuatorOutput out;
    out.torque_b = (m_mag * axis_b).cross(B_body);
    return out;
}
 
Real Magnetorquer::alignment_with_field(const Vec3& B_body) const
{
    const Real Bnorm = B_body.norm();
    return (Bnorm < Real(1e-15)) ? Real(0)
                                 : std::abs(axis_b.dot(B_body / Bnorm));
}
 
// ─────────────────────────────────────────────────────────────────────────────
// MagnetorquerAssembly
// ─────────────────────────────────────────────────────────────────────────────
ActuatorOutput MagnetorquerAssembly::apply_torque_cmd(const Vec3& tau_cmd_b,
                                                      const Vec3& B_body) const
{
    const Real Bsq = B_body.squaredNorm();
    if (Bsq < Real(1e-20)) return {};
 
    const Vec3 m_desired = B_body.cross(tau_cmd_b) / Bsq;
 
    ActuatorOutput total;
    for (int i = 0; i < 3; ++i) {
        const Real cmd = m_desired.dot(mtqs[i].axis_b) / mtqs[i].max_dipole;
        total += mtqs[i].apply(cmd, B_body);
    }
    return total;
}
 
ActuatorOutput MagnetorquerAssembly::apply_dipole_cmd(const Vec3& dipole_cmd_b,
                                                      const Vec3& B_body) const
{
    ActuatorOutput total;
    for (int i = 0; i < 3; ++i) {
        const Real cmd = dipole_cmd_b.dot(mtqs[i].axis_b) / mtqs[i].max_dipole;
        total += mtqs[i].apply(cmd, B_body);
    }
    return total;
}
 
MagnetorquerAssembly MagnetorquerAssembly::make_xyz(Real max_dipole_Am2)
{
    return { {{
        { max_dipole_Am2, Vec3::UnitX() },
        { max_dipole_Am2, Vec3::UnitY() },
        { max_dipole_Am2, Vec3::UnitZ() }
    }} };
}

Quaternion QuaternionPDController::error_quaternion(const AttitudeState& state,
                                                    const Quaternion&    q_desired)
{
    // q_e = q_desired* ⊗ q_current
    Quaternion q_e = q_desired.conjugate() * state.q;
 
    // Always take the short arc: if scalar part is negative, negate the whole
    // quaternion (represents the same rotation, just the long-way round)
    if (q_e.w < Real(0)) {
        q_e.w = -q_e.w;
        q_e.x = -q_e.x;
        q_e.y = -q_e.y;
        q_e.z = -q_e.z;
    }
 
    return q_e;
}
 
Vec3 QuaternionPDController::compute_torque(const AttitudeState& state,
                                            const Quaternion&    q_desired) const
{
    const Quaternion q_e = error_quaternion(state, q_desired);
 
    // Vector part of error quaternion ≈ axis * sin(θ/2) ≈ axis * θ/2 for small errors
    const Vec3 q_e_vec(q_e.x, q_e.y, q_e.z);
 
    return -k_p * q_e_vec - k_d * state.omega;
}

Quaternion make_desired_quaternion(const Vec3& primary_axis_b,
                                   const Vec3& primary_target,
                                   const Vec3& secondary_axis_b,
                                   const Vec3& secondary_target)
{
    // --- Build desired DCM column-by-column (body axes expressed in ECI) ---
 
    // First column of R_eci_from_body: primary body axis maps to primary target
    const Vec3 col1 = primary_target.normalized();
 
    // Orthogonalise secondary target against primary, then normalise
    const Vec3 sec_orth = secondary_target - secondary_target.dot(col1) * col1;
    if (sec_orth.norm() < Real(1e-10))
        throw std::runtime_error(
            "make_desired_quaternion: secondary target is parallel to primary");
 
    // Third body axis (cross of primary and secondary body axes, in ECI)
    // Determine which ECI direction the (primary × secondary) body axis maps to
    const Vec3 tertiary_axis_b = primary_axis_b.cross(secondary_axis_b);
 
    const Vec3 col2 = sec_orth.normalized();                 // secondary ECI direction
    const Vec3 col3 = col1.cross(col2);                      // tertiary ECI direction
 
    // Map body axes to ECI directions:
    //   R_eci_from_body * primary_axis_b   = col1
    //   R_eci_from_body * secondary_axis_b = col2  (approx, after orthogonalisation)
    //   R_eci_from_body * tertiary_axis_b  = col3
    //
    // Assemble R by solving: R * B = E  →  R = E * B^-1 = E * B^T  (B is orthonormal)
    Eigen::Matrix<Real,3,3> B, E;
    B.col(0) = primary_axis_b;
    B.col(1) = secondary_axis_b;
    B.col(2) = tertiary_axis_b;
    E.col(0) = col1;
    E.col(1) = col2;
    E.col(2) = col3;
 
    const Eigen::Matrix<Real,3,3> R = E * B.transpose();
 
    // --- DCM → quaternion (Shepperd's method) ---
    const Real trace = R.trace();
    Real w, x, y, z;
 
    if (trace > Real(0)) {
        const Real s = Real(0.5) / std::sqrt(trace + Real(1));
        w = Real(0.25) / s;
        x = (R(2,1) - R(1,2)) * s;
        y = (R(0,2) - R(2,0)) * s;
        z = (R(1,0) - R(0,1)) * s;
    } else if (R(0,0) > R(1,1) && R(0,0) > R(2,2)) {
        const Real s = Real(2) * std::sqrt(Real(1) + R(0,0) - R(1,1) - R(2,2));
        w = (R(2,1) - R(1,2)) / s;
        x = Real(0.25) * s;
        y = (R(0,1) + R(1,0)) / s;
        z = (R(0,2) + R(2,0)) / s;
    } else if (R(1,1) > R(2,2)) {
        const Real s = Real(2) * std::sqrt(Real(1) + R(1,1) - R(0,0) - R(2,2));
        w = (R(0,2) - R(2,0)) / s;
        x = (R(0,1) + R(1,0)) / s;
        y = Real(0.25) * s;
        z = (R(1,2) + R(2,1)) / s;
    } else {
        const Real s = Real(2) * std::sqrt(Real(1) + R(2,2) - R(0,0) - R(1,1));
        w = (R(1,0) - R(0,1)) / s;
        x = (R(0,2) + R(2,0)) / s;
        y = (R(1,2) + R(2,1)) / s;
        z = Real(0.25) * s;
    }
 
    Quaternion q{w, x, y, z};
    q.normalize();
    return q;
}
 
} // namespace orb