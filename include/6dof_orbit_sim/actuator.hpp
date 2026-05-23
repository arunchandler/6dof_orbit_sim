#pragma once
 
/**
 * actuator.hpp
 *
 * Attitude actuator models for 6DoF spacecraft simulation.
 * Includes: reaction wheels, thrusters, magnetorquers, and a
 * tilted-dipole magnetic field model (ECI frame).
 *
 * Conventions:
 *   - All vectors body-frame unless suffixed _eci
 *   - SI units throughout (kg, m, s, N, T, A·m²)
 *   - Actuator omega_dot is ADDITIVE to disturbance omega_dot:
 *
 *       AttitudeStateDot d_dist = computeTotalAttitudeDynamics(...);
 *       AttitudeStateDot d_act  = computeActuatorAttitudeDynamics(...);
 *       omega_dot_total = d_dist.omega_dot + d_act.omega_dot   ← just Vec3 +=
 *       q_dot_total     = d_dist.q_dot                         ← unchanged
 *
 *     This decomposition is exact because:
 *       I dω/dt = (τ_dist − ω×Iω)  +  (τ_act − ω×h_rw)
 *                 ↑ existing funcs       ↑ this file
 */
 
#include <array>
#include <cmath>
#include <stdexcept>
 
#include "types.hpp"
#include "attitude_state.hpp"
 
namespace orb {
 
// ─────────────────────────────────────────────────────────────────────────────
// Magnetic field model  —  IGRF tilted co-rotating dipole  (ECI output)
//
//   r_eci  : spacecraft ECI position [m]  ← pass eci_state.position
//   t_sec  : seconds since J2000.0
//   return : magnetic field vector in ECI [T]
// ─────────────────────────────────────────────────────────────────────────────
Vec3 dipole_field_eci(const Vec3& r_eci, Real t_sec);
 
// ─────────────────────────────────────────────────────────────────────────────
// Actuator output  —  torque (and optionally force) produced in one step.
// Accumulate with += across all subsystems before applying.
// ─────────────────────────────────────────────────────────────────────────────
struct ActuatorOutput {
    Vec3 torque_b = Vec3::Zero();  // Body-frame torque [N·m]
    Vec3 force_b  = Vec3::Zero();  // Body-frame force  [N]  (thrusters only)
 
    ActuatorOutput& operator+=(const ActuatorOutput& rhs) noexcept {
        torque_b += rhs.torque_b;
        force_b  += rhs.force_b;
        return *this;
    }
};
 
// ─────────────────────────────────────────────────────────────────────────────
// REACTION WHEEL
//
// State: omega_w — wheel angular velocity [rad/s]
// Input: tau_cmd — commanded motor torque [N·m]
//
// Wheel EOM:  I_w · dω_w/dt = τ_cmd − τ_friction
// Body torque: −τ_cmd  (Newton's 3rd law)
//
// The stored angular momentum h_rw must be included in the gyroscopic term
// of Euler's equation — see computeActuatorAttitudeDynamics() below.
// ─────────────────────────────────────────────────────────────────────────────
struct ReactionWheel {
    Real inertia;       // Wheel moment of inertia [kg·m²]
    Real max_torque;    // Peak motor torque [N·m]
    Real max_speed;     // Saturation speed [rad/s]
    Real friction_coef; // Viscous friction coefficient [N·m·s/rad]
    Vec3 spin_axis_b;   // Unit spin axis in body frame (normalise before use)
 
    Real omega_w = 0.0; // Current wheel speed [rad/s]
 
    // Step the wheel and return the reaction torque on the spacecraft body.
    // dt should match your outer integrator timestep.
    ActuatorOutput update(Real tau_cmd, Real dt);
 
    Vec3 momentum()     const;
    bool is_saturated() const;
};
 
// ─────────────────────────────────────────────────────────────────────────────
// REACTION WHEEL ARRAY  (N wheels, arbitrary geometry)
//
// Typical configs: 3-wheel orthogonal, 4-wheel pyramid, 6-wheel redundant.
// A pseudo-inverse distributes a desired 3-vector body torque across all wheels.
//
// NOTE: Template — fully defined here; cannot be split into a .cpp.
// ─────────────────────────────────────────────────────────────────────────────
template<std::size_t N>
struct ReactionWheelArray {
    std::array<ReactionWheel, N> wheels;
 
    // 3×N matrix whose columns are the spin axes of each wheel
    MatX axis_matrix() const
    {
        MatX A(3, static_cast<int>(N));
        for (std::size_t i = 0; i < N; ++i)
            A.col(i) = wheels[i].spin_axis_b;
        return A;
    }
 
    // Distribute desired body torque via Moore-Penrose pseudo-inverse.
    // Achieved torque may differ if individual wheels saturate.
    ActuatorOutput update(const Vec3& tau_cmd_b, Real dt)
    {
        const MatX A     = axis_matrix();
        const MatX Apinv = A.completeOrthogonalDecomposition().pseudoInverse();
        const VecX tau_w = Apinv * tau_cmd_b;
 
        ActuatorOutput total;
        for (std::size_t i = 0; i < N; ++i)
            total += wheels[i].update(tau_w(static_cast<int>(i)), dt);
        return total;
    }
 
    Vec3 total_momentum() const
    {
        Vec3 h = Vec3::Zero();
        for (const auto& w : wheels) h += w.momentum();
        return h;
    }
 
    bool any_saturated() const
    {
        for (const auto& w : wheels)
            if (w.is_saturated()) return true;
        return false;
    }
};
 
// ─────────────────────────────────────────────────────────────────────────────
// THRUSTER
//
// A single throttleable thruster at a fixed body-frame offset from the CoM.
// Input: duty_cycle ∈ [0, 1]  (0 = off, 1 = full thrust; bang-bang → 0 or 1)
//
// Force  = thrust · direction_b
// Torque = position_b × Force
//
// Force is accumulated in ActuatorOutput.force_b for coupled orbit/attitude.
// ─────────────────────────────────────────────────────────────────────────────
struct Thruster {
    Real max_thrust;  // Peak thrust [N]
    Real min_impulse; // Minimum impulse bit [N·s]; 0 = no check
    Vec3 position_b;  // Position in body frame from CoM [m]
    Vec3 direction_b; // Thrust direction unit vector (body frame)
 
    ActuatorOutput fire(Real duty_cycle) const;
};
 
// Thruster pair: two opposing thrusters that produce pure torque about one axis
// with no net force. Pass a signed torque command; the correct side fires.
struct ThrusterPair {
    Thruster pos_thruster; // fires when torque_cmd > 0
    Thruster neg_thruster; // fires when torque_cmd < 0
 
    ActuatorOutput fire(Real torque_cmd) const;
};
 
// ─────────────────────────────────────────────────────────────────────────────
// MAGNETORQUER
//
// Generates a magnetic dipole moment m [A·m²] along a fixed body axis.
// Torque: τ = m × B_body
//
// One axis is always uncontrollable (the component parallel to B).
// Full 3-axis authority requires orbital variation of B over time.
//
// Input: dipole_cmd ∈ [−1, 1] (fraction of max_dipole)
//        B_body               [T] — body-frame magnetic field
// ─────────────────────────────────────────────────────────────────────────────
struct Magnetorquer {
    Real max_dipole; // Maximum dipole moment [A·m²]
    Vec3 axis_b;     // Unit rod axis in body frame
 
    ActuatorOutput apply(Real dipole_cmd, const Vec3& B_body) const;
 
    // Returns how aligned this rod axis is with B (0 = perpendicular, 1 = parallel).
    // Values near 1 mean this rod has near-zero torque authority right now.
    Real alignment_with_field(const Vec3& B_body) const;
};
 
// ─────────────────────────────────────────────────────────────────────────────
// MAGNETORQUER ASSEMBLY  (3 orthogonal rods — standard CubeSat layout)
// ─────────────────────────────────────────────────────────────────────────────
struct MagnetorquerAssembly {
    std::array<Magnetorquer, 3> mtqs;
 
    // Torque-command interface — solves τ = m × B for the minimum-norm dipole:
    //   m = (B × τ_desired) / |B|²
    // The component of τ parallel to B is unachievable and is silently dropped.
    ActuatorOutput apply_torque_cmd(const Vec3& tau_cmd_b,
                                    const Vec3& B_body) const;
 
    // Raw dipole-command interface — used for B-dot detumbling:
    //   dipole_cmd_b = −k_bdot · B_dot_body   (pass the 3-vector directly)
    ActuatorOutput apply_dipole_cmd(const Vec3& dipole_cmd_b,
                                    const Vec3& B_body) const;
 
    // Convenience factory: standard X/Y/Z aligned rods, all same max dipole
    static MagnetorquerAssembly make_xyz(Real max_dipole_Am2);
};
 
// ─────────────────────────────────────────────────────────────────────────────
// ACTUATOR SUITE  —  all subsystems in one place
// NOTE: Template — fully defined here; cannot be split into a .cpp.
// ─────────────────────────────────────────────────────────────────────────────
template<std::size_t N_RWA = 4, std::size_t N_THR_PAIRS = 3>
struct ActuatorSuite {
    ReactionWheelArray<N_RWA>             rwa;
    MagnetorquerAssembly                  mtq;
    std::array<ThrusterPair, N_THR_PAIRS> thruster_pairs;
};
 
// ─────────────────────────────────────────────────────────────────────────────
// ACTUATOR COMMANDS  —  inputs delivered every integration step
// NOTE: Template — fully defined here; cannot be split into a .cpp.
// ─────────────────────────────────────────────────────────────────────────────
template<std::size_t N_RWA = 4, std::size_t N_THR_PAIRS = 3>
struct ActuatorCommands {
    Vec3 tau_rwa_cmd = Vec3::Zero(); // Desired body torque from RWA [N·m]
    Vec3 tau_mtq_cmd = Vec3::Zero(); // Desired body torque from MTQs [N·m]
    std::array<Real, N_THR_PAIRS> thr_cmds = {}; // Signed torque per pair [N·m]
};
 
// ─────────────────────────────────────────────────────────────────────────────
// computeActuatorAttitudeDynamics()
//
// Returns the AttitudeStateDot contribution from all actuators. Add its
// omega_dot field to the omega_dot returned by computeTotalAttitudeDynamics()
// before integrating. The q_dot field is always zero (actuators do not
// directly affect kinematics).
//
// The split is exact:
//   I dω/dt = (τ_dist − ω×Iω)    ← computeTotalAttitudeDynamics
//           + (τ_act  − ω×h_rw)   ← this function
//
// Parameters:
//   attitude_state      Current attitude state (q, omega)
//   r_eci               Spacecraft ECI position [m]  ← eci_state.position
//   inertia_tensor      Inertia tensor in body frame [kg·m²]
//   inertia_tensor_inv  Inverse inertia tensor
//   suite               Actuator suite — wheel speeds are integrated in-place
//   cmds                Commanded inputs for this step
//   t_sec               Simulation time since J2000 [s] (for B-field model)
//   dt                  Integration timestep [s] (for wheel dynamics)
//
// NOTE: Template — fully defined here; cannot be split into a .cpp.
// ─────────────────────────────────────────────────────────────────────────────
template<std::size_t N_RWA, std::size_t N_THR_PAIRS>
AttitudeStateDot computeActuatorAttitudeDynamics(
    const AttitudeState&                         attitude_state,
    const Vec3&                                  r_eci,
    const Mat3&                                  inertia_tensor,
    const Mat3&                                  inertia_tensor_inv,
    ActuatorSuite<N_RWA, N_THR_PAIRS>&           suite,
    const ActuatorCommands<N_RWA, N_THR_PAIRS>&  cmds,
    Real                                         t_sec,
    Real                                         dt)
{
    // Rotate B from ECI into body frame using the attitude quaternion
    const Vec3 B_eci  = dipole_field_eci(r_eci, t_sec);
    const Vec3 B_body = attitude_state.q.rotate(B_eci);
 
    // Accumulate torques from all subsystems
    ActuatorOutput total;
    total += suite.rwa.update(cmds.tau_rwa_cmd, dt);
    total += suite.mtq.apply_torque_cmd(cmds.tau_mtq_cmd, B_body);
    for (std::size_t i = 0; i < N_THR_PAIRS; ++i)
        total += suite.thruster_pairs[i].fire(cmds.thr_cmds[i]);
 
    // Gyroscopic term from stored wheel momentum (separate from ω×Iω in dynamics)
    const Vec3& omega     = attitude_state.omega;
    const Vec3  h_rw      = suite.rwa.total_momentum();
    const Vec3  omega_dot = inertia_tensor_inv
                          * (total.torque_b - omega.cross(h_rw));
 
    return AttitudeStateDot{
        quaternion_kinematics(attitude_state.q, omega),
        omega_dot
    };
}

struct QuaternionPDController {
    Real k_p;  // Proportional gain [N·m]
    Real k_d;  // Derivative gain   [N·m·s/rad]
 
    Vec3 compute_torque(const AttitudeState& state,
                        const Quaternion&    q_desired) const;
 
    // Error quaternion (q_desired* ⊗ q_current), sign-corrected for short arc.
    static Quaternion error_quaternion(const AttitudeState& state,
                                       const Quaternion&    q_desired);
};

Quaternion make_desired_quaternion(const Vec3& primary_axis_b,
                                   const Vec3& primary_target,
                                   const Vec3& secondary_axis_b,
                                   const Vec3& secondary_target);


 
} // namespace orb