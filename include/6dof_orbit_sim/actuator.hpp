#pragma once
#include "types.hpp"
#include <array>
#include <cmath>
#include <stdexcept>

// ─────────────────────────────────────────────────────────────────────────────
// Magnetic field model  —  IGRF tilted co-rotating dipole (ECI output)
//
// Uses a simple titled-dipole model that rotates with the Earth.
// Accuracy: ~few hundred nT at LEO — good enough for control design and
// hardware-in-the-loop, but replace with full IGRF for precision work.
//
//   r_eci  : spacecraft position in ECI [m]
//   t_sec  : seconds since J2000.0 epoch (TT or UTC both fine for this model)
//   returns: magnetic field vector in ECI [T]
// ─────────────────────────────────────────────────────────────────────────────
inline Vec3 dipole_field_eci(const Vec3& r_eci, Real t_sec)
{
    // --- Earth dipole axis in ECI ---
    // Rotates with Earth starting from its J2000 longitude reference
    const Real theta = constants::EARTH_ROT_RATE * t_sec + constants::MAG_DIPOLE_LON_J2000;
    const Real sin_tilt = std::sin(constants::MAG_DIPOLE_TILT);
    const Real cos_tilt = std::cos(constants::MAG_DIPOLE_TILT);

    // Unit dipole axis in ECI (co-rotating with Earth)
    const Vec3 m_hat(
        sin_tilt * std::cos(theta),
        sin_tilt * std::sin(theta),
        cos_tilt
    );

    // --- Dipole field formula: B = (μ₀/4π) * m_E / r³ * [3(m̂·r̂)r̂ − m̂] ---
    const Real r = r_eci.norm();
    if (r < 1.0) throw std::domain_error("dipole_field_eci: |r_eci| is near zero");

    const Vec3 r_hat = r_eci / r;
    const Real scale = constants::MU0_OVER_4PI * constants::EARTH_MAG_MOM / (r * r * r);
    const Real mdotr = m_hat.dot(r_hat);

    return scale * (3.0 * mdotr * r_hat - m_hat);
}

// ─────────────────────────────────────────────────────────────────────────────
// Common output structure — torque (and optionally force) produced in one
// integration step.  Accumulate across all actuators before applying.
// ─────────────────────────────────────────────────────────────────────────────
struct ActuatorOutput {
    Vec3 torque_b = Vec3::Zero(); // Body-frame torque [N·m]
    Vec3 force_b  = Vec3::Zero(); // Body-frame force  [N]  (thrusters only)
};

// ─────────────────────────────────────────────────────────────────────────────
// REACTION WHEEL
//
// A single momentum wheel spinning about a fixed body-frame axis.
//
// State:  omega_w  — wheel angular velocity [rad/s]
// Input:  tau_cmd  — commanded motor torque [N·m]  (positive = accelerate wheel)
//
// Physics:
//   I_w * dω_w/dt  =  τ_cmd − τ_friction        (wheel EOM)
//   τ_on_body      = −τ_cmd                       (Newton 3rd law)
//
// The angular momentum stored in the wheel array must be tracked separately
// and included in Euler's equation as:
//   I_sc * dω/dt = τ_ext − ω × (I_sc * ω + h_rw_total)
// where h_rw_total = Σ I_w,i * ω_w,i * â_i
// ─────────────────────────────────────────────────────────────────────────────
struct ReactionWheel {
    // ---- Configuration (set once) ----
    Real          inertia;       // Wheel moment of inertia [kg·m²]
    Real          max_torque;    // Peak motor torque [N·m]
    Real          max_speed;     // Saturation speed [rad/s]
    Real          friction_coef; // Viscous friction coefficient [N·m·s/rad]
    Vec3 spin_axis_b;   // Unit spin axis in body frame

    // ---- State (integrated each step) ----
    Real omega_w = 0.0;          // Current wheel speed [rad/s]

    // Clamp and apply commanded torque; returns torque exerted ON the spacecraft body.
    // Call update() every integration step, then add output.torque_b to tau_total.
    ActuatorOutput update(Real tau_cmd, Real dt)
    {
        // Clamp commanded torque
        const Real tau = std::clamp(tau_cmd, -max_torque, max_torque);

        // Viscous friction opposes wheel motion
        const Real tau_friction = friction_coef * omega_w;
        const Real tau_net_wheel = tau - tau_friction;

        // Integrate wheel speed (simple Euler — match to your integrator)
        omega_w += (tau_net_wheel / inertia) * dt;

        // Hard saturation: if the wheel hits max speed it can't accelerate further
        // (real wheels would trigger desaturation logic upstream)
        omega_w = std::clamp(omega_w, -max_speed, max_speed);

        ActuatorOutput out;
        // Reaction torque on body is equal and opposite to motor torque on wheel
        out.torque_b = -tau * spin_axis_b;
        return out;
    }

    // Angular momentum currently stored in this wheel [N·m·s]
    Vec3 momentum() const { return inertia * omega_w * spin_axis_b; }

    // True when the wheel is saturated (body torque authority = 0 in this axis)
    bool is_saturated() const { return std::abs(omega_w) >= max_speed * 0.999; }
};

// ─────────────────────────────────────────────────────────────────────────────
// REACTION WHEEL ARRAY  (convenience wrapper for N wheels)
//
// Typical configurations:
//   3-wheel orthogonal, 4-wheel pyramid, 6-wheel (pairs per axis)
//
// Usage:
//   ReactionWheelArray<4> rwa;
//   rwa.wheels[i] = { ... };          // configure
//   auto out = rwa.update(tau_cmd_vec, dt);   // 3-vector command → outputs
// ─────────────────────────────────────────────────────────────────────────────
template<std::size_t N>
struct ReactionWheelArray {
    std::array<ReactionWheel, N> wheels;

    // Configuration matrix A (3×N): each column is the spin axis of wheel i.
    // Build once from the configured wheel axes.
    MatX axis_matrix() const
    {
        MatX A(3, N);
        for (std::size_t i = 0; i < N; ++i)
            A.col(i) = wheels[i].spin_axis_b;
        return A;
    }

    // Distribute a desired body torque vector across N wheels via
    // pseudo-inverse allocation.  Returns the achieved body torque
    // (may differ if wheels saturate).
    //
    // tau_cmd_b : desired 3-vector body torque [N·m]
    // dt        : timestep [s]
    ActuatorOutput update(const Vec3& tau_cmd_b, Real dt)
    {
        // Moore-Penrose pseudo-inverse: τ_w = A† τ_body
        const MatX A     = axis_matrix();
        const MatX Apinv = A.completeOrthogonalDecomposition().pseudoInverse();
        const VecX tau_w = Apinv * tau_cmd_b;

        ActuatorOutput total;
        for (std::size_t i = 0; i < N; ++i) {
            auto out = wheels[i].update(tau_w(i), dt);
            total.torque_b += out.torque_b;
        }
        return total;
    }

    // Total angular momentum stored across all wheels [N·m·s]
    Vec3 total_momentum() const
    {
        Vec3 h = Vec3::Zero();
        for (const auto& w : wheels) h += w.momentum();
        return h;
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// THRUSTER
//
// A single on/off (or throttleable) thruster at a fixed body-frame location.
//
// Input: duty_cycle ∈ [0, 1]  (0 = off, 1 = full thrust)
//        For bang-bang attitude control, pass 0 or 1.
//        For continuous approximation, pass a fractional value.
//
// Produced torque is r × F where r is measured from the spacecraft CoM.
// Produced force is also accumulated — relevant for coupled attitude/orbit.
// ─────────────────────────────────────────────────────────────────────────────
struct Thruster {
    // ---- Configuration ----
    Real          max_thrust;   // Peak thrust [N]
    Real          min_impulse;  // Minimum impulse bit [N·s] — 0 disables check
    Vec3 position_b;   // Thruster position in body frame [m] from CoM
    Vec3 direction_b;  // Thrust direction unit vector (body frame)
                                  // (pointing FROM nozzle, i.e. direction of force)

    ActuatorOutput fire(Real duty_cycle) const
    {
        const Real dc = std::clamp(duty_cycle, 0.0, 1.0);
        const Real thrust = max_thrust * dc;
        const Vec3 F = thrust * direction_b;

        ActuatorOutput out;
        out.force_b  = F;
        out.torque_b = position_b.cross(F);  // τ = r × F
        return out;
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// THRUSTER PAIR  (couple for pure torque, no net force)
//
// Two thrusters fired simultaneously, producing torque about a body axis.
// ─────────────────────────────────────────────────────────────────────────────
struct ThrusterPair {
    Thruster pos_thruster;  // fires when torque_cmd > 0
    Thruster neg_thruster;  // fires when torque_cmd < 0

    // torque_cmd : signed desired torque magnitude [N·m] along the pair's axis
    ActuatorOutput fire(Real torque_cmd) const
    {
        if (torque_cmd > 0.0) return pos_thruster.fire(1.0);
        if (torque_cmd < 0.0) return neg_thruster.fire(1.0);
        return {};  // zero command → no fire
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// MAGNETORQUER
//
// Generates a magnetic dipole moment m [A·m²] along a fixed body axis.
// Torque produced: τ = m × B_body
//
// Note: magnetorquers cannot produce torque parallel to B — there is always
// one uncontrollable axis at any instant.  Three-axis control requires that
// the orbit provides sufficient variation in B direction over time.
//
// Input: dipole_cmd ∈ [-1, 1]  (fraction of max dipole)
//        B_body               : magnetic field in body frame [T]
// ─────────────────────────────────────────────────────────────────────────────
struct Magnetorquer {
    // ---- Configuration ----
    Real          max_dipole;  // Maximum magnetic dipole moment [A·m²]
    Vec3 axis_b;      // Unit axis of dipole in body frame

    ActuatorOutput apply(Real dipole_cmd, const Vec3& B_body) const
    {
        const Real m_mag = max_dipole * std::clamp(dipole_cmd, -1.0, 1.0);
        const Vec3 m_vec = m_mag * axis_b;

        ActuatorOutput out;
        out.torque_b = m_vec.cross(B_body);  // τ = m × B
        return out;
    }

    // Null space check: returns how aligned the commanded axis is with B.
    // Values near 1.0 mean this torquer is nearly ineffective right now.
    double alignment_with_field(const Vec3& B_body) const
    {
        const double Bnorm = B_body.norm();
        if (Bnorm < 1e-15) return 0.0;
        return std::abs(axis_b.dot(B_body / Bnorm));
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// MAGNETORQUER ASSEMBLY  (3 orthogonal rods — standard cubesat config)
//
// Provides a B-dot or cross-product control law interface.
// ─────────────────────────────────────────────────────────────────────────────
struct MagnetorquerAssembly {
    std::array<Magnetorquer, 3> mtqs;

    // Distribute a desired body torque into dipole commands via
    // τ_desired = m × B  →  solve for m via cross-product inverse.
    //
    // The exact solution: m = (B × τ_desired) / |B|²
    // This is the minimum-norm dipole that produces τ_desired.
    // The component of τ along B cannot be produced — it is silently dropped.
    //
    // Returns the body torque that will actually be produced.
    ActuatorOutput apply_torque_cmd(const Vec3& tau_cmd_b,
                                    const Vec3& B_body) const
    {
        const double Bsq = B_body.squaredNorm();
        if (Bsq < 1e-20) return {};  // no field → no torque

        // Minimum-norm dipole solution
        const Eigen::Vector3d m_desired = B_body.cross(tau_cmd_b) / Bsq;

        ActuatorOutput total;
        for (int i = 0; i < 3; ++i) {
            // Project desired dipole onto each rod axis
            const double cmd_Am2 = m_desired.dot(mtqs[i].axis_b);
            const double cmd_norm = cmd_Am2 / mtqs[i].max_dipole;  // normalize
            total.torque_b += mtqs[i].apply(cmd_norm, B_body).torque_b;
        }
        return total;
    }

    // Raw dipole command — useful for B-dot detumbling:
    //   dipole_cmd_vec = -k * B_dot  (each component normalized to [-1,1])
    ActuatorOutput apply_dipole_cmd(const Eigen::Vector3d& dipole_cmd_b,
                                    const Eigen::Vector3d& B_body) const
    {
        ActuatorOutput total;
        for (int i = 0; i < 3; ++i) {
            const double cmd = dipole_cmd_b.dot(mtqs[i].axis_b) / mtqs[i].max_dipole;
            total.torque_b += mtqs[i].apply(cmd, B_body).torque_b;
        }
        return total;
    }

    // Build a standard 3-axis assembly with body-aligned rods (X, Y, Z)
    static MagnetorquerAssembly make_xyz(double max_dipole_Am2)
    {
        MagnetorquerAssembly asm_;
        asm_.mtqs[0] = { max_dipole_Am2, Eigen::Vector3d::UnitX() };
        asm_.mtqs[1] = { max_dipole_Am2, Eigen::Vector3d::UnitY() };
        asm_.mtqs[2] = { max_dipole_Am2, Eigen::Vector3d::UnitZ() };
        return asm_;
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// FULL ACTUATOR SUITE  (aggregate for a single spacecraft)
//
// Example usage in your propagator's state-derivative function:
//
//   // 1. Get B in body frame
//   Eigen::Vector3d B_eci = dipole_field_eci(r_eci, t);
//   Eigen::Vector3d B_b   = R_eci_to_body * B_eci;
//
//   // 2. Compute commanded actuator outputs
//   auto rwa_out = suite.rwa.update(tau_rwa_cmd, dt);
//   auto mtq_out = suite.mtq.apply_torque_cmd(tau_mtq_cmd, B_b);
//   auto thr_out = suite.thruster_pairs[0].fire(torque_cmd_x);
//
//   // 3. Sum total torque (and force if needed)
//   Eigen::Vector3d tau_total = rwa_out.torque_b
//                             + mtq_out.torque_b
//                             + thr_out.torque_b;
//
//   // 4. Euler's equation (include stored RWA momentum)
//   //    I_sc * dw/dt = tau_total - w x (I_sc*w + rwa.total_momentum())
//   Eigen::Vector3d h_total = I_sc * omega + suite.rwa.total_momentum();
//   Eigen::Vector3d domega  = I_sc_inv * (tau_total - omega.cross(h_total));
// ─────────────────────────────────────────────────────────────────────────────
template<std::size_t N_RWA = 4, std::size_t N_THR_PAIRS = 3>
struct ActuatorSuite {
    ReactionWheelArray<N_RWA>               rwa;
    MagnetorquerAssembly                    mtq;
    std::array<ThrusterPair, N_THR_PAIRS>   thruster_pairs;

    // Helper: compute total body torque from all active actuators.
    // Pass zero vectors for actuators not in use this step.
    ActuatorOutput compute_total(
        const Eigen::Vector3d& tau_rwa_cmd,          // desired RWA body torque
        const Eigen::Vector3d& tau_mtq_cmd,          // desired MTQ body torque
        const Eigen::Vector3d& B_body,               // magnetic field in body frame [T]
        const std::array<double, N_THR_PAIRS>& thr_cmds,  // signed torque cmd per pair
        double dt)
    {
        ActuatorOutput total;

        total.torque_b += rwa.update(tau_rwa_cmd, dt).torque_b;
        total.torque_b += mtq.apply_torque_cmd(tau_mtq_cmd, B_body).torque_b;

        for (std::size_t i = 0; i < N_THR_PAIRS; ++i) {
            auto out = thruster_pairs[i].fire(thr_cmds[i]);
            total.torque_b += out.torque_b;
            total.force_b  += out.force_b;
        }

        return total;
    }
};