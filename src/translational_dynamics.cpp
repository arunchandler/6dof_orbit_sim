#include "6dof_orbit_sim/translational_dynamics.hpp"
#include <cmath>
#include <stdexcept>

namespace orb {

ECIStateDot accTwoBody(const ECIState& state, Real mu) {
    const Vec3 r_vec = state.position;
    const Real r = r_vec.norm();
    if (r == 0.0)
        throw std::invalid_argument("Position vector cannot be zero for two-body force calculation.");

    ECIStateDot derivative;
    derivative.velocity = state.velocity; // dr/dt = v
    derivative.acceleration = -mu * r_vec / (r * r * r);
    return derivative;
}

ECIStateDot accJ2(const ECIState& state, Real mu, Real J2, Real R) {
    const Vec3 r_vec = state.position;
    const Real r = r_vec.norm();
    if (r == 0.0)
        throw std::invalid_argument("Position vector cannot be zero for J2 force calculation.");

    const Real z2 = r_vec.z() * r_vec.z();
    const Real r2 = r * r;
    const Real factor = (3.0/2.0) * J2 * mu * R * R / (r2 * r2 * r);

    ECIStateDot derivative;
    derivative.velocity = state.velocity;
    derivative.acceleration << factor * r_vec.x() * (5.0*z2/r2 - 1.0),
                               factor * r_vec.y() * (5.0*z2/r2 - 1.0),
                               factor * r_vec.z() * (5.0*z2/r2 - 3.0);
    return derivative;
}

// Simple exponential atmosphere model for drag force
ECIStateDot accDrag(const ECIState& state, Real mu, Real Cd, Real A, Real m) {
    const Vec3 r_vec = state.position;
    const Vec3 v_vec = state.velocity;

    const Real r = r_vec.norm();
    const Real v = v_vec.norm();
    if (r == 0.0)
        throw std::invalid_argument("Position vector cannot be zero for drag force calculation.");

    // Simple exponential atmosphere model
    // const Real rho0 = 1.225; // kg/m³ at sea level
    // const Real H = 8500000.0;    // Scale height in meters
    // const Real rho = rho0 * std::exp(-r / H);
    const Real rho = 1e-10; // Placeholder - replace with actual atmosphere model

    ECIStateDot derivative;
    derivative.velocity = state.velocity;
    derivative.acceleration = -0.5 * Cd * A / m * rho * v * v_vec;
    return derivative;
}

ECIStateDot accSRP(const ECIState& state, Real mu, Real Cr, Real A, Real m) {
    const Vec3 r_vec = state.position;

    const Real r = r_vec.norm();
    if (r == 0.0)
        throw std::invalid_argument("Position vector cannot be zero for SRP force calculation.");

    const Real L = 3.846e26; // Solar luminosity in watts
    const Real P = L / (4.0 * M_PI * constants::AU * constants::AU * constants::C_LIGHT); // Solar radiation pressure at distance AU for simplicity & ignoring Earth's shadow

    ECIStateDot derivative;
    derivative.velocity = state.velocity;
    derivative.acceleration = -P * Cr * A / m * r_vec.normalized();
    return derivative;
}

ECIStateDot accTotal(const ECIState& state, const ForceModelConfig& config) {
    ECIStateDot totalAcc = accTwoBody(state, constants::MU_EARTH);
    if (config.useJ2)
        totalAcc.acceleration += accJ2(state, constants::MU_EARTH, constants::J2, constants::R_EARTH).acceleration;
    if (config.useDrag)
        totalAcc.acceleration += accDrag(state, constants::MU_EARTH, config.Cd, config.A, config.m).acceleration;
    if (config.useSRP)
        totalAcc.acceleration += accSRP(state, constants::MU_EARTH, config.Cr, config.A, config.m).acceleration;
    return totalAcc;
}

OrbitalElements dynamicsGVE(const OrbitalElements& oe, const ForceModelConfig& config, Real mu) {

    // ── 1. Convert elements → ECI state ──────────────────────────────────────
    ECIState state   = elementsToECI(oe, mu);
    Vec3 r_vec   = state.position;
    Vec3 v_vec   = state.velocity;

    // ── 2. Isolate perturbing acceleration in ECI ─────────────────────────────
    // Total minus two-body leaves only the perturbations
    Vec3 a_total = accTotal(state, config).acceleration;
    Vec3 a_2b    = accTwoBody(state, mu).acceleration;
    Vec3 a_pert  = a_total - a_2b;

    // ── 3. Build RSW frame ────────────────────────────────────────────────────
    Vec3 R_hat = r_vec.normalized();                 // radial (away from Earth)
    Vec3 W_hat = r_vec.cross(v_vec).normalized();    // normal (along h)
    Vec3 S_hat = W_hat.cross(R_hat).normalized();    // along-track

    // ── 4. Project perturbation into RSW ──────────────────────────────────────
    Real f_R = a_pert.dot(R_hat);
    Real f_S = a_pert.dot(S_hat);
    Real f_W = a_pert.dot(W_hat);

    // ── 5. Derived orbital quantities ─────────────────────────────────────────
    const Real e   = oe.ecc;
    const Real nu  = oe.ta;
    const Real inc = oe.inc;
    const Real aop = oe.aop;
    const Real p   = oe.sma * (1.0 - e * e);               // semi-latus rectum
    const Real h   = std::sqrt(mu * p);                    // specific angular momentum
    const Real r   = p / (1.0 + e * std::cos(nu));         // orbital radius

    // Singularity guards
    const Real eps_e   = 1e-10;   // circular orbit guard
    const Real eps_i   = 1e-10;   // equatorial orbit guard
    const Real sin_i   = std::fabs(std::sin(inc)) < eps_i
                         ? eps_i : std::sin(inc);

    // ── 6. Gauss Variational Equations ────────────────────────────────────────
    OrbitalElements oe_dot;

    // ȧ = (2a²/h) [e sinν f_R + (p/r) f_S]
    oe_dot.sma = (2.0 * oe.sma * oe.sma / h)
                 * (e * std::sin(nu) * f_R + (p / r) * f_S);

    // ė = (1/h) [p sinν f_R + ((p+r)cosν + re) f_S]
    oe_dot.ecc = (1.0 / h)
                 * (p * std::sin(nu) * f_R
                 + ((p + r) * std::cos(nu) + r * e) * f_S);

    // i̇ = (r cos(ω+ν) / h) f_W
    oe_dot.inc = (r * std::cos(aop + nu) / h) * f_W;

    // Ω̇ = (r sin(ω+ν)) / (h sini) f_W
    oe_dot.raan = (r * std::sin(aop + nu)) / (h * sin_i) * f_W;

    // ω̇ = (1/he)[-p cosν f_R + (p+r) sinν f_S] - (r sin(ω+ν) cosi)/(h sini) f_W
    if (e < eps_e) {
        oe_dot.aop = 0.0;   // undefined for circular orbits
    } else {
        oe_dot.aop = (1.0 / (h * e))
                     * (-p * std::cos(nu) * f_R
                     +  (p + r) * std::sin(nu) * f_S)
                     - (r * std::sin(aop + nu) * std::cos(inc))
                     / (h * sin_i) * f_W;
    }

    // ν̇ = h/r² + (1/he)[p cosν f_R - (p+r) sinν f_S]
    oe_dot.ta = h / (r * r);
    if (e > eps_e) {
        oe_dot.ta += (1.0 / (h * e))
                     * (p * std::cos(nu) * f_R
                     - (p + r) * std::sin(nu) * f_S);
    }

    return oe_dot;
}

} // namespace orb