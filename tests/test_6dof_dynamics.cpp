#include "6dof_orbit_sim/6dof_dynamics.hpp"
#include <iostream>
#include <cmath>
#include <string>

using namespace orb;

// ─────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────

static int g_passed = 0;
static int g_failed = 0;

static void check(bool condition, const std::string& label) {
    if (condition) {
        std::cout << "PASS  " << label << "\n";
        ++g_passed;
    } else {
        std::cout << "FAIL  " << label << "\n";
        ++g_failed;
    }
}

static Quaternion identityQuat() {
    return {1.0, 0.0, 0.0, 0.0};
}

static Mat3 diagI(Real Ix, Real Iy, Real Iz) {
    Mat3 I = Mat3::Zero();
    I(0,0)=Ix; I(1,1)=Iy; I(2,2)=Iz;
    return I;
}

static Mat3 diagIinv(Real Ix, Real Iy, Real Iz) {
    return diagI(1.0/Ix, 1.0/Iy, 1.0/Iz);
}

// ─────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────

static void test6DoF_matchesSeparateModels() {
    SixDoFState state;

    state.eci.position << 6.778e6, 0.0, 0.0;
    state.eci.velocity << 0.0, 7.8e3, 0.0;

    state.attitude.q = identityQuat();
    state.attitude.omega << 0.0, 0.0, 0.0;

    Mat3 I    = diagI(10.0, 20.0, 30.0);
    Mat3 Iinv = diagIinv(10.0, 20.0, 30.0);

    ForceModelConfig force_cfg;
    force_cfg.useJ2   = true;
    force_cfg.useDrag = true;
    force_cfg.useSRP  = false;
    force_cfg.Cd      = 2.2;
    force_cfg.A       = 1.0;
    force_cfg.m       = 100.0;

    TorqueModelConfig torque_cfg;
    torque_cfg.useGravityGradient = true;
    torque_cfg.useDrag = true;
    torque_cfg.useSRP  = false;
    torque_cfg.center_of_mass = Vec3(0.1, 0.0, 0.0);
    torque_cfg.Cd   = 2.2;
    torque_cfg.area = 1.0;

    SixDoFStateDot sixdof = sixDoFDynamics(state, I, Iinv, force_cfg, torque_cfg);

    ECIStateDot eci_ref = computeTotalTranslationalDynamics(state.eci, force_cfg);
    AttitudeStateDot att_ref = computeTotalAttitudeDynamics(state.attitude, state.eci,
                                                           I, Iinv, torque_cfg);

    bool pass =
        nearlyEqual((sixdof.eci_dot.velocity - eci_ref.velocity).norm(), 0.0, 1e-12) &&
        nearlyEqual((sixdof.eci_dot.acceleration - eci_ref.acceleration).norm(), 0.0, 1e-12) &&
        nearlyEqual((sixdof.attitude_dot.omega_dot - att_ref.omega_dot).norm(), 0.0, 1e-12);

    check(pass,
          "sixDoFDynamics: matches translational + attitude models");
}

static void test6DoF_zeroOptionalModelsStillRuns() {
    SixDoFState state;

    state.eci.position << 7.0e6, 0.0, 0.0;
    state.eci.velocity << 0.0, 7.5e3, 0.0;

    state.attitude.q = identityQuat();
    state.attitude.omega << 0.01, 0.0, 0.0;

    Mat3 I    = diagI(5.0, 6.0, 7.0);
    Mat3 Iinv = diagIinv(5.0, 6.0, 7.0);

    ForceModelConfig force_cfg;
    force_cfg.useJ2   = false;
    force_cfg.useDrag = false;
    force_cfg.useSRP  = false;

    TorqueModelConfig torque_cfg;
    torque_cfg.useGravityGradient = false;
    torque_cfg.useDrag = false;
    torque_cfg.useSRP  = false;

    SixDoFStateDot result = sixDoFDynamics(state, I, Iinv, force_cfg, torque_cfg);

    check(std::isfinite(result.eci_dot.acceleration.norm()) &&
          std::isfinite(result.attitude_dot.omega_dot.norm()),
          "sixDoFDynamics: runs with all optional models disabled");
}

static void test6DoF_attitudeRespondsToTorque() {
    SixDoFState state;

    state.eci.position << 6.778e6, 0.0, 0.0;
    state.eci.velocity << 0.0, 7.8e3, 0.0;

    state.attitude.q = identityQuat();
    state.attitude.omega << 0.0, 0.0, 0.0;

    Mat3 I    = diagI(10.0, 20.0, 30.0);
    Mat3 Iinv = diagIinv(10.0, 20.0, 30.0);

    ForceModelConfig force_cfg;

    TorqueModelConfig torque_cfg;
    torque_cfg.useGravityGradient = false;
    torque_cfg.useDrag = true;
    torque_cfg.useSRP  = false;
    torque_cfg.center_of_mass = Vec3(0.1, 0.0, 0.0);
    torque_cfg.Cd   = 2.2;
    torque_cfg.area = 1.0;

    auto result = sixDoFDynamics(state, I, Iinv, force_cfg, torque_cfg);

    check(result.attitude_dot.omega_dot.norm() > 0.0,
          "sixDoFDynamics: nonzero torque produces angular acceleration");
}

// ─────────────────────────────────────────────
// Main
// ─────────────────────────────────────────────

int main() {
    std::cout << "\n=== 6DoF Dynamics ===\n";

    test6DoF_matchesSeparateModels();
    test6DoF_zeroOptionalModelsStillRuns();
    test6DoF_attitudeRespondsToTorque();

    std::cout << "\n────────────────────────────\n";
    std::cout << "Results: " << g_passed << " passed, "
              << g_failed << " failed\n\n";

    return g_failed > 0 ? 1 : 0;
}