#include "6dof_orbit_sim/attitude_dynamics.hpp"
#include <iostream>
#include <cmath>
#include <string>

using namespace orb;

// ─────────────────────────────────────────────
//  Helpers
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

/// Identity quaternion (no rotation)
static Quaternion identityQuat() {
    return {1.0, 0.0, 0.0, 0.0};
}

/// Build a diagonal 3x3 inertia tensor from principal moments
static Mat3 diagI(Real Ix, Real Iy, Real Iz) {
    Mat3 I = Mat3::Zero();
    I(0, 0) = Ix;
    I(1, 1) = Iy;
    I(2, 2) = Iz;
    return I;
}

/// Inverse of a diagonal inertia tensor
static Mat3 diagIinv(Real Ix, Real Iy, Real Iz) {
    return diagI(1.0 / Ix, 1.0 / Iy, 1.0 / Iz);
}

// ─────────────────────────────────────────────
//  1. Gravity-Gradient Attitude Dynamics
// ─────────────────────────────────────────────

static void testGgrad_symmetricInertiaZeroTorque() {
    // Spherical inertia: r̂ × (I·r̂) = r̂ × (s·r̂) = 0 for any r̂
    AttitudeState att  = {identityQuat(), Vec3(0.0, 0.0, 0.0)};
    ECIState      eci;
    eci.position << 0.0, 0.0, 6.778e6;
    eci.velocity << 0.0, 0.0, 0.0;
    Mat3 I    = diagI(10.0, 10.0, 10.0);
    Mat3 Iinv = diagIinv(10.0, 10.0, 10.0);

    AttitudeStateDot result = computeGravGradAttitudeDynamics(att, eci, I, Iinv);

    check(nearlyEqual(result.omega_dot(0), 0.0, 1e-10) &&
          nearlyEqual(result.omega_dot(1), 0.0, 1e-10) &&
          nearlyEqual(result.omega_dot(2), 0.0, 1e-10),
          "ggrad: spherical inertia produces zero angular acceleration");
}

static void testGgrad_principalAxisAlignedZeroTorque() {
    // With diagonal I and r̂_body = ẑ: I·ẑ = Iz·ẑ, so ẑ × Iz·ẑ = 0
    AttitudeState att  = {identityQuat(), Vec3(0.0, 0.0, 0.0)};
    ECIState      eci;
    eci.position << 0.0, 0.0, 6.778e6;
    eci.velocity << 0.0, 0.0, 0.0;

    Mat3 I    = diagI(10.0, 20.0, 30.0);
    Mat3 Iinv = diagIinv(10.0, 20.0, 30.0);

    AttitudeStateDot result = computeGravGradAttitudeDynamics(att, eci, I, Iinv);

    check(nearlyEqual(result.omega_dot(0), 0.0, 1e-10) &&
          nearlyEqual(result.omega_dot(1), 0.0, 1e-10) &&
          nearlyEqual(result.omega_dot(2), 0.0, 1e-10),
          "ggrad: r_hat aligned with principal axis produces zero torque");
}

static void testGgrad_inverseR3Scaling() {
    AttitudeState att = {identityQuat(), Vec3(0,0,0)};
    Real r1 = 6.778e6;
    Real r2 = 2.0*r1;

    Real s = 1.0/std::sqrt(2.0);

    Mat3 I    = diagI(10,20,30);
    Mat3 Iinv = diagIinv(10,20,30);

    ECIState eci1, eci2;
    eci1.position << r1*s, r1*s, 0;
    eci2.position << r2*s, r2*s, 0;

    eci1.velocity.setZero();
    eci2.velocity.setZero();

    auto d1 = computeGravGradAttitudeDynamics(att, eci1, I, Iinv);
    auto d2 = computeGravGradAttitudeDynamics(att, eci2, I, Iinv);

    Real ratio = d1.omega_dot.norm() / d2.omega_dot.norm();

    check(nearlyEqual(ratio, 8.0, 1e-6),
          "ggrad: angular acceleration scales as 1/r^3");
}

static void testGgrad_qdotZeroForZeroOmega() {
    // q̇ = 0.5 q ⊗ [0, ω]; with ω = 0 → all components of q_dot are zero
    AttitudeState att = {identityQuat(), Vec3(0.0, 0.0, 0.0)};
    ECIState      eci;
    eci.position << 0.0, 0.0, 6.778e6;
    eci.velocity << 0.0, 0.0, 0.0;

    Mat3 I    = diagI(10.0, 20.0, 30.0);
    Mat3 Iinv = diagIinv(10.0, 20.0, 30.0);

    AttitudeStateDot result = computeGravGradAttitudeDynamics(att, eci, I, Iinv);

    check(nearlyEqual(result.q_dot.w, 0.0, 1e-12) &&
          nearlyEqual(result.q_dot.x, 0.0, 1e-12) &&
          nearlyEqual(result.q_dot.y, 0.0, 1e-12) &&
          nearlyEqual(result.q_dot.z, 0.0, 1e-12),
          "ggrad: q_dot is zero when omega is zero");
}

static void testGgrad_qdotNonzeroForSpinning() {
    // Any non-zero omega should produce a non-zero q_dot
    AttitudeState att = {identityQuat(), Vec3(0.0, 0.0, 0.1)};
    ECIState      eci;
    eci.position << 0.0, 0.0, 6.778e6;
    eci.velocity << 0.0, 0.0, 0.0;

    Mat3 I    = diagI(10.0, 20.0, 30.0);
    Mat3 Iinv = diagIinv(10.0, 20.0, 30.0);

    AttitudeStateDot result = computeGravGradAttitudeDynamics(att, eci, I, Iinv);

    Real qdot_mag = std::sqrt(result.q_dot.w * result.q_dot.w +
                              result.q_dot.x * result.q_dot.x +
                              result.q_dot.y * result.q_dot.y +
                              result.q_dot.z * result.q_dot.z);
    check(qdot_mag > 1e-12,
          "ggrad: q_dot is nonzero when omega is nonzero");
}

static void testGgrad_spinOnPrincipalAxisNoPrecession() {
    // ω along x̂, r̂ along x̂ → torque = 0 and ω × Iω = 0 → omega_dot = 0
    AttitudeState att = {identityQuat(), Vec3(0.1, 0.0, 0.0)};
    ECIState      eci;
    eci.position << 6.778e6, 0.0, 0.0;
    eci.velocity << 0.0, 0.0, 0.0;

    Mat3 I    = diagI(10.0, 20.0, 30.0);
    Mat3 Iinv = diagIinv(10.0, 20.0, 30.0);

    AttitudeStateDot result = computeGravGradAttitudeDynamics(att, eci, I, Iinv);

    check(nearlyEqual(result.omega_dot(0), 0.0, 1e-8) &&
          nearlyEqual(result.omega_dot(1), 0.0, 1e-8) &&
          nearlyEqual(result.omega_dot(2), 0.0, 1e-8),
          "ggrad: spin about principal axis aligned with r_hat produces no precession");
}

// ─────────────────────────────────────────────
//  2. Drag Attitude Dynamics
// ─────────────────────────────────────────────

static void testDrag_comAtOriginZeroTorque() {
    // Force applied at CoM → zero moment arm → zero torque
    AttitudeState att = {identityQuat(), Vec3(0.0, 0.0, 0.0)};
    ECIState      eci;
    eci.position << 6.778e6, 0.0, 0.0;
    eci.velocity << 0.0, 7.8e3, 0.0;

    Vec3 com(0.0, 0.0, 0.0);
    Mat3 I    = diagI(10.0, 20.0, 30.0);
    Mat3 Iinv = diagIinv(10.0, 20.0, 30.0);

    AttitudeStateDot result = computeDragAttitudeDynamics(att, eci, com, I, Iinv, 2.2, 1.0);

    check(nearlyEqual(result.omega_dot(0), 0.0, 1e-20) &&
          nearlyEqual(result.omega_dot(1), 0.0, 1e-20) &&
          nearlyEqual(result.omega_dot(2), 0.0, 1e-20),
          "drag: CoM at origin produces zero angular acceleration");
}

static void testDrag_torqueSignWithComOffset() {
    // v along +y, drag force along -y in body frame; CoM offset along +x
    // torque = (-com) × force_body = (-0.1, 0, 0) × (0, -F, 0) = (0, 0, +0.1F)
    // → omega_dot_z > 0
    AttitudeState att = {identityQuat(), Vec3(0.0, 0.0, 0.0)};
    ECIState      eci;
    eci.position << 6.778e6, 0.0, 0.0;
    eci.velocity << 0.0, 7.8e3, 0.0;

    Vec3 com(0.1, 0.0, 0.0);
    Mat3 I    = diagI(10.0, 20.0, 30.0);
    Mat3 Iinv = diagIinv(10.0, 20.0, 30.0);

    AttitudeStateDot result = computeDragAttitudeDynamics(att, eci, com, I, Iinv, 2.2, 1.0);

    check(result.omega_dot(2) > 0.0,
          "drag: CoM offset along +x with velocity along +y gives omega_dot_z > 0");
}

static void testDrag_scalesWithCd() {
    // Doubling Cd should double omega_dot
    AttitudeState att = {identityQuat(), Vec3(0.0, 0.0, 0.0)};
    ECIState      eci;
    eci.position << 6.778e6, 0.0, 0.0;
    eci.velocity << 0.0, 7.8e3, 0.0;

    Vec3 com(0.1, 0.0, 0.0);
    Mat3 I    = diagI(10.0, 20.0, 30.0);
    Mat3 Iinv = diagIinv(10.0, 20.0, 30.0);

    AttitudeStateDot d1 = computeDragAttitudeDynamics(att, eci, com, I, Iinv, 2.2, 1.0);
    AttitudeStateDot d2 = computeDragAttitudeDynamics(att, eci, com, I, Iinv, 4.4, 1.0);

    Real ratio = d2.omega_dot(2) / d1.omega_dot(2);
    check(nearlyEqual(ratio, 2.0, 1e-9),
          "drag: angular acceleration scales linearly with Cd");
}

static void testDrag_scalesWithArea() {
    // Doubling area should double omega_dot
    AttitudeState att = {identityQuat(), Vec3(0.0, 0.0, 0.0)};
    ECIState      eci;
    eci.position << 6.778e6, 0.0, 0.0;
    eci.velocity << 0.0, 7.8e3, 0.0;

    Vec3 com(0.1, 0.0, 0.0);
    Mat3 I    = diagI(10.0, 20.0, 30.0);
    Mat3 Iinv = diagIinv(10.0, 20.0, 30.0);

    AttitudeStateDot d1 = computeDragAttitudeDynamics(att, eci, com, I, Iinv, 2.2, 1.0);
    AttitudeStateDot d2 = computeDragAttitudeDynamics(att, eci, com, I, Iinv, 2.2, 2.0);

    Real ratio = d2.omega_dot(2) / d1.omega_dot(2);
    check(nearlyEqual(ratio, 2.0, 1e-9),
          "drag: angular acceleration scales linearly with area");
}

static void testDrag_zeroVelocityZeroTorque() {
    // Zero relative velocity → zero drag force → zero torque even with CoM offset
    AttitudeState att = {identityQuat(), Vec3(0.0, 0.0, 0.0)};
    ECIState      eci;
    eci.position << 6.778e6, 0.0, 0.0;
    eci.velocity << 0.0, 0.0, 0.0;

    Vec3 com(0.1, 0.0, 0.0);
    Mat3 I    = diagI(10.0, 20.0, 30.0);
    Mat3 Iinv = diagIinv(10.0, 20.0, 30.0);

    AttitudeStateDot result = computeDragAttitudeDynamics(att, eci, com, I, Iinv, 2.2, 1.0);

    check(nearlyEqual(result.omega_dot(0), 0.0, 1e-20) &&
          nearlyEqual(result.omega_dot(1), 0.0, 1e-20) &&
          nearlyEqual(result.omega_dot(2), 0.0, 1e-20),
          "drag: zero velocity produces zero angular acceleration");
}

static void testDrag_qdotZeroForZeroOmega() {
    AttitudeState att = {identityQuat(), Vec3(0.0, 0.0, 0.0)};
    ECIState      eci;
    eci.position << 6.778e6, 0.0, 0.0;
    eci.velocity << 0.0, 7.8e3, 0.0;

    Vec3 com(0.1, 0.0, 0.0);
    Mat3 I    = diagI(10.0, 20.0, 30.0);
    Mat3 Iinv = diagIinv(10.0, 20.0, 30.0);

    AttitudeStateDot result = computeDragAttitudeDynamics(att, eci, com, I, Iinv, 2.2, 1.0);

    check(nearlyEqual(result.q_dot.w, 0.0, 1e-12) &&
          nearlyEqual(result.q_dot.x, 0.0, 1e-12) &&
          nearlyEqual(result.q_dot.y, 0.0, 1e-12) &&
          nearlyEqual(result.q_dot.z, 0.0, 1e-12),
          "drag: q_dot is zero when omega is zero");
}

// ─────────────────────────────────────────────
//  3. SRP Attitude Dynamics
// ─────────────────────────────────────────────

static void testSRP_comAtOriginZeroTorque() {
    AttitudeState att = {identityQuat(), Vec3(0.0, 0.0, 0.0)};
    ECIState      eci;
    eci.position << 6.778e6, 0.0, 0.0;
    eci.velocity << 0.0, 7.8e3, 0.0;

    Vec3 sun_dir(1.0, 0.0, 0.0);
    Vec3 com(0.0, 0.0, 0.0);
    Mat3 I    = diagI(10.0, 20.0, 30.0);
    Mat3 Iinv = diagIinv(10.0, 20.0, 30.0);

    AttitudeStateDot result = computeSRPAttitudeDynamics(att, eci, sun_dir, com, I, Iinv, 1.5, 2.0, 4.56e-6);

    check(nearlyEqual(result.omega_dot(0), 0.0, 1e-20) &&
          nearlyEqual(result.omega_dot(1), 0.0, 1e-20) &&
          nearlyEqual(result.omega_dot(2), 0.0, 1e-20),
          "SRP: CoM at origin produces zero angular acceleration");
}

static void testSRP_torqueSignWithComOffset() {
    // sun along +x → force_body along -x; CoM offset along +y
    // torque = (-com) × force_body = (0, -0.05, 0) × (-F, 0, 0)
    //   k: (0)(0) - (-0.05)(-F) = -0.05F  → omega_dot_z < 0
    AttitudeState att = {identityQuat(), Vec3(0.0, 0.0, 0.0)};
    ECIState      eci;
    eci.position << 6.778e6, 0.0, 0.0;
    eci.velocity << 0.0, 7.8e3, 0.0;

    Vec3 sun_dir(1.0, 0.0, 0.0);
    Vec3 com(0.0, 0.05, 0.0);
    Mat3 I    = diagI(10.0, 20.0, 30.0);
    Mat3 Iinv = diagIinv(10.0, 20.0, 30.0);

    AttitudeStateDot result = computeSRPAttitudeDynamics(att, eci, sun_dir, com, I, Iinv, 1.5, 2.0, 4.56e-6);

    check(result.omega_dot(2) < 0.0,
          "SRP: CoM offset along +y with sun along +x gives omega_dot_z < 0");
}

static void testSRP_analyticalMagnitude() {
    // |τ_z| = P_srp * Cr * area * com_y  (from cross product above)
    // omega_dot_z = τ_z / Iz
    Real P     = 4.56e-6;
    Real Cr    = 1.5;
    Real area  = 2.0;
    Real com_y = 0.05;
    Real Iz    = 30.0;

    Real expected = -(P * Cr * area * com_y) / Iz;

    AttitudeState att = {identityQuat(), Vec3(0.0, 0.0, 0.0)};
    ECIState      eci;
    eci.position << 6.778e6, 0.0, 0.0;
    eci.velocity << 0.0, 7.8e3, 0.0;

    Vec3 sun_dir(1.0, 0.0, 0.0);
    Vec3 com(0.0, com_y, 0.0);
    Mat3 I    = diagI(10.0, 20.0, Iz);
    Mat3 Iinv = diagIinv(10.0, 20.0, Iz);

    AttitudeStateDot result = computeSRPAttitudeDynamics(att, eci, sun_dir, com, I, Iinv, Cr, area, P);

    check(nearlyEqual(result.omega_dot(2), expected, 1e-18),
          "SRP: omega_dot_z matches analytical value");
}

static void testSRP_scalesWithP() {
    AttitudeState att = {identityQuat(), Vec3(0.0, 0.0, 0.0)};
    ECIState      eci;
    eci.position << 6.778e6, 0.0, 0.0;
    eci.velocity << 0.0, 7.8e3, 0.0;

    Vec3 sun_dir(1.0, 0.0, 0.0);
    Vec3 com(0.0, 0.05, 0.0);
    Mat3 I    = diagI(10.0, 20.0, 30.0);
    Mat3 Iinv = diagIinv(10.0, 20.0, 30.0);

    AttitudeStateDot d1 = computeSRPAttitudeDynamics(att, eci, sun_dir, com, I, Iinv, 1.5, 2.0, 4.56e-6);
    AttitudeStateDot d2 = computeSRPAttitudeDynamics(att, eci, sun_dir, com, I, Iinv, 1.5, 2.0, 9.12e-6);

    Real ratio = d2.omega_dot(2) / d1.omega_dot(2);
    check(nearlyEqual(ratio, 2.0, 1e-9),
          "SRP: angular acceleration scales linearly with P_srp");
}

static void testSRP_scalesWithCr() {
    AttitudeState att = {identityQuat(), Vec3(0.0, 0.0, 0.0)};
    ECIState      eci;
    eci.position << 6.778e6, 0.0, 0.0;
    eci.velocity << 0.0, 7.8e3, 0.0;

    Vec3 sun_dir(1.0, 0.0, 0.0);
    Vec3 com(0.0, 0.05, 0.0);
    Mat3 I    = diagI(10.0, 20.0, 30.0);
    Mat3 Iinv = diagIinv(10.0, 20.0, 30.0);

    AttitudeStateDot d1 = computeSRPAttitudeDynamics(att, eci, sun_dir, com, I, Iinv, 1.5, 2.0, 4.56e-6);
    AttitudeStateDot d2 = computeSRPAttitudeDynamics(att, eci, sun_dir, com, I, Iinv, 3.0, 2.0, 4.56e-6);

    Real ratio = d2.omega_dot(2) / d1.omega_dot(2);
    check(nearlyEqual(ratio, 2.0, 1e-9),
          "SRP: angular acceleration scales linearly with Cr");
}

static void testSRP_scalesWithArea() {
    AttitudeState att = {identityQuat(), Vec3(0.0, 0.0, 0.0)};
    ECIState      eci;
    eci.position << 6.778e6, 0.0, 0.0;
    eci.velocity << 0.0, 7.8e3, 0.0;

    Vec3 sun_dir(1.0, 0.0, 0.0);
    Vec3 com(0.0, 0.05, 0.0);
    Mat3 I    = diagI(10.0, 20.0, 30.0);
    Mat3 Iinv = diagIinv(10.0, 20.0, 30.0);

    AttitudeStateDot d1 = computeSRPAttitudeDynamics(att, eci, sun_dir, com, I, Iinv, 1.5, 2.0, 4.56e-6);
    AttitudeStateDot d2 = computeSRPAttitudeDynamics(att, eci, sun_dir, com, I, Iinv, 1.5, 4.0, 4.56e-6);

    Real ratio = d2.omega_dot(2) / d1.omega_dot(2);
    check(nearlyEqual(ratio, 2.0, 1e-9),
          "SRP: angular acceleration scales linearly with area");
}

static void testSRP_reverseSunDirReversesTorque() {
    AttitudeState att = {identityQuat(), Vec3(0.0, 0.0, 0.0)};
    ECIState      eci;
    eci.position << 6.778e6, 0.0, 0.0;
    eci.velocity << 0.0, 7.8e3, 0.0;

    Vec3 sun_pos( 1.0, 0.0, 0.0);
    Vec3 sun_neg(-1.0, 0.0, 0.0);
    Vec3 com(0.0, 0.05, 0.0);
    Mat3 I    = diagI(10.0, 20.0, 30.0);
    Mat3 Iinv = diagIinv(10.0, 20.0, 30.0);

    AttitudeStateDot d1 = computeSRPAttitudeDynamics(att, eci, sun_pos, com, I, Iinv, 1.5, 2.0, 4.56e-6);
    AttitudeStateDot d2 = computeSRPAttitudeDynamics(att, eci, sun_neg, com, I, Iinv, 1.5, 2.0, 4.56e-6);

    check(nearlyEqual(d1.omega_dot(2), -d2.omega_dot(2), 1e-18),
          "SRP: reversing sun direction reverses the torque");
}

static void testSRP_qdotZeroForZeroOmega() {
    AttitudeState att = {identityQuat(), Vec3(0.0, 0.0, 0.0)};
    ECIState      eci;
    eci.position << 6.778e6, 0.0, 0.0;
    eci.velocity << 0.0, 7.8e3, 0.0;

    Vec3 sun_dir(1.0, 0.0, 0.0);
    Vec3 com(0.0, 0.05, 0.0);
    Mat3 I    = diagI(10.0, 20.0, 30.0);
    Mat3 Iinv = diagIinv(10.0, 20.0, 30.0);

    AttitudeStateDot result = computeSRPAttitudeDynamics(att, eci, sun_dir, com, I, Iinv, 1.5, 2.0, 4.56e-6);

    check(nearlyEqual(result.q_dot.w, 0.0, 1e-12) &&
          nearlyEqual(result.q_dot.x, 0.0, 1e-12) &&
          nearlyEqual(result.q_dot.y, 0.0, 1e-12) &&
          nearlyEqual(result.q_dot.z, 0.0, 1e-12),
          "SRP: q_dot is zero when omega is zero");
}

// ─────────────────────────────────────────────
//  4. Total Torque Model
// ─────────────────────────────────────────────

static void testTorqueTotal_matchesGravityGradientOnly() {
    AttitudeState att = {identityQuat(), Vec3(0.0, 0.0, 0.0)};
    ECIState eci;
    eci.position << 4.79e6, 4.79e6, 0.0;   // non-principal direction
    eci.velocity << 0.0, 0.0, 0.0;

    Mat3 I    = diagI(10.0, 20.0, 30.0);
    Mat3 Iinv = diagIinv(10.0, 20.0, 30.0);

    TorqueModelConfig cfg;
    cfg.useGravityGradient = true;
    cfg.useDrag = false;
    cfg.useSRP = false;

    auto ref = computeGravGradAttitudeDynamics(att, eci, I, Iinv, 0.0);
    auto tot = computeTotalAttitudeDynamics(att, eci, I, Iinv, cfg, 0.0);

    check(nearlyEqual((tot.omega_dot - ref.omega_dot).norm(), 0.0, 1e-15),
          "torqueTotal: gravity-gradient only matches ggradAttitudeDynamics");
}

static void testTorqueTotal_matchesDragOnly() {
    AttitudeState att = {identityQuat(), Vec3(0.0, 0.0, 0.0)};
    ECIState eci;
    eci.position << 6.778e6, 0.0, 0.0;
    eci.velocity << 0.0, 7.8e3, 0.0;

    Mat3 I    = diagI(10.0, 20.0, 30.0);
    Mat3 Iinv = diagIinv(10.0, 20.0, 30.0);

    TorqueModelConfig cfg;
    cfg.useGravityGradient = false;
    cfg.useDrag = true;
    cfg.useSRP = false;
    cfg.center_of_mass = Vec3(0.1, 0.0, 0.0);
    cfg.Cd = 2.2;
    cfg.area = 1.0;

    auto ref = computeDragAttitudeDynamics(att, eci, cfg.center_of_mass, I, Iinv,
                                           cfg.Cd, cfg.area);
    auto tot = computeTotalAttitudeDynamics(att, eci, I, Iinv, cfg, 0.0);

    check(nearlyEqual((tot.omega_dot - ref.omega_dot).norm(), 0.0, 1e-15),
          "torqueTotal: drag only matches dragAttitudeDynamics");
}

static void testTorqueTotal_matchesSRPOnly() {
    AttitudeState att = {identityQuat(), Vec3(0.0, 0.0, 0.0)};
    ECIState eci;
    eci.position << 6.778e6, 0.0, 0.0;
    eci.velocity << 0.0, 0.0, 0.0;

    Mat3 I    = diagI(10.0, 20.0, 30.0);
    Mat3 Iinv = diagIinv(10.0, 20.0, 30.0);

    TorqueModelConfig cfg;
    cfg.useGravityGradient = false;
    cfg.useDrag = false;
    cfg.useSRP = true;
    cfg.center_of_mass = Vec3(0.0, 0.05, 0.0);
    cfg.sun_dir_eci = Vec3(1.0, 0.0, 0.0);
    cfg.Cr = 1.5;
    cfg.area = 2.0;
    cfg.P_srp = 4.56e-6;

    auto ref = computeSRPAttitudeDynamics(att, eci, cfg.sun_dir_eci,
                                          cfg.center_of_mass, I, Iinv,
                                          cfg.Cr, cfg.area, cfg.P_srp);
    auto tot = computeTotalAttitudeDynamics(att, eci, I, Iinv, cfg, 0.0);

    check(nearlyEqual((tot.omega_dot - ref.omega_dot).norm(), 0.0, 1e-15),
          "torqueTotal: SRP only matches srpAttitudeDynamics");
}

static void testTorqueTotal_combinedEqualsSum() {
    AttitudeState att = {identityQuat(), Vec3(0.0, 0.0, 0.0)};
    ECIState eci;
    eci.position << 4.79e6, 4.79e6, 0.0;
    eci.velocity << 0.0, 7.8e3, 0.0;

    Mat3 I    = diagI(10.0, 20.0, 30.0);
    Mat3 Iinv = diagIinv(10.0, 20.0, 30.0);

    TorqueModelConfig cfg;
    cfg.useGravityGradient = true;
    cfg.useDrag = true;
    cfg.useSRP = true;
    cfg.center_of_mass = Vec3(0.1, 0.05, 0.0);
    cfg.sun_dir_eci = Vec3(1.0, 0.0, 0.0);
    cfg.Cd = 2.2;
    cfg.Cr = 1.5;
    cfg.area = 2.0;
    cfg.P_srp = 4.56e-6;

    auto gg  = computeGravGradAttitudeDynamics(att, eci, I, Iinv, 0.0);
    auto drg = computeDragAttitudeDynamics(att, eci, cfg.center_of_mass,
                                           I, Iinv, cfg.Cd, cfg.area);
    auto srp = computeSRPAttitudeDynamics(att, eci, cfg.sun_dir_eci,
                                           cfg.center_of_mass, I, Iinv,
                                           cfg.Cr, cfg.area, cfg.P_srp);

    Vec3 expected = gg.omega_dot + drg.omega_dot + srp.omega_dot;

    auto tot = computeTotalAttitudeDynamics(att, eci, I, Iinv, cfg, 0.0);

    check(nearlyEqual((tot.omega_dot - expected).norm(), 0.0, 1e-15),
          "attitudeTotal: combined model equals sum of component models");
}

static void testTorqueTotal_noneEnabledZeroOmegaDot() {
    AttitudeState att = {identityQuat(), Vec3(0.0, 0.0, 0.0)};
    ECIState eci;
    eci.position << 6.778e6, 0.0, 0.0;
    eci.velocity << 0.0, 7.8e3, 0.0;

    Mat3 I    = diagI(10.0, 20.0, 30.0);
    Mat3 Iinv = diagIinv(10.0, 20.0, 30.0);

    TorqueModelConfig cfg;
    cfg.useGravityGradient = false;
    cfg.useDrag = false;
    cfg.useSRP = false;

    auto tot = computeTotalAttitudeDynamics(att, eci, I, Iinv, cfg, 0.0);

    check(nearlyEqual(tot.omega_dot.norm(), 0.0, 1e-18),
          "attitudeTotal: no enabled torques gives zero omega_dot");
}

// ─────────────────────────────────────────────
//  Main
// ─────────────────────────────────────────────

int main() {
    std::cout << "\n=== Gravity-Gradient Attitude Dynamics ===\n";
    testGgrad_symmetricInertiaZeroTorque();
    testGgrad_principalAxisAlignedZeroTorque();
    testGgrad_inverseR3Scaling();
    testGgrad_qdotZeroForZeroOmega();
    testGgrad_qdotNonzeroForSpinning();
    testGgrad_spinOnPrincipalAxisNoPrecession();

    std::cout << "\n=== Drag Attitude Dynamics ===\n";
    testDrag_comAtOriginZeroTorque();
    testDrag_torqueSignWithComOffset();
    testDrag_scalesWithCd();
    testDrag_scalesWithArea();
    testDrag_zeroVelocityZeroTorque();
    testDrag_qdotZeroForZeroOmega();

    std::cout << "\n=== SRP Attitude Dynamics ===\n";
    testSRP_comAtOriginZeroTorque();
    testSRP_torqueSignWithComOffset();
    testSRP_analyticalMagnitude();
    testSRP_scalesWithP();
    testSRP_scalesWithCr();
    testSRP_scalesWithArea();
    testSRP_reverseSunDirReversesTorque();
    testSRP_qdotZeroForZeroOmega();

    std::cout << "\n=== Total Torque Model ===\n";
    testTorqueTotal_matchesGravityGradientOnly();
    testTorqueTotal_matchesDragOnly();
    testTorqueTotal_matchesSRPOnly();
    testTorqueTotal_combinedEqualsSum();
    testTorqueTotal_noneEnabledZeroOmegaDot();

    std::cout << "\n────────────────────────────\n";
    std::cout << "Results: " << g_passed << " passed, "
                             << g_failed << " failed\n\n";
    return g_failed > 0 ? 1 : 0;
}