#include "6dof_orbit_sim/translational_dynamics.hpp"
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

/// Build a simple circular LEO state at the given radius [m]
static ECIState circularLEO(Real r) {
    Real v = std::sqrt(constants::MU_EARTH / r);
    ECIState state;
    state.position << r, 0.0, 0.0;
    state.velocity << 0.0, v, 0.0;
    return state;
}

// ─────────────────────────────────────────────
//  1. Two-body gravity
// ─────────────────────────────────────────────

static void testTwoBody_magnitude() {
    // For circular orbit, |a| = mu/r²
    Real r_mag  = 7000000.0;
    ECIState state  = circularLEO(r_mag);
    ECIStateDot deriv  = accTwoBody(state);

    Vec3 acc    = deriv.acceleration;
    Real a_mag  = acc.norm();
    Real a_exp  = constants::MU_EARTH / (r_mag * r_mag);

    check(nearlyEqual(a_mag, a_exp, 1e-3),
          "twoBody: acceleration magnitude matches mu/r^2");
}

static void testTwoBody_direction() {
    // Acceleration must point toward origin (anti-parallel to r)
    ECIState state = circularLEO(7000000.0);
    ECIStateDot deriv = accTwoBody(state);

    Vec3 r_hat = state.position.normalized();
    Vec3 a_hat = deriv.acceleration.normalized();

    // dot product should be -1
    Real dot = r_hat.dot(a_hat);
    check(nearlyEqual(dot, -1.0, 1e-9),
          "twoBody: acceleration points toward Earth center");
}

static void testTwoBody_energyConservation() {
    // In two-body, specific energy E = v²/2 - mu/r is constant.
    // Check the time derivative dE/dt = v·a - mu*(r·v)/r³ = 0
    ECIState state = circularLEO(7000000.0);
    ECIStateDot deriv = accTwoBody(state);
    Vec3  r     = state.position;
    Vec3  v     = state.velocity;
    Vec3  a     = deriv.acceleration;
    Real  r_mag = r.norm();

    Real dEdt = v.dot(a) + constants::MU_EARTH * r.dot(v) / std::pow(r_mag, 3);
    check(nearlyEqual(dEdt, 0.0, 1e-4),
          "twoBody: specific orbital energy is conserved (dE/dt = 0)");
}

static void testTwoBody_customMu() {
    // Doubling mu should double the acceleration magnitude
    Real r_mag  = 7000000.0;
    ECIState state  = circularLEO(r_mag);
    ECIStateDot d1     = accTwoBody(state, constants::MU_EARTH);
    ECIStateDot d2     = accTwoBody(state, 2.0 * constants::MU_EARTH);

    Real ratio  = d2.acceleration.norm() / d1.acceleration.norm();
    check(nearlyEqual(ratio, 2.0, 1e-9),
          "twoBody: acceleration scales linearly with mu");
}

// ─────────────────────────────────────────────
//  2. J2 perturbation
// ─────────────────────────────────────────────

static void testJ2_smallerThanTwoBody() {
    // J2 perturbation should be ~1000x smaller than two-body
    ECIState state  = circularLEO(7000000.0);
    Real a_2b   = accTwoBody(state).acceleration.norm();
    Real a_j2   = accJ2(state).acceleration.norm();

    check(a_j2 < a_2b * 1e-2,
          "J2: perturbation is at least 100x smaller than two-body");
}

static void testJ2_equatorialSymmetry() {
    // On the equatorial plane (z=0), J2 acceleration has no z-component
    Real r = 7000000.0;
    ECIState state;
    state.position << r, 0.0, 0.0;
    state.velocity << 0.0, std::sqrt(constants::MU_EARTH / r), 0.0;

    ECIStateDot deriv = accJ2(state);
    check(nearlyEqual(deriv.acceleration(2), 0.0, 1e-6),
          "J2: no z-acceleration for equatorial orbit");
}

static void testJ2_polarEnhancement() {
    // At the pole (r along z), radial J2 acceleration is larger than equatorial
    Real r = 7000000.0;

    ECIState eq_state;
    eq_state.position << r, 0.0, 0.0;
    eq_state.velocity << 0.0, 0.0, 0.0;

    ECIState pol_state;
    pol_state.position << 0.0, 0.0, r;
    pol_state.velocity << 0.0, 0.0, 0.0;

    Real a_eq  = accJ2(eq_state).acceleration.norm();
    Real a_pol = accJ2(pol_state).acceleration.norm();

    check(a_pol > a_eq,
          "J2: radial acceleration is stronger at poles than equator");
}

static void testJ2_vanishesWithZeroJ2() {
    ECIState state = circularLEO(7000000.0);
    ECIStateDot deriv = accJ2(state, constants::MU_EARTH, 0.0, constants::R_EARTH);

    check(nearlyEqual(deriv.acceleration.norm(), 0.0, 1e-12),
          "J2: acceleration is zero when J2 = 0");
}

// ─────────────────────────────────────────────
//  3. Atmospheric drag
// ─────────────────────────────────────────────

static void testDrag_outputStructure() {
    ECIState state = circularLEO(7000000.0);
    ECIStateDot deriv = accDrag(state);

    check(nearlyEqual(deriv.acceleration(0), 0.0) &&
          nearlyEqual(deriv.acceleration(1), 0.0) &&
          nearlyEqual(deriv.acceleration(2), 0.0),
          "drag: position rows of derivative are zero");
}

static void testDrag_opposesVelocity() {
    // Drag must oppose the velocity vector
    ECIState state = circularLEO(7000000.0);
    ECIStateDot deriv = accDrag(state);

    Vec3 v_hat = state.velocity.normalized();
    Vec3 a_hat = deriv.acceleration.normalized();

    Real dot = v_hat.dot(a_hat);
    check(dot < 0.0,
          "drag: acceleration opposes velocity");
}

static void testDrag_smallerThanTwoBody() {
    ECIState state = circularLEO(7000000.0);
    ECIStateDot deriv_2b = accTwoBody(state);
    ECIStateDot deriv_drag = accDrag(state);

    Real a_2b = deriv_2b.acceleration.norm();
    Real a_drag = deriv_drag.acceleration.norm();

    check(a_drag < a_2b * 1e-3,
          "drag: perturbation is at least 1000x smaller than two-body");
}

static void testDrag_scalesWithArea() {
    // Doubling area should double drag acceleration
    ECIState state = circularLEO(7000000.0);
    ECIStateDot d1    = accDrag(state, constants::MU_EARTH, 2.2, 1.0, 100.0);
    ECIStateDot d2    = accDrag(state, constants::MU_EARTH, 2.2, 2.0, 100.0);

    Real ratio = d2.acceleration.norm() / d1.acceleration.norm();
    check(nearlyEqual(ratio, 2.0, 1e-6),
          "drag: acceleration scales linearly with area");
}

static void testDrag_scalesInverselyWithMass() {
    // Doubling mass should halve drag acceleration
    ECIState state = circularLEO(7000000.0);
    ECIStateDot d1    = accDrag(state, constants::MU_EARTH, 2.2, 1.0, 100.0);
    ECIStateDot d2    = accDrag(state, constants::MU_EARTH, 2.2, 1.0, 200.0);

    Real ratio = d2.acceleration.norm() / d1.acceleration.norm();
    check(nearlyEqual(ratio, 0.5, 1e-6),
          "drag: acceleration scales inversely with mass");
}

// ─────────────────────────────────────────────
//  4. Solar radiation pressure
// ─────────────────────────────────────────────

static void testSRP_outputStructure() {
    ECIState state = circularLEO(7000000.0);
    ECIStateDot deriv = accSRP(state);

    check(nearlyEqual(deriv.acceleration(0), 0.0) &&
          nearlyEqual(deriv.acceleration(1), 0.0) &&
          nearlyEqual(deriv.acceleration(2), 0.0),
          "SRP: position rows of derivative are zero");
}

static void testSRP_smallerThanTwoBody() {
    ECIState state = circularLEO(7000000.0);
    ECIStateDot deriv_2b = accTwoBody(state);
    ECIStateDot deriv_srp = accSRP(state);

    Real a_2b = deriv_2b.acceleration.norm();
    Real a_srp = deriv_srp.acceleration.norm();

    check(a_srp < a_2b * 1e-4,
          "SRP: perturbation is much smaller than two-body");
}

static void testSRP_scalesWithArea() {
    ECIState state = circularLEO(7000000.0);
    ECIStateDot d1    = accSRP(state, constants::MU_EARTH, 1.5, 1.0, 100.0);
    ECIStateDot d2    = accSRP(state, constants::MU_EARTH, 1.5, 2.0, 100.0);

    Real ratio = d2.acceleration.norm() / d1.acceleration.norm();
    check(nearlyEqual(ratio, 2.0, 1e-6),
          "SRP: acceleration scales linearly with area");
}

static void testSRP_scalesInverselyWithMass() {
    ECIState state = circularLEO(7000000.0);
    ECIStateDot d1    = accSRP(state, constants::MU_EARTH, 1.5, 1.0, 100.0);
    ECIStateDot d2    = accSRP(state, constants::MU_EARTH, 1.5, 1.0, 200.0);

    Real ratio = d2.acceleration.norm() / d1.acceleration.norm();
    check(nearlyEqual(ratio, 0.5, 1e-6),
          "SRP: acceleration scales inversely with mass");
}

// ─────────────────────────────────────────────
//  5. Total acceleration / ForceModelConfig
// ─────────────────────────────────────────────

static void testTotal_twoBodyOnly() {
    // With J2/drag/SRP off, total should equal two-body exactly
    ForceModelConfig cfg;
    cfg.useJ2   = false;
    cfg.useDrag = false;
    cfg.useSRP  = false;

    ECIState state  = circularLEO(7000000.0);
    ECIStateDot d_2b   = accTwoBody(state);
    ECIStateDot d_tot  = accTotal(state, cfg);

    Real diff = (d_tot.acceleration - d_2b.acceleration).norm();
    check(nearlyEqual(diff, 0.0, 1e-9),
          "total: two-body-only config matches accTwoBody exactly");
}

static void testTotal_j2AddsToTwoBody() {
    ForceModelConfig cfg_no_j2, cfg_j2;
    cfg_no_j2.useJ2 = false;  cfg_no_j2.useDrag = false;  cfg_no_j2.useSRP = false;
    cfg_j2.useJ2    = true;   cfg_j2.useDrag    = false;  cfg_j2.useSRP    = false;

    ECIState state   = circularLEO(7000000.0);
    ECIStateDot d_no_j2 = accTotal(state, cfg_no_j2);
    ECIStateDot d_j2    = accTotal(state, cfg_j2);

    Real a_no_j2 = d_no_j2.acceleration.norm();
    Real a_j2    = d_j2.acceleration.norm();

    check(a_j2 != a_no_j2,
          "total: enabling J2 changes the acceleration magnitude");
}

static void testTotal_dragReducesEnergy() {
    // With drag on, the acceleration component along velocity should be negative
    ForceModelConfig cfg;
    cfg.useJ2   = false;
    cfg.useDrag = true;
    cfg.useSRP  = false;

    ECIState state  = circularLEO(7000000.0);
    ECIStateDot deriv  = accTotal(state, cfg);
    Vec3 v_hat  = state.velocity.normalized();
    Real along_v = deriv.acceleration.dot(v_hat);

    check(along_v < 0.0,
          "total: drag config produces negative along-track acceleration");
}

static void testTotal_allModels() {
    // Smoke test — all models on, no crash, acceleration is nonzero
    ForceModelConfig cfg;
    cfg.useJ2   = true;
    cfg.useDrag = true;
    cfg.useSRP  = true;

    ECIState state  = circularLEO(7000000.0);
    ECIStateDot deriv  = accTotal(state, cfg);
    Real a_mag  = deriv.acceleration.norm();

    check(a_mag > 0.0,
          "total: all models enabled produces nonzero acceleration");
}

// ─────────────────────────────────────────────
//  6. Gauss's Variational Equations
// ───────────────────────────────────────────── 

static void testGVE_twoBodyOnly_keplerian() {
    ForceModelConfig cfg;
    cfg.useJ2   = false;
    cfg.useDrag = false;
    cfg.useSRP  = false;

    OrbitalElements oe;
    oe.sma  = 7000000.0;
    oe.ecc  = 0.01;
    oe.inc  = 0.1;
    oe.raan = 0.2;
    oe.aop  = 0.3;
    oe.ta   = 0.4;

    OrbitalElements dot = dynamicsGVE(oe, cfg, constants::MU_EARTH);

    check(nearlyEqual(dot.sma,  0.0, 1e-9) &&
          nearlyEqual(dot.ecc,  0.0, 1e-9) &&
          nearlyEqual(dot.inc,  0.0, 1e-9) &&
          nearlyEqual(dot.raan, 0.0, 1e-9) &&
          nearlyEqual(dot.aop,  0.0, 1e-9) &&
          dot.ta > 0.0,
          "GVE: elements constant under pure two-body");

}

static void testGVE_circularOrbit_noNaNs() {
    ForceModelConfig cfg;
    cfg.useJ2 = true;  // include perturbation to exercise logic

    OrbitalElements oe;
    oe.sma  = 7000000.0;
    oe.ecc  = 0.0;     // circular
    oe.inc  = 0.3;
    oe.raan = 0.2;
    oe.aop  = 0.0;
    oe.ta   = 1.0;

    OrbitalElements dot = dynamicsGVE(oe, cfg, constants::MU_EARTH);

    check(std::isfinite(dot.sma) &&
          std::isfinite(dot.ecc) &&
          std::isfinite(dot.inc) &&
          std::isfinite(dot.raan) &&
          std::isfinite(dot.aop) &&
          std::isfinite(dot.ta),
          "GVE: circular orbit produces finite derivatives");

    check(nearlyEqual(dot.aop, 0.0, 1e-12),
          "GVE: argument of perigee frozen for circular orbit");
}

static void testGVE_equatorialOrbit_noNaNs() {
    ForceModelConfig cfg;
    cfg.useJ2 = true;

    OrbitalElements oe;
    oe.sma  = 7000000.0;
    oe.ecc  = 0.01;
    oe.inc  = 0.0;   // equatorial
    oe.raan = 0.0;
    oe.aop  = 0.3;
    oe.ta   = 0.5;

    OrbitalElements dot = dynamicsGVE(oe, cfg, constants::MU_EARTH);

    check(std::isfinite(dot.raan),
          "GVE: equatorial orbit RAAN derivative is finite");
}

static void testGVE_drag_reducesSMA() {
    ForceModelConfig cfg;
    cfg.useDrag = true;
    cfg.useJ2   = false;
    cfg.useSRP  = false;

    OrbitalElements oe;
    oe.sma  = 7000000.0;
    oe.ecc  = 0.001;
    oe.inc  = 0.1;
    oe.raan = 0.0;
    oe.aop  = 0.0;
    oe.ta   = 0.0;

    OrbitalElements dot = dynamicsGVE(oe, cfg, constants::MU_EARTH);

    check(dot.sma < 0.0,
          "GVE: drag decreases semi-major axis");
}

static void testGVE_J2_affectsRAAN() {
    ForceModelConfig cfg;
    cfg.useJ2   = true;
    cfg.useDrag = false;
    cfg.useSRP  = false;

    OrbitalElements oe;
    oe.sma  = 7000000.0;
    oe.ecc  = 0.01;
    oe.inc  = 0.75;   // non-equatorial
    oe.raan = 0.0;
    oe.aop  = 0.3;
    oe.ta   = 0.4;

    OrbitalElements dot = dynamicsGVE(oe, cfg, constants::MU_EARTH);

    check(std::fabs(dot.raan) > 0.0,
          "GVE: J2 produces nodal precession (RAAN rate nonzero)");
}

// ─────────────────────────────────────────────
//  Main
// ─────────────────────────────────────────────

int main() {
    std::cout << "\n=== Two-Body ===\n";
    testTwoBody_magnitude();
    testTwoBody_direction();
    testTwoBody_energyConservation();
    testTwoBody_customMu();

    std::cout << "\n=== J2 ===\n";
    testJ2_smallerThanTwoBody();
    testJ2_equatorialSymmetry();
    testJ2_polarEnhancement();
    testJ2_vanishesWithZeroJ2();

    std::cout << "\n=== Drag ===\n";
    testDrag_outputStructure();
    testDrag_opposesVelocity();
    testDrag_smallerThanTwoBody();
    testDrag_scalesWithArea();
    testDrag_scalesInverselyWithMass();

    std::cout << "\n=== SRP ===\n";
    testSRP_outputStructure();
    testSRP_smallerThanTwoBody();
    testSRP_scalesWithArea();
    testSRP_scalesInverselyWithMass();

    std::cout << "\n=== Total Acceleration ===\n";
    testTotal_twoBodyOnly();
    testTotal_j2AddsToTwoBody();
    testTotal_dragReducesEnergy();
    testTotal_allModels();

    std::cout << "\n=== GVE ===\n";
    testGVE_twoBodyOnly_keplerian();
    testGVE_circularOrbit_noNaNs();
    testGVE_equatorialOrbit_noNaNs();
    testGVE_drag_reducesSMA();
    testGVE_J2_affectsRAAN();

    std::cout << "\n────────────────────────────\n";
    std::cout << "Results: " << g_passed << " passed, "
                             << g_failed << " failed\n\n";
    return g_failed > 0 ? 1 : 0;
}