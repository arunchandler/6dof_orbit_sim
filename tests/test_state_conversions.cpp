#include "6dof_orbit_sim/state_conversions.hpp"
#include <iostream>
#include <cmath>

using namespace orb;

static void testRoundTrip() {
    OrbitalElements oe_in;
    oe_in.sma  = 7000000.0;   // 7000 km
    oe_in.ecc  = 0.01;
    oe_in.inc  = 0.5;         // ~28.6 deg
    oe_in.raan = 1.1;
    oe_in.aop  = 0.8;
    oe_in.ta   = 2.0;

    ECIState state        = elementsToECI(oe_in);
    OrbitalElements oe_out = eciToElements(state);

    bool pass =
        nearlyEqual(oe_in.sma,  oe_out.sma,  1e-3) &&
        nearlyEqual(oe_in.ecc,  oe_out.ecc,  1e-9) &&
        nearlyEqual(oe_in.inc,  oe_out.inc,  1e-9) &&
        nearlyEqual(oe_in.raan, oe_out.raan, 1e-9) &&
        nearlyEqual(oe_in.aop,  oe_out.aop,  1e-9) &&
        nearlyEqual(oe_in.ta,   oe_out.ta,   1e-9);

    std::cout << (pass ? "PASS" : "FAIL") << "  round-trip elements->ECI->elements\n";
}

static void testKeplerEquation() {
    Real M    = 1.2;
    Real ecc  = 0.3;
    Real E    = meanToEccentric(M, ecc);
    Real M2   = eccentricToMean(E, ecc);

    bool pass = nearlyEqual(M, M2, 1e-12);
    std::cout << (pass ? "PASS" : "FAIL") << "  Kepler equation round-trip\n";
}

int main() {
    testRoundTrip();
    testKeplerEquation();
    return 0;
}