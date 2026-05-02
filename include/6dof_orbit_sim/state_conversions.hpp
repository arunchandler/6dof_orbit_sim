#pragma once
#include "elements.hpp"
#include "types.hpp"
#include "translational_state.hpp"

namespace orb {

/// Convert classical orbital elements -> ECI Cartesian state [m, m/s].
/// Returns an ECIState: [position, velocity] in SI units.
ECIState elementsToECI(const OrbitalElements& oe,
                   Real mu = constants::MU_EARTH);

/// Convert ECI Cartesian state -> classical orbital elements.
/// Handles circular and equatorial edge cases.
OrbitalElements eciToElements(const ECIState& state,
                               Real mu = constants::MU_EARTH);

/// True anomaly -> eccentric anomaly
Real trueToEccentric(Real ta, Real ecc);

/// Eccentric anomaly -> mean anomaly (Kepler)
Real eccentricToMean(Real E, Real ecc);

/// Mean anomaly -> eccentric anomaly (Newton-Raphson)
Real meanToEccentric(Real M, Real ecc, Real tol = 1e-12, int maxIter = 100);

} // namespace orb