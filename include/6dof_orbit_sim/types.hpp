#pragma once
#include <Eigen/Dense>
#include <cmath>

namespace orb {

// Scalar and vector aliases
typedef double Real;
typedef Eigen::Matrix<Real, 3, 1> Vec3;
typedef Eigen::Matrix<Real, 6, 1> Vec6;
typedef Eigen::Matrix<Real, 3, 3> Mat3;

// Physical constants (WGS-84 / standard)
namespace constants {
    static const Real MU_EARTH = 3.986004418e14; // m³/s²
    static const Real R_EARTH  = 6378137.0;      // m
    static const Real J2       = 1.08262668e-3;  // dimensionless
    static const Real TWO_PI   = 2.0 * M_PI;
    static const Real AU       = 149597870700.0; // m
    static const Real C_LIGHT  = 299792458.0;     // m/s
}

// unit testing utilities
inline bool nearlyEqual(Real a, Real b, Real tol = 1e-6) {
    return std::fabs(a - b) < tol;
}

} // namespace orb