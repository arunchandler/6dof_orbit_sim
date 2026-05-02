#pragma once
#include "types.hpp"
#include <stdexcept>

namespace orb {

/// Classical Keplerian orbital elements (SI units — meters and radians)
struct OrbitalElements {
    Real sma;   ///< Semi-major axis          [m]
    Real ecc;   ///< Eccentricity             [-]
    Real inc;   ///< Inclination              [rad]
    Real raan;  ///< Right ascension of AN    [rad]
    Real aop;   ///< Argument of periapsis    [rad]
    Real ta;    ///< True anomaly             [rad]

    void validate() const {
        if (sma <= 0.0)
            throw std::invalid_argument("sma must be positive");
        if (ecc < 0.0 || ecc >= 1.0)
            throw std::invalid_argument("ecc must be in [0, 1) for elliptic orbits");
        if (inc < 0.0 || inc > M_PI)
            throw std::invalid_argument("inc must be in [0, pi]");
    }
};

/// Wrap an angle to [0, 2pi)
inline Real wrapTwoPI(Real angle) {
    angle = std::fmod(angle, constants::TWO_PI);
    return angle < 0.0 ? angle + constants::TWO_PI : angle;
}

/// Clamp a value to [lo, hi]
inline Real clamp(Real val, Real lo, Real hi) {
    return val < lo ? lo : (val > hi ? hi : val);
}

} // namespace orb