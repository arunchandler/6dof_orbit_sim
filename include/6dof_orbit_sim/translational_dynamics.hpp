#pragma once
#include "types.hpp"
#include "elements.hpp"
#include "state_conversions.hpp"
#include "translational_state.hpp"

namespace orb {

/// @brief Two-body gravitational acceleration.
/// @param state Current state vector [rx, ry, rz, vx, vy, vz]
/// @param mu Gravitational parameter (default: Earth's mu)
/// @return  State derivative due to two-body gravity [0, 0, 0, ax, ay, az]
ECIStateDot accTwoBody(const ECIState& state, Real mu = constants::MU_EARTH);

/// @brief Perturbation due to Earth's oblateness (J2 effect).
/// @param state Current state vector [rx, ry, rz, vx, vy, vz]
/// @param mu Gravitational parameter (default: Earth's mu)
/// @param J2 Dimensionless J2 coefficient (default: Earth's J2)
/// @param R Reference radius (default: Earth's mean radius)
/// @return  State derivative due to J2 perturbation [0, 0, 0, ax, ay, az]
ECIStateDot accJ2(const ECIState& state, Real mu = constants::MU_EARTH, Real J2 = constants::J2, Real R = constants::R_EARTH);

/// @brief Drag force due to atmospheric drag.
/// @param state Current state vector [rx, ry, rz, vx, vy, vz]
/// @param mu Gravitational parameter (default: Earth's mu)
/// @param Cd Drag coefficient (default: 2.2)
/// @param A Reference area (default: 1.0)
/// @param m Mass (default: 100.0)
/// @return  State derivative due to drag force [0, 0, 0, ax, ay, az]
ECIStateDot accDrag(const ECIState& state, Real mu = constants::MU_EARTH, Real Cd = 2.2, Real A = 1.0, Real m = 100.0);

/// @brief Solar radiation pressure force.
/// @param state Current state vector [rx, ry, rz, vx, vy, vz]
/// @param mu Gravitational parameter (default: Earth's mu)
/// @param Cr Reflectivity coefficient (default: 1.5)
/// @param A Reference area (default: 1.0)
/// @param m Mass (default: 100.0)
/// @return  State derivative due to solar radiation pressure [0, 0, 0, ax, ay, az]
ECIStateDot accSRP(const ECIState& state, Real mu = constants::MU_EARTH, Real Cr = 1.5, Real A = 1.0, Real m = 100.0);

/// @brief Configuration for selecting which forces to include in the total acceleration calculation.
struct ForceModelConfig {
    bool useJ2  = true;
    bool useDrag = false;
    bool useSRP  = false;
    Real Cd = 2.2;
    Real Cr = 1.5;
    Real A = 1.0;   // m²
    Real m = 100.0; // kg
};

/// @brief Compute total acceleration from selected force model components.
/// @param state Current state vector [rx, ry, rz, vx, vy, vz]
/// @param config Force model configuration
/// @return  State derivative due to selected forces [0, 0, 0, ax, ay, az]
ECIStateDot accTotal(const ECIState& state, const ForceModelConfig& config);

/// @brief Compute the time derivative of orbital elements using Gauss's variational equations.
/// @param oe Current orbital elements
/// @param config Force model configuration
/// @param mu Gravitational parameter (default: Earth's mu)
/// @return  Time derivative of orbital elements [sma_dot, ecc_dot, inc_dot, raan_dot, aop_dot, ta_dot]
OrbitalElements dynamicsGVE(const OrbitalElements& oe, const ForceModelConfig& config, Real mu = constants::MU_EARTH);

} // namespace orb