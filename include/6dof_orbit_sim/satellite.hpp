#pragma once
#include "types.hpp"
#include "attitude_state.hpp"
#include "translational_state.hpp"
#include "elements.hpp"

namespace orb {

class Satellite {
public:
    Real mu = constants::MU_EARTH; // Gravitational parameter (default: Earth's mu)
    ECIState state;           // Current state in ECI frame
    AttitudeState attitude;   // Current attitude state (quaternion and angular velocity)
    OrbitalElements elements; // Current orbital elements
    Real mass;             // Mass of the satellite [kg]
    Real dry_mass;         // Dry mass of the satellite [kg]
    Real area;             // Reference area for drag and SRP [m²]
    Real Cd;               // Drag coefficient
    Real Cr;               // Reflectivity coefficient for SRP
    Mat3 inertia_tensor;   // Inertia tensor in body frame [kg·m²]
    Mat3 inertia_tensor_inv; // Inverse of inertia tensor for dynamics calculations
    Vec3 center_of_mass;  // Center of mass offset from geometric center in body frame [m]
                            // Going to assume center of pressure is at the geometric center for simplicity

    void updateElements() {
        elements = eciToElements(state, mu);
    }

    void updateECIState() {
        state = elementsToECI(elements, mu);
    }

};

}