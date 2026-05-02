#include "6dof_orbit_sim/state_conversions.hpp"
#include <stdexcept>
#include <cmath>

namespace orb {

ECIState elementsToECI(const OrbitalElements& oe, Real mu) {
    oe.validate();

    const Real e  = oe.ecc;
    const Real a  = oe.sma;
    const Real ta = oe.ta;

    // --- Perifocal (PQW) frame ---
    const Real p      = a * (1.0 - e * e);
    const Real r      = p / (1.0 + e * std::cos(ta));
    const Real sqMuP  = std::sqrt(mu / p);

    Vec3 r_pqw;
    r_pqw << r * std::cos(ta), r * std::sin(ta), 0.0;

    Vec3 v_pqw;
    v_pqw << -sqMuP * std::sin(ta), sqMuP * (e + std::cos(ta)), 0.0;

    // --- Rotation matrix PQW -> ECI (3-1-3: RAAN, inc, AoP) ---
    const Real cO = std::cos(oe.raan), sO = std::sin(oe.raan);
    const Real ci = std::cos(oe.inc),  si = std::sin(oe.inc);
    const Real cw = std::cos(oe.aop),  sw = std::sin(oe.aop);

    Mat3 R;
    R << cO*cw - sO*sw*ci,  -cO*sw - sO*cw*ci,  sO*si,
         sO*cw + cO*sw*ci,  -sO*sw + cO*cw*ci, -cO*si,
         sw*si,               cw*si,              ci;

    ECIState state;
    state.position = R * r_pqw;
    state.velocity = R * v_pqw;
    return state;
}

OrbitalElements eciToElements(const ECIState& state, Real mu) {
    const Vec3 r_vec = state.position;
    const Vec3 v_vec = state.velocity;

    const Real r = r_vec.norm();
    const Real v = v_vec.norm();

    const Vec3 h_vec = r_vec.cross(v_vec);
    const Real h     = h_vec.norm();

    Vec3 z_hat;
    z_hat << 0.0, 0.0, 1.0;
    const Vec3 n_vec = z_hat.cross(h_vec);
    const Real n     = n_vec.norm();

    const Vec3 e_vec = ((v*v - mu/r) * r_vec - r_vec.dot(v_vec) * v_vec) / mu;

    OrbitalElements oe;
    oe.ecc = e_vec.norm();
    oe.sma = 1.0 / (2.0/r - v*v/mu);
    oe.inc = std::acos(clamp(h_vec(2) / h, -1.0, 1.0));

    const Real eps = 1e-10;

    // RAAN
    if (n > eps) {
        Real val = clamp(n_vec(0) / n, -1.0, 1.0);
        oe.raan = wrapTwoPI(std::acos(val) * (n_vec(1) >= 0.0 ? 1.0 : -1.0));
    } else {
        oe.raan = 0.0;
    }

    // AoP
    if (n > eps && oe.ecc > eps) {
        Real val = clamp(n_vec.dot(e_vec) / (n * oe.ecc), -1.0, 1.0);
        oe.aop = wrapTwoPI(std::acos(val) * (e_vec(2) >= 0.0 ? 1.0 : -1.0));
    } else {
        oe.aop = 0.0;
    }

    // True anomaly
    if (oe.ecc > eps) {
        Real val = clamp(e_vec.dot(r_vec) / (oe.ecc * r), -1.0, 1.0);
        oe.ta = wrapTwoPI(std::acos(val) * (r_vec.dot(v_vec) >= 0.0 ? 1.0 : -1.0));
    } else {
        Real val = clamp(n_vec.dot(r_vec) / (n * r), -1.0, 1.0);
        oe.ta = wrapTwoPI(std::acos(val) * (r_vec(2) >= 0.0 ? 1.0 : -1.0));
    }

    return oe;
}

Real trueToEccentric(Real ta, Real ecc) {
    return 2.0 * std::atan2(
        std::sqrt(1.0 - ecc) * std::sin(ta / 2.0),
        std::sqrt(1.0 + ecc) * std::cos(ta / 2.0)
    );
}

Real eccentricToMean(Real E, Real ecc) {
    return E - ecc * std::sin(E);
}

Real meanToEccentric(Real M, Real ecc, Real tol, int maxIter) {
    Real E = (ecc > 0.8) ? M_PI : M;
    for (int i = 0; i < maxIter; ++i) {
        const Real dE = (M - E + ecc * std::sin(E)) / (1.0 - ecc * std::cos(E));
        E += dE;
        if (std::fabs(dE) < tol) return E;
    }
    throw std::runtime_error("Kepler equation did not converge");
}

} // namespace orb