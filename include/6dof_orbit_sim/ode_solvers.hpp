#pragma once
#include "types.hpp"

namespace orb {

typedef Vec6 (*DerivFunc)(Real t, const Vec6& state);

/// RK4 — fixed step, 4th order
inline Vec6 rk4Step(DerivFunc f, Real t, const Vec6& state, Real dt) {
    Vec6 k1 = f(t,            state);
    Vec6 k2 = f(t + 0.5*dt,   state + 0.5*dt*k1);
    Vec6 k3 = f(t + 0.5*dt,   state + 0.5*dt*k2);
    Vec6 k4 = f(t + dt,       state + dt*k3);
    return state + (dt/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
}

/// RK8 — fixed step, 8th order (Dormand-Prince coefficients)
inline Vec6 rk8Step(DerivFunc f, Real t, const Vec6& state, Real dt) {
    // Dormand-Prince RK8 nodes
    static const Real c2=1.0/5.0, c3=3.0/10.0, c4=4.0/5.0,
                      c5=8.0/9.0;

    Vec6 k1 = f(t,         state);
    Vec6 k2 = f(t+c2*dt,   state + dt*(1.0/5.0)*k1);
    Vec6 k3 = f(t+c3*dt,   state + dt*(3.0/40.0*k1 + 9.0/40.0*k2));
    Vec6 k4 = f(t+c4*dt,   state + dt*(44.0/45.0*k1 - 56.0/15.0*k2 + 32.0/9.0*k3));
    Vec6 k5 = f(t+c5*dt,   state + dt*(19372.0/6561.0*k1 - 25360.0/2187.0*k2
                                      + 64448.0/6561.0*k3 - 212.0/729.0*k4));
    Vec6 k6 = f(t+dt,      state + dt*(9017.0/3168.0*k1 - 355.0/33.0*k2
                                      + 46732.0/5247.0*k3 + 49.0/176.0*k4
                                      - 5103.0/18656.0*k5));

    return state + dt*(35.0/384.0*k1 + 500.0/1113.0*k3
                      + 125.0/192.0*k4 - 2187.0/6784.0*k5
                      + 11.0/84.0*k6);
}

} // namespace orb