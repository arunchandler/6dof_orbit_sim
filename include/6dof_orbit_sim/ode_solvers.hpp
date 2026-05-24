#pragma once
#include "types.hpp"
#include "translational_state.hpp"
#include "attitude_state.hpp"
#include "6dof_state.hpp"
#include "state_conversions.hpp"

namespace orb {

/// RK4 — fixed step, 4th order
inline VecX rk4Step(DerivFunc f, Real t, const VecX& state, Real dt) {
    VecX k1 = f(t,            state);
    VecX k2 = f(t + 0.5*dt,   state + 0.5*dt*k1);
    VecX k3 = f(t + 0.5*dt,   state + 0.5*dt*k2);
    VecX k4 = f(t + dt,       state + dt*k3);
    return state + (dt/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
}

inline VecX eulerStep(DerivFunc f, Real t, const VecX& state, Real dt) {
    return state + dt * f(t, state);
}

} // namespace orb