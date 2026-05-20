#pragma once
#include "types.hpp"
#include "ode_solvers.hpp"
#include "6dof_state.hpp"
#include "state_conversions.hpp"

namespace orb {
    
SixDoFStateMat propagate6DoF(const SixDoFState& initial_state, Real t0, Real tf, Real dt, DerivFunc deriv_func) {
    int steps = static_cast<int>((tf - t0) / dt) + 1;
    SixDoFStateMat states(SIX_DOF_STATE_DIM, steps);
    states.col(0) = pack6DoF(initial_state);
    for (int i = 1; i < steps; ++i) {
        VecX next = rk4Step(deriv_func, t0 + (i-1)*dt, states.col(i-1), dt);
        next.segment<4>(6).normalize(); // Normalize quaternion part to prevent drift
        states.col(i) = next;
    }
    return states;
}

TranslationalStateMat propagateTranslationalState(const ECIState& initial_state, Real t0, Real tf, Real dt, DerivFunc deriv_func) {
    int steps = static_cast<int>((tf - t0) / dt) + 1;
    TranslationalStateMat states(TRANLATIONAL_STATE_DIM, steps);
    states.col(0) = packTranslational(initial_state);
    for (int i = 1; i < steps; ++i) {
        states.col(i) = rk4Step(deriv_func, t0 + (i-1)*dt, states.col(i-1), dt);
    }
    return states;
}

AttitudeStateMat propagateAttitudeState(const AttitudeState& initial_state, Real t0, Real tf, Real dt, DerivFunc deriv_func) {
    int steps = static_cast<int>((tf - t0) / dt) + 1;
    AttitudeStateMat states(ATTITUDE_STATE_DIM, steps);
    states.col(0) = packAttitude(initial_state);
    for (int i = 1; i < steps; ++i) {
        VecX next = rk4Step(deriv_func, t0 + (i-1)*dt, states.col(i-1), dt);
        next.segment<4>(6).normalize(); // Normalize quaternion part to prevent drift
        states.col(i) = next;
    }
    return states;
}

} // namespace orb
