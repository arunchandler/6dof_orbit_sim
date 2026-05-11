#include "6dof_orbit_sim/6dof_dynamics.hpp"
#include <cmath>
#include <stdexcept>

namespace orb {

SixDoFStateDot sixDoFDynamics(const SixDoFState& state,
                             Mat3& inertia_tensor,
                             Mat3& inertia_tensor_inv,
                             const ForceModelConfig& force_config,
                             const TorqueModelConfig& torque_config) {

    ECIStateDot eci_dot = computeTotalTranslationalDynamics(state.eci, force_config);
    AttitudeStateDot attitude_dot = computeTotalAttitudeDynamics(state.attitude, state.eci, inertia_tensor, inertia_tensor_inv, torque_config);

    return SixDoFStateDot{eci_dot, attitude_dot};

}

DerivFunc makeSixDoFDeriv(Mat3& inertia, Mat3& inertia_inv,
                           const ForceModelConfig& force_config,
                           const TorqueModelConfig& torque_config) {
    return [&](Real t, const VecX& state_vec) -> VecX {
        SixDoFState state = unpack6DoF(state_vec);
        SixDoFStateDot dot = sixDoFDynamics(state, inertia, inertia_inv,
                                             force_config, torque_config);
        return pack6DoFDot(dot);
    };
}

} // namespace orb