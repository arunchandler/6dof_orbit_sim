#include "6dof_orbit_sim/6dof_dynamics.hpp"
#include <cmath>
#include <stdexcept>

namespace orb {

SixDoFStateDot sixDoFDynamics(const SixDoFState& state,
                             Mat3& inertia_tensor,
                             Mat3& inertia_tensor_inv,
                             const ForceModelConfig& force_config,
                             const TorqueModelConfig& torque_config,
                             ActuatorSuite<4,3>& suite,
                             ActuatorCommands<4,3>& cmds,
                             const Quaternion& q_desired,
                             const QuaternionPDController& ctrl,
                             Real t_sec,
                             Real dt) {

    Vec3 control_torque = ctrl.compute_torque(state.attitude, q_desired);
    cmds.tau_rwa_cmd = control_torque; // For simplicity, send all control torque to RWAs
    const Vec3 h_rw = suite.rwa.total_momentum();
    Real k_desat = 0.0001; // Desaturation gain (tunable)
    cmds.tau_mtq_cmd = -k_desat * h_rw;
    ECIStateDot eci_dot = computeTotalTranslationalDynamics(state.eci, force_config);
    AttitudeStateDot attitude_dot = computeTotalAttitudeDynamics(state.attitude, state.eci, inertia_tensor, inertia_tensor_inv, torque_config, suite, cmds, t_sec, dt);

    return SixDoFStateDot{eci_dot, attitude_dot};

}

DerivFunc makeSixDoFDeriv(Mat3& inertia, Mat3& inertia_inv,
                           const ForceModelConfig& force_config,
                           const TorqueModelConfig& torque_config,
                           ActuatorSuite<4,3>& suite,
                           ActuatorCommands<4,3>& cmds,
                           const Quaternion& q_desired,
                           const QuaternionPDController& ctrl,
                           Real dt) {
    return [&, dt](Real t, const VecX& state_vec) -> VecX {
        SixDoFState state = unpack6DoF(state_vec);
        
        SixDoFStateDot dot = sixDoFDynamics(state, inertia, inertia_inv,
                                             force_config, torque_config, suite, cmds, q_desired, ctrl, t, dt);
        return pack6DoFDot(dot);
    };
}

} // namespace orb