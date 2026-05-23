#include "6dof_orbit_sim/6dof_orbit_sim.hpp"
#include <iostream>
#include <cmath>
#include <string>
#include <iomanip>

using namespace orb;

int main() {

    // Sun synchronous orbit elements
    OrbitalElements oe_init;
    oe_init.sma = 6.7781e6;
    oe_init.ecc = 0.001;
    oe_init.inc = 97.04 * M_PI / 180.0;
    oe_init.raan = 0.0;
    oe_init.aop = 0.0;
    oe_init.ta = 0.0;

    // initial attitude state: identity quaternion, zero angular velocity
    // intertial to body
    AttitudeState att_init;
    att_init.q.w = 1.0;
    att_init.q.x = 0.0;
    att_init.q.y = 0.0;
    att_init.q.z = 0.0;
    att_init.omega = Vec3(0.0, 0.0, 0.0);

    Real mu = constants::MU_EARTH;
    ECIState state = elementsToECI(oe_init, mu);

    SixDoFState sixdof_init;
    sixdof_init.eci = state;
    sixdof_init.attitude = att_init;

    TorqueModelConfig torque_config;
    torque_config.useGravityGradient = true;
    torque_config.useDrag = true;
    torque_config.useSRP = true;
    torque_config.useActuators = true;
    torque_config.center_of_mass = Vec3(0.1, 0.05, 0.0);

    ForceModelConfig force_config;
    force_config.useDrag = true;
    force_config.useJ2 = true;
    force_config.useSRP = true;

    // controller
    Vec3 sun_dir_eci = Vec3(1.0, 0.0, 0.0); // Simplified constant sun direction for testing
    torque_config.sun_dir_eci = sun_dir_eci;
    Quaternion q_d = make_desired_quaternion(Vec3::UnitX(), sun_dir_eci, Vec3::UnitY(), Vec3::UnitZ());
    QuaternionPDController ctrl{
        // .k_p = 0.027,
        // .k_d = 1.35
        .k_p = 0.0,
        .k_d = 0.0,
    };

    // actuators 
    ActuatorSuite<4,3> suite;
    const Real beta = std::atan(std::sqrt(2.0)); // 54.74 deg from z-axis
    const Real sb   = std::sin(beta);             // ≈ 0.8165
    const Real cb   = std::cos(beta);             // ≈ 0.5774

    ReactionWheel rw_proto;
    rw_proto.inertia       = 0.05;    // kg·m²   (medium-class wheel)
    rw_proto.max_torque    = 0.5;     // N·m
    rw_proto.max_speed     = 500.0;   // rad/s   (~4800 RPM)
    rw_proto.friction_coef = 1.0e-4;  // N·m·s/rad
    rw_proto.omega_w       = 0.0;     // start at rest

    // Axes at 45°, 135°, 225°, 315° around z
    suite.rwa.wheels[0]             = rw_proto;
    suite.rwa.wheels[0].spin_axis_b = Vec3(sb * std::cos(1*M_PI/4),
                                            sb * std::sin(1*M_PI/4), cb).normalized();

    suite.rwa.wheels[1]             = rw_proto;
    suite.rwa.wheels[1].spin_axis_b = Vec3(sb * std::cos(3*M_PI/4),
                                            sb * std::sin(3*M_PI/4), cb).normalized();

    suite.rwa.wheels[2]             = rw_proto;
    suite.rwa.wheels[2].spin_axis_b = Vec3(sb * std::cos(5*M_PI/4),
                                            sb * std::sin(5*M_PI/4), cb).normalized();

    suite.rwa.wheels[3]             = rw_proto;
    suite.rwa.wheels[3].spin_axis_b = Vec3(sb * std::cos(7*M_PI/4),
                                            sb * std::sin(7*M_PI/4), cb).normalized();

    suite.mtq = MagnetorquerAssembly::make_xyz(10.0); // [A·m²]

    ActuatorCommands<4,3> cmds;
    cmds.tau_rwa_cmd = Vec3::Zero();
    cmds.tau_mtq_cmd = Vec3::Zero();
    cmds.thr_cmds    = {};

    Mat3 inertia_tensor;
    inertia_tensor << 18.5, -0.3,  0.1,
                    -0.3, 15.9, -0.5,
                    0.1, -0.5, 17.2;

    Mat3 inertia_tensor_inv = inertia_tensor.inverse();

    //timing parameters
    Real t0 = 0.0;
    Real tf = 3600; // simulate time in seconds
    Real dt = 0.1; // time step
    VecX t_vec = Eigen::VectorXd::LinSpaced(static_cast<int>((tf - t0) / dt) + 1, t0, tf);

    DerivFunc deriv = makeSixDoFDeriv(inertia_tensor, inertia_tensor_inv, force_config, torque_config, suite, cmds, q_d, ctrl, 0.1);

    SixDoFStateMat state_history = propagate6DoF(sixdof_init, t0, tf, dt, deriv);

    save_state_history(state_history);
    TranslationalStateMat trans_history = state_history.topRows(6);
    OEStateMat oe_history = translationalMatToOEMat(trans_history);
    save_oe_history(oe_history);

    // Helper lambdas
    auto print_header = [](const std::string& title) {
        std::cout << "\n╔══════════════════════════════════════╗\n";
        std::cout << "║  " << std::left << std::setw(36) << title << "║\n";
        std::cout << "╚══════════════════════════════════════╝\n";
    };

    auto print_field = [](const std::string& label, double value, const std::string& unit = "") {
        std::cout << "  " << std::left << std::setw(22) << label
                << ": " << std::fixed << std::setprecision(6) << value;
        if (!unit.empty()) std::cout << " " << unit;
        std::cout << "\n";
    };

    // Reusable print lambda — takes a 13-element source via index
    auto print_state = [&](auto get) {
        print_field("Pos X",  get(0),  "m");
        print_field("Pos Y",  get(1),  "m");
        print_field("Pos Z",  get(2),  "m");
        print_field("Vel X",  get(3),  "m/s");
        print_field("Vel Y",  get(4),  "m/s");
        print_field("Vel Z",  get(5),  "m/s");
        print_field("q.w",    get(6));
        print_field("q.x",    get(7));
        print_field("q.y",    get(8));
        print_field("q.z",    get(9));
        print_field("ω X",    get(10), "rad/s");
        print_field("ω Y",    get(11), "rad/s");
        print_field("ω Z",    get(12), "rad/s");
    };

    // Initial state — map struct fields to indices
    print_header("Initial State");
    print_state([&](int i) -> double {
        switch(i) {
            case 0:  return sixdof_init.eci.position.x();
            case 1:  return sixdof_init.eci.position.y();
            case 2:  return sixdof_init.eci.position.z();
            case 3:  return sixdof_init.eci.velocity.x();
            case 4:  return sixdof_init.eci.velocity.y();
            case 5:  return sixdof_init.eci.velocity.z();
            case 6:  return sixdof_init.attitude.q.w;
            case 7:  return sixdof_init.attitude.q.x;
            case 8:  return sixdof_init.attitude.q.y;
            case 9:  return sixdof_init.attitude.q.z;
            case 10: return sixdof_init.attitude.omega.x();
            case 11: return sixdof_init.attitude.omega.y();
            case 12: return sixdof_init.attitude.omega.z();
            default: return 0.0;
        }
    });

    // Final state — pull directly from Eigen column
    print_header("Final State");
    Eigen::VectorXd final_state = state_history.col(state_history.cols() - 1);
    print_state([&](int i) -> double { return final_state(i); });

    std::cout << "\n";

}