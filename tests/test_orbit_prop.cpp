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
    torque_config.useGravityGradient = false;
    torque_config.useDrag = true;
    torque_config.useSRP = false;
    torque_config.center_of_mass = Vec3(0.1, 0.05, 0.0);

    ForceModelConfig force_config;
    force_config.useDrag = true;
    force_config.useJ2 = true;
    force_config.useSRP = true;

    Mat3 inertia_tensor;
    inertia_tensor << 18.5, -0.3,  0.1,
                    -0.3, 15.9, -0.5,
                    0.1, -0.5, 17.2;

    Mat3 inertia_tensor_inv = inertia_tensor.inverse();

    DerivFunc deriv = makeSixDoFDeriv(inertia_tensor, inertia_tensor_inv, force_config, torque_config);

    //timing parameters
    Real t0 = 0.0;
    Real tf = 3600; // simulate for one hour
    Real dt = 0.1; // time step of 0.1 seconds
    VecX t_vec = Eigen::VectorXd::LinSpaced(static_cast<int>((tf - t0) / dt) + 1, t0, tf);

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