// include/6dof_orbit_sim/plot.hpp
#pragma once
#include <fstream>
#include <Eigen/Dense>

inline void save_state_history(const Eigen::MatrixXd& state_history, 
                                const std::string& filename = "output/state_history.csv") {
    std::ofstream file(filename);
    file << "pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,q_w,q_x,q_y,q_z,wx,wy,wz\n";
    for (int i = 0; i < state_history.cols(); ++i) {
        for (int j = 0; j < state_history.rows(); ++j) {
            file << state_history(j, i);
            if (j < state_history.rows() - 1) file << ",";
        }
        file << "\n";
    }
}

inline void save_oe_history(const Eigen::MatrixXd& oe_history, 
                                const std::string& filename = "output/oe_history.csv") {
    std::ofstream file(filename);
    file << "a,e,i,RAAN,argp,nu\n";
    for (int i = 0; i < oe_history.cols(); ++i) {
        for (int j = 0; j < oe_history.rows(); ++j) {
            file << oe_history(j, i);
            if (j < oe_history.rows() - 1) file << ",";
        }        file << "\n";
    }
}