#pragma once
#include <array>
#include <cmath>
#include "types.hpp"

namespace orb {

struct Quaternion {
    Real w, x, y, z;

    // Normalize the quaternion to unit length
    void normalize() {
        Real norm = std::sqrt(w*w + x*x + y*y + z*z);
        if (norm > 0.0) {
            w /= norm;
            x /= norm;
            y /= norm;
            z /= norm;
        }
    }

    Quaternion operator*(const Quaternion& rhs) const noexcept;

    Quaternion conjugate() const noexcept { return {w, -x, -y, -z}; }

    Vec3 rotate(const Vec3& v) const noexcept {
        Quaternion v_quat{0.0, v.x(), v.y(), v.z()};
        Quaternion result = (*this) * v_quat * this->conjugate();
        return Vec3(result.x, result.y, result.z);
    }

    Real norm() const noexcept { return std::sqrt(w*w + x*x + y*y + z*z); }

};

struct AttitudeState {
    Quaternion q; // Attitude quaternion
    Vec3 omega;   // Angular velocity in body frame [rad/s]
};

struct AttitudeStateDot {
    Quaternion q_dot; // Time derivative of quaternion
    Vec3 omega_dot;   // Time derivative of angular velocity [rad/s²]
};

} // namespace orb
