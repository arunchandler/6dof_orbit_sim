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

    Quaternion operator*(const Quaternion& rhs) const noexcept {
        return {
            w*rhs.w - x*rhs.x - y*rhs.y - z*rhs.z,
            w*rhs.x + x*rhs.w + y*rhs.z - z*rhs.y,
            w*rhs.y - x*rhs.z + y*rhs.w + z*rhs.x,
            w*rhs.z + x*rhs.y - y*rhs.x + z*rhs.w
        };
    }

    Quaternion conjugate() const noexcept { return {w, -x, -y, -z}; }

    Vec3 rotate(const Vec3& v) const noexcept {
        Quaternion v_quat{0.0, v.x(), v.y(), v.z()};
        Quaternion result = (*this) * v_quat * this->conjugate();
        return Vec3(result.x, result.y, result.z);
    }

    Real norm() const noexcept { return std::sqrt(w*w + x*x + y*y + z*z); }

    Quaternion operator*(Real s) const noexcept { return {w*s, x*s, y*s, z*s}; }
    friend Quaternion operator*(Real s, const Quaternion& q) noexcept { return q * s; }

};

struct QuaternionDot {
    Real w, x, y, z;

    // Only operations valid on a derivative:
    QuaternionDot operator+(const QuaternionDot& rhs) const noexcept {
        return {w+rhs.w, x+rhs.x, y+rhs.y, z+rhs.z};
    }

    QuaternionDot& operator+=(const QuaternionDot& rhs) noexcept {
        w += rhs.w;
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    QuaternionDot operator*(Real s) const noexcept { return {w*s, x*s, y*s, z*s}; }
    friend QuaternionDot operator*(Real s, const QuaternionDot& qd) noexcept { return qd * s; }

};

struct AttitudeState {
    Quaternion q; // Attitude quaternion
    Vec3 omega;   // Angular velocity in body frame [rad/s]
};

struct AttitudeStateDot {
    QuaternionDot q_dot; // Time derivative of quaternion
    Vec3 omega_dot;   // Time derivative of angular velocity [rad/s²]
};

inline QuaternionDot quaternion_kinematics(const Quaternion& q, const Vec3& omega) noexcept {
    return {
        0.5 * (-q.x*omega.x() - q.y*omega.y() - q.z*omega.z()),
        0.5 * ( q.w*omega.x() + q.y*omega.z() - q.z*omega.y()),
        0.5 * ( q.w*omega.y() - q.x*omega.z() + q.z*omega.x()),
        0.5 * ( q.w*omega.z() + q.x*omega.y() - q.y*omega.x())
    };
}

} // namespace orb
