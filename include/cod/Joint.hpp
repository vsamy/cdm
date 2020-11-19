#pragma once

#include "cod/API.hpp"
#include "cod/typedefs.hpp"

namespace cod {

class COD_DLLAPI Joint {
public:
    enum class Type {
        Free,
        Spherical,
        Revolute,
        Prismatic,
        Fixed
    };

public:
    Joint() = default;
    Joint(const std::string& name, Type type, const Eigen::Vector3d& axis = Eigen::Vector3d::Zero());

    Type type() const noexcept;
    int nParam() const noexcept;
    int dof() const noexcept;
    const std::string& name() const noexcept;
    const MotionSubspace& S() const noexcept;

    friend bool operator==(const Joint& lhs, const Joint& rhs) { return lhs.m_name == rhs.m_name; }
    friend bool operator!=(const Joint& lhs, const Joint& rhs) { return !(lhs == rhs); }

private:
    void makeJoint(Type type, const Eigen::Vector3d& axis);

private:
    std::string m_name;
    Type m_type;
    MotionSubspace m_S;
    int m_params;
    int m_dof;
};

} // namespace cod
