#pragma once

#include "cdm/API.hpp"
#include "cdm/typedefs.hpp"

namespace cdm {

class CDM_DLLAPI Joint {
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
    Index nParam() const noexcept;
    Index dof() const noexcept;
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
    Index m_params;
    Index m_dof;
};

} // namespace cdm
