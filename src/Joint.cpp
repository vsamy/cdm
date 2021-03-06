/*
 * Copyright 2020-2021 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include "cdm/Joint.hpp"

namespace cdm {

Joint::Joint(const std::string& name, Type type, const Eigen::Vector3d& axis)
    : m_name(name)
    , m_type(type)
{
    makeJoint(type, axis);
}

Joint::Type Joint::type() const noexcept
{
    return m_type;
}

Index Joint::nParam() const noexcept
{
    return m_nParams;
}

Index Joint::dof() const noexcept
{
    return m_dof;
}

const std::string& Joint::name() const noexcept
{
    return m_name;
}

const MotionSubspace& Joint::S() const noexcept
{
    return m_S;
}

void Joint::makeJoint(Type type, const Eigen::Vector3d& axis)
{
    switch (type) {
    case Type::Free:
        m_S = MotionSubspace{ Eigen::Matrix6d::Identity() };
        m_nParams = 7;
        m_dof = 6;
        break;
    case Type::Spherical:
        m_S = MotionSubspace{ (Eigen::Matrix<double, 6, 3>() << Eigen::Matrix3d::Identity(), Eigen::Matrix3d::Zero()).finished() };
        m_nParams = 4;
        m_dof = 3;
        break;
    case Type::Revolute:
        m_S = MotionSubspace{ (Eigen::Vector6d() << axis, Eigen::Vector3d::Zero()).finished() };
        m_nParams = 1;
        m_dof = 1;
        break;
    case Type::Prismatic:
        m_S = MotionSubspace{ (Eigen::Vector6d() << Eigen::Vector3d::Zero(), axis).finished() };
        m_nParams = 1;
        m_dof = 1;
        break;
    default:
    case Type::Fixed:
        m_S = MotionSubspace{ Eigen::Matrix<double, 6, 0>{} };
        m_nParams = 0;
        m_dof = 0;
        break;
    }
}

} // namespace cdm
