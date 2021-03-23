#pragma once

#include "cdm/API.hpp"
#include "cdm/typedefs.hpp"

namespace cdm {

/*! \brief Representation of a model joint.
 * Only basic type joints are available. \see Joint::Type
 */
class CDM_DLLAPI Joint {
public:
    /*! \brief Joint types */
    enum class Type {
        Free, /*!< Free flyer joint */
        Spherical, /*!< Spherical or ball joint */
        Revolute, /*!< Revolute or hinge joint */
        Prismatic, /*!< Prismatic or linear joint */
        Fixed /*!< Fixed joint */
    };

public:
    /*! \brief Default constructor. */
    Joint() = default;
    /*! \brief Joint constructor.
     * \param name Name of the joint.
     * \param type Type of the joint.
     * \param axis Axis of the joint. Used for type Revolute or Prismatic.
     */
    Joint(const std::string& name, Type type, const Eigen::Vector3d& axis = Eigen::Vector3d::Zero());

    /*! \brief Get type of the joint. */
    Type type() const noexcept;
    /*! \brief Get number of parameter of the joint. */
    Index nParam() const noexcept;
    /*! \brief Get number of DoF of the joint. */
    Index dof() const noexcept;
    /*! \brief Get the name of the joint. */
    const std::string& name() const noexcept;
    /*! \brief Get the motion subspace of the joint. */
    const MotionSubspace& S() const noexcept;

    friend bool operator==(const Joint& lhs, const Joint& rhs) { return lhs.m_name == rhs.m_name; }
    friend bool operator!=(const Joint& lhs, const Joint& rhs) { return !(lhs == rhs); }

private:
    void makeJoint(Type type, const Eigen::Vector3d& axis);

private:
    std::string m_name; /*!< Name of the joint */
    Type m_type; /*!< Type of the joint */
    MotionSubspace m_S; /*!< Motion subspace of the joint */
    Index m_nParams; /*!< Number of joint parameters */
    Index m_dof; /*!< Number of joint DoF */
};

} // namespace cdm
