#pragma once

#include "cdm/API.hpp"
#include "cdm/typedefs.hpp"

namespace cdm {

/*! \brief Representation of a body.
 * A body is composed of a name and a spatial inertia.
 */
class CDM_DLLAPI Body {
public:
    /*! \brief Default constructor. */
    Body() = default;
    /*! \brief Constructor.
     * \param name Name of the body.
     * \param spatialInertia Spatial inertia of the body at its origin.
     */
    Body(const std::string& name, const Inertia& spatialInertia);
    /*! \brief Constructor.
     * \param name Name of the body.
     * \param mass Mass of the body.
     * \param momentum Momentum of the body at its origin.
     * \param inertia Inertia of the body at its origin.
     */
    Body(const std::string& name, double mass, const Eigen::Vector3d& momentum, const Eigen::Matrix3d& inertia);

    /*! \brief Return name of body. */
    const std::string& name() const noexcept;
    /*! \brief Return spatial inertia of body. */
    const Inertia& inertia() const noexcept;

    friend bool operator==(const Body& lhs, const Body& rhs) { return lhs.m_name == rhs.m_name; }
    friend bool operator!=(const Body& lhs, const Body& rhs) { return !(lhs == rhs); }

private:
    std::string m_name; /*!< Name of the body */
    Inertia m_inertia; /*!< Spatial inertia of the body */
};

} // namespace cdm
