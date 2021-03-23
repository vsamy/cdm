#pragma once

#include "cdm/API.hpp"
#include "cdm/typedefs.hpp"

namespace cdm {

/*! \brief Representation of a link.
 * A link is composed of a name and a spatial inertia.
 */
class CDM_DLLAPI Link {
public:
    /*! \brief Default constructor. */
    Link() = default;
    /*! \brief Constructor.
     * \param name Name of the link.
     * \param spatialInertia Spatial inertia of the link at its origin.
     */
    Link(const std::string& name, const Inertia& spatialInertia);
    /*! \brief Constructor.
     * \param name Name of the link.
     * \param mass Mass of the link.
     * \param momentum Momentum of the link at its origin.
     * \param inertia Inertia of the link at its origin.
     */
    Link(const std::string& name, double mass, const Eigen::Vector3d& momentum, const Eigen::Matrix3d& inertia);

    /*! \brief Return name of link. */
    const std::string& name() const noexcept;
    /*! \brief Return spatial inertia of link. */
    const Inertia& inertia() const noexcept;

    friend bool operator==(const Link& lhs, const Link& rhs) { return lhs.m_name == rhs.m_name; }
    friend bool operator!=(const Link& lhs, const Link& rhs) { return !(lhs == rhs); }

private:
    std::string m_name; /*!< Name of the link */
    Inertia m_inertia; /*!< Spatial inertia of the link */
};

} // namespace cdm
