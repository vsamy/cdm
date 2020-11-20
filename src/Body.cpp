#include "cdm/Body.hpp"

namespace cdm {

Body::Body(const std::string& name, const Inertia& spatialInertia)
    : m_name(name)
    , m_inertia(spatialInertia)
{}

Body::Body(const std::string& name, double mass, const Eigen::Vector3d& momentum, const Eigen::Matrix3d& inertia)
    : m_name(name)
    , m_inertia(mass, momentum, inertia)
{}

const std::string& Body::name() const noexcept
{
    return m_name;
}

const Inertia& Body::inertia() const noexcept
{
    return m_inertia;
}

} // namespace cdm
