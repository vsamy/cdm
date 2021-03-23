#include "cdm/Link.hpp"

namespace cdm {

Link::Link(const std::string& name, const Inertia& spatialInertia)
    : m_name(name)
    , m_inertia(spatialInertia)
{}

Link::Link(const std::string& name, double mass, const Eigen::Vector3d& momentum, const Eigen::Matrix3d& inertia)
    : m_name(name)
    , m_inertia(mass, momentum, inertia)
{}

const std::string& Link::name() const noexcept
{
    return m_name;
}

const Inertia& Link::inertia() const noexcept
{
    return m_inertia;
}

} // namespace cdm
