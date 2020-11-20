#pragma once

#include "cdm/API.hpp"
#include "cdm/typedefs.hpp"

namespace cdm {

class CDM_DLLAPI Body {
public:
    Body() = default;
    Body(const std::string& name, const Inertia& spatialInertia);
    Body(const std::string& name, double mass, const Eigen::Vector3d& momentum, const Eigen::Matrix3d& inertia);

    const std::string& name() const noexcept;
    const Inertia& inertia() const noexcept;

    friend bool operator==(const Body& lhs, const Body& rhs) { return lhs.m_name == rhs.m_name; }
    friend bool operator!=(const Body& lhs, const Body& rhs) { return !(lhs == rhs); }

private:
    std::string m_name;
    Inertia m_inertia;
};

} // namespace cdm
