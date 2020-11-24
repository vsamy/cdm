#pragma once

#include "cdm/Body.hpp"
#include "cdm/Joint.hpp"

namespace cdm {

class CDM_DLLAPI Model {
public:
    Model() = default;
    Model(std::vector<Joint> joints, std::vector<Body> bodies, std::vector<Index> jointParents,
        std::vector<Index> jointChilds, std::vector<Transform> A0);

    Index nLinks() const noexcept;
    Index nParam() const noexcept;
    Index nDof() const noexcept;
    Index jointParent(Index jointIndex) const noexcept;
    Index jointParentAt(Index jointIndex) const;
    Index jointChild(Index jointIndex) const noexcept;
    Index jointChildAt(Index jointIndex) const;
    Index jointPosInParam(Index jointIndex) const;
    Index jointPosInDof(Index jointIndex) const;

    const Joint& joint(Index jointIndex) const noexcept;
    const Joint& jointAt(Index jointIndex) const;
    const Body& body(Index bodyIndex) const noexcept;
    const Body& bodyAt(Index bodyIndex) const;
    const Transform& A0(Index linkIndex) const noexcept;
    const Transform& A0At(Index linkIndex) const;

    const std::vector<Joint>& joints() const noexcept;
    const std::vector<Body>& bodies() const noexcept;
    const std::vector<Index>& jointParents() const noexcept;
    const std::vector<Index>& jointChilds() const noexcept;
    const std::vector<Transform>& A0s() const noexcept;

private:
    Index m_nLinks;
    Index m_nParam;
    Index m_nDof;

    std::vector<Joint> m_joints;
    std::vector<Body> m_bodies;
    std::vector<Index> m_jointParents;
    std::vector<Index> m_jointChilds;
    std::vector<Transform> m_A0;
    std::vector<Index> m_jointPosInParam;
    std::vector<Index> m_jointPosInDof;
    std::unordered_map<std::string, Index> m_jointIndexByName;
    std::unordered_map<std::string, Index> m_bodyIndexByName;
};

} // namespace cdm
