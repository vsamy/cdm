#pragma once

#include "cdm/Body.hpp"
#include "cdm/Joint.hpp"

namespace cdm {

class CDM_DLLAPI Model {
public:
    Model() = default;
    Model(std::vector<Joint> joints, std::vector<Body> bodies, std::vector<int> jointParents,
        std::vector<int> jointChilds, std::vector<Transform> A0);

    int nLinks() const noexcept;
    int nParam() const noexcept;
    int nDof() const noexcept;
    int jointParent(int jointIndex) const noexcept;
    int jointParentAt(int jointIndex) const;
    int jointChild(int jointIndex) const noexcept;
    int jointChildAt(int jointIndex) const;
    int jointPosInParam(int jointIndex) const;
    int jointPosInDof(int jointIndex) const;

    const Joint& joint(int jointIndex) const noexcept;
    const Joint& jointAt(int jointIndex) const;
    const Body& body(int bodyIndex) const noexcept;
    const Body& bodyAt(int bodyIndex) const;
    const Transform& A0(int linkIndex) const noexcept;
    const Transform& A0At(int linkIndex) const;

    const std::vector<Joint>& joints() const noexcept;
    const std::vector<Body>& bodies() const noexcept;
    const std::vector<int>& jointParents() const noexcept;
    const std::vector<int>& jointChilds() const noexcept;
    const std::vector<Transform>& A0s() const noexcept;

private:
    int m_nLinks;
    int m_nParam;
    int m_nDof;

    std::vector<Joint> m_joints;
    std::vector<Body> m_bodies;
    std::vector<int> m_jointParents;
    std::vector<int> m_jointChilds;
    std::vector<Transform> m_A0;
    std::vector<int> m_jointPosInParam;
    std::vector<int> m_jointPosInDof;
    std::unordered_map<std::string, int> m_jointIndexByName;
    std::unordered_map<std::string, int> m_bodyIndexByName;
};

} // namespace cdm
