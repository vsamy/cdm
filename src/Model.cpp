#include "cdm/Model.hpp"

namespace cdm {

Model::Model(std::vector<Joint> joints, std::vector<Body> bodies, std::vector<Index> jointParents,
    std::vector<Index> jointChilds, std::vector<Transform> A0)
    : m_nLinks(static_cast<Index>(joints.size()))
    , m_nParam(0)
    , m_nDof(0)
    , m_joints(std::move(joints))
    , m_bodies(std::move(bodies))
    , m_jointParents(std::move(jointParents))
    , m_jointChilds(std::move(jointChilds))
    , m_A0(std::move(A0))
    , m_jointPosInParam(m_nLinks)
    , m_jointPosInDof(m_nLinks)
{
    for (Index i = 0; i < m_nLinks; ++i) {
        m_bodyIndexByName[m_bodies[i].name()] = i;
        m_jointIndexByName[m_joints[i].name()] = i;
        m_jointPosInParam[i] = m_nParam;
        m_jointPosInDof[i] = m_nDof;
        m_nParam += m_joints[i].nParam();
        m_nDof += m_joints[i].dof();
    }
}

Index Model::nLinks() const noexcept
{
    return m_nLinks;
}

Index Model::nParam() const noexcept
{
    return m_nParam;
}

Index Model::nDof() const noexcept
{
    return m_nDof;
}

Index Model::jointParent(Index jointIndex) const noexcept
{
    return m_jointParents[jointIndex];
}

Index Model::jointParentAt(Index jointIndex) const
{
    return m_jointParents.at(jointIndex);
}

Index Model::jointChild(Index jointIndex) const noexcept
{
    return m_jointChilds[jointIndex];
}

Index Model::jointChildAt(Index jointIndex) const
{
    return m_jointChilds.at(jointIndex);
}

Index Model::jointPosInParam(Index jointIndex) const
{
    return m_jointPosInParam.at(jointIndex);
}

Index Model::jointPosInDof(Index jointIndex) const
{
    return m_jointPosInDof.at(jointIndex);
}

const Joint& Model::joint(Index jointIndex) const noexcept
{
    return m_joints[jointIndex];
}

const Joint& Model::jointAt(Index jointIndex) const
{
    return m_joints.at(jointIndex);
}

const Body& Model::body(Index bodyIndex) const noexcept
{
    return m_bodies[bodyIndex];
}

const Body& Model::bodyAt(Index bodyIndex) const
{
    return m_bodies.at(bodyIndex);
}

const Transform& Model::A0(Index linkIndex) const noexcept
{
    return m_A0[linkIndex];
}

const Transform& Model::A0At(Index linkIndex) const
{
    return m_A0.at(linkIndex);
}

const std::vector<Joint>& Model::joints() const noexcept
{
    return m_joints;
}

const std::vector<Body>& Model::bodies() const noexcept
{
    return m_bodies;
}

const std::vector<Index>& Model::jointParents() const noexcept
{
    return m_jointParents;
}

const std::vector<Index>& Model::jointChilds() const noexcept
{
    return m_jointChilds;
}

const std::vector<Transform>& Model::A0s() const noexcept
{
    return m_A0;
}

} // namespace cdm
