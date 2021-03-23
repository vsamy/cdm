#include "cdm/Model.hpp"

namespace cdm {

Model::Model(std::vector<Joint> joints, std::vector<Body> bodies, std::vector<Index> jointParents,
    std::vector<Index> jointChilds, std::vector<Transform> T0)
    : m_nLinks(static_cast<Index>(joints.size()))
    , m_nParam(0)
    , m_nDof(0)
    , m_joints(std::move(joints))
    , m_bodies(std::move(bodies))
    , m_jointParents(std::move(jointParents))
    , m_jointChilds(std::move(jointChilds))
    , m_T0(std::move(T0))
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

const Body& Model::body(Index linkIndex) const noexcept
{
    return m_bodies[linkIndex];
}

const Body& Model::bodyAt(Index linkIndex) const
{
    return m_bodies.at(linkIndex);
}

const Transform& Model::T0(Index jointIndex) const noexcept
{
    return m_T0[jointIndex];
}

const Transform& Model::T0At(Index jointIndex) const
{
    return m_T0.at(jointIndex);
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

const std::vector<Index>& Model::jointPosInParam() const noexcept
{
    return m_jointPosInParam;
}

const std::vector<Index>& Model::jointPosInDof() const noexcept
{
    return m_jointPosInDof;
}

const std::vector<Transform>& Model::T0s() const noexcept
{
    return m_T0;
}

} // namespace cdm
