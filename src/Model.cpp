#include "cdm/Model.hpp"

namespace cdm {

Model::Model(std::vector<Joint> joints, std::vector<Body> bodies, std::vector<int> jointParents,
    std::vector<int> jointChilds, std::vector<Transform> A0)
    : m_nLinks(static_cast<int>(joints.size()))
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
    for (int i = 0; i < m_nLinks; ++i) {
        m_bodyIndexByName[m_bodies[i].name()] = i;
        m_jointIndexByName[m_joints[i].name()] = i;
        m_jointPosInParam[i] = m_nParam;
        m_jointPosInDof[i] = m_nDof;
        m_nParam += m_joints[i].nParam();
        m_nDof += m_joints[i].dof();
    }
}

int Model::nLinks() const noexcept
{
    return m_nLinks;
}

int Model::nParam() const noexcept
{
    return m_nParam;
}

int Model::nDof() const noexcept
{
    return m_nDof;
}

int Model::jointParent(int jointIndex) const noexcept
{
    return m_jointParents[jointIndex];
}

int Model::jointParentAt(int jointIndex) const
{
    return m_jointParents.at(jointIndex);
}

int Model::jointChild(int jointIndex) const noexcept
{
    return m_jointChilds[jointIndex];
}

int Model::jointChildAt(int jointIndex) const
{
    return m_jointChilds.at(jointIndex);
}

int Model::jointPosInParam(int jointIndex) const
{
    return m_jointPosInParam.at(jointIndex);
}

int Model::jointPosInDof(int jointIndex) const
{
    return m_jointPosInDof.at(jointIndex);
}

const Joint& Model::joint(int jointIndex) const noexcept
{
    return m_joints[jointIndex];
}

const Joint& Model::jointAt(int jointIndex) const
{
    return m_joints.at(jointIndex);
}

const Body& Model::body(int bodyIndex) const noexcept
{
    return m_bodies[bodyIndex];
}

const Body& Model::bodyAt(int bodyIndex) const
{
    return m_bodies.at(bodyIndex);
}

const Transform& Model::A0(int linkIndex) const noexcept
{
    return m_A0[linkIndex];
}

const Transform& Model::A0At(int linkIndex) const
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

const std::vector<int>& Model::jointParents() const noexcept
{
    return m_jointParents;
}

const std::vector<int>& Model::jointChilds() const noexcept
{
    return m_jointChilds;
}

const std::vector<Transform>& Model::A0s() const noexcept
{
    return m_A0;
}

} // namespace cdm
