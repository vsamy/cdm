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

} // namespace cdm
