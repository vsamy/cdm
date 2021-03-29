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
    , m_jointPosInParam(static_cast<size_t>(m_nLinks))
    , m_jointPosInDof(static_cast<size_t>(m_nLinks))
{
    for (Index i = 0; i < m_nLinks; ++i) {
        size_t ui = static_cast<size_t>(i);
        m_bodyIndexByName[m_bodies[ui].name()] = i;
        m_jointIndexByName[m_joints[ui].name()] = i;
        m_jointPosInParam[ui] = m_nParam;
        m_jointPosInDof[ui] = m_nDof;
        m_nParam += m_joints[ui].nParam();
        m_nDof += m_joints[ui].dof();
    }
}

} // namespace cdm
