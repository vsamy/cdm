#include "cdm/ModelConstructor.hpp"

namespace cdm {

void ModelConstructor::addLink(const Joint& joint, const Body& body)
{
    m_jointIndexFromName[joint.name()] = static_cast<Index>(m_jointBases.size());
    m_bodyIndexFromName[body.name()] = static_cast<Index>(m_bodyBases.size());
    m_jointBases.push_back(JointBase{ joint, -1 });
    m_bodyBases.push_back(LinkBase{ body, {} });
}

void ModelConstructor::connectLink(const std::string& linkName, const std::string& nextJointName, const Transform& T_l_j)
{
    assert(m_bodyIndexFromName.find(linkName) != m_bodyIndexFromName.end() && "linkName is unknown");
    assert(m_jointIndexFromName.find(nextJointName) != m_jointIndexFromName.end() && "nextJointName is unknown");
    Index bInd = m_bodyIndexFromName[linkName];
    Index jInd = m_jointIndexFromName[nextJointName];

    m_bodyBases[bInd].next.emplace_back(jInd, T_l_j);
    m_jointBases[jInd].previous = bInd;
    // TODO: Check for infinite loop that could occur due to recursion (only open kinematic chain for now)
}

Model ModelConstructor::build(const std::string& rootName, const Transform& A_world_root) const
{
    assert(m_jointIndexFromName.find(rootName) != m_jointIndexFromName.end() && "rootName is unknown");
    // TODO: warn if root is not first joint
    std::vector<Joint> joints;
    std::vector<Body> bodies;
    std::vector<Index> jointParents;
    std::vector<Index> jointChilds;
    std::vector<Transform> T0;
    std::unordered_map<std::string, Index> jointIndexByName;
    std::unordered_map<std::string, Index> bodyIndexByName;

    Index curPos = 0;
    Index curInd = m_jointIndexFromName.at(rootName);

    std::function<void(Index, const Transform&)> addNext;
    addNext = [&](Index linkIndex, const Transform& T) {
        const auto& jb = m_jointBases[linkIndex];
        const auto& bb = m_bodyBases[linkIndex];
        joints.push_back(jb.joint);
        bodies.push_back(bb.body);
        Index previousLink = jb.previous == -1 ? -1 : bodyIndexByName[m_bodyBases[jb.previous].body.name()];
        jointParents.push_back(previousLink);
        jointChilds.push_back(curPos);
        T0.push_back(T);
        jointIndexByName[jb.joint.name()] = curPos;
        bodyIndexByName[bb.body.name()] = curPos;
        curPos++;
        for (const auto& child : bb.next) {
            addNext(child.first, child.second);
        }
    };

    addNext(curInd, A_world_root);
    return { std::move(joints), std::move(bodies), std::move(jointParents), std::move(jointChilds), std::move(T0) };
}

} // namespace cdm
