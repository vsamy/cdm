#include "cdm/ModelConstructor.hpp"

namespace cdm {

void ModelConstructor::addLink(const Joint& joint, const Body& body)
{
    m_jointIndexFromName[joint.name()] = static_cast<Index>(m_jointBases.size());
    m_bodyIndexFromName[body.name()] = static_cast<Index>(m_bodyBases.size());
    m_jointBases.push_back(JointBase{ joint, -1 });
    m_bodyBases.push_back(BodyBase{ body, {} });
}

void ModelConstructor::link(const std::string& bodyName, const std::string& nextJointName, const Transform& A_b_j)
{
    assert(m_bodyIndexFromName.find(bodyName) != m_bodyIndexFromName.end() && "bodyName is unknown");
    assert(m_jointIndexFromName.find(nextJointName) != m_jointIndexFromName.end() && "nextJointName is unknown");
    Index bInd = m_bodyIndexFromName[bodyName];
    Index jInd = m_jointIndexFromName[nextJointName];

    m_bodyBases[bInd].next.emplace_back(jInd, A_b_j);
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
    std::vector<Transform> A0;
    std::unordered_map<std::string, Index> jointIndexByName;
    std::unordered_map<std::string, Index> bodyIndexByName;

    Index curPos = 0;
    Index curInd = m_jointIndexFromName.at(rootName);

    std::function<void(Index, const Transform&)> addNext;
    addNext = [&](Index linkIndex, const Transform& A) {
        const auto& jb = m_jointBases[linkIndex];
        const auto& bb = m_bodyBases[linkIndex];
        joints.push_back(jb.joint);
        bodies.push_back(bb.body);
        Index previousBody = jb.previous == -1 ? -1 : bodyIndexByName[m_bodyBases[jb.previous].body.name()];
        jointParents.push_back(previousBody);
        jointChilds.push_back(curPos);
        A0.push_back(A);
        jointIndexByName[jb.joint.name()] = curPos;
        bodyIndexByName[bb.body.name()] = curPos;
        curPos++;
        for (const auto& child : bb.next) {
            addNext(child.first, child.second);
        }
    };

    addNext(curInd, A_world_root);
    return { std::move(joints), std::move(bodies), std::move(jointParents), std::move(jointChilds), std::move(A0) };
}

} // namespace cdm
