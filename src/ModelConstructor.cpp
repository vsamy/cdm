#include "cdm/ModelConstructor.hpp"

namespace cdm {

void ModelConstructor::addPair(const Joint& joint, const Link& link)
{
    m_jointIndexFromName[joint.name()] = static_cast<Index>(m_jointBases.size());
    m_linkIndexFromName[link.name()] = static_cast<Index>(m_linkBases.size());
    m_jointBases.push_back(JointBase{ joint, -1 });
    m_linkBases.push_back(LinkBase{ link, {} });
}

void ModelConstructor::connect(const std::string& linkName, const std::string& nextJointName, const Transform& T_l_j)
{
    assert(m_linkIndexFromName.find(linkName) != m_linkIndexFromName.end() && "linkName is unknown");
    assert(m_jointIndexFromName.find(nextJointName) != m_jointIndexFromName.end() && "nextJointName is unknown");
    Index bInd = m_linkIndexFromName[linkName];
    Index jInd = m_jointIndexFromName[nextJointName];

    m_linkBases[bInd].next.emplace_back(jInd, T_l_j);
    m_jointBases[jInd].previous = bInd;
    // TODO: Check for infinite loop that could occur due to recursion (only open kinematic chain for now)
}

Model ModelConstructor::build(const std::string& rootName, const Transform& A_world_root) const
{
    assert(m_jointIndexFromName.find(rootName) != m_jointIndexFromName.end() && "rootName is unknown");
    // TODO: warn if root is not first joint
    std::vector<Joint> joints;
    std::vector<Link> links;
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
        const auto& bb = m_linkBases[linkIndex];
        joints.push_back(jb.joint);
        links.push_back(bb.link);
        Index previousLink = jb.previous == -1 ? -1 : bodyIndexByName[m_linkBases[jb.previous].link.name()];
        jointParents.push_back(previousLink);
        jointChilds.push_back(curPos);
        T0.push_back(T);
        jointIndexByName[jb.joint.name()] = curPos;
        bodyIndexByName[bb.link.name()] = curPos;
        curPos++;
        for (const auto& child : bb.next) {
            addNext(child.first, child.second);
        }
    };

    addNext(curInd, A_world_root);
    return { std::move(joints), std::move(links), std::move(jointParents), std::move(jointChilds), std::move(T0) };
}

} // namespace cdm
