#pragma once

#include "cdm/Model.hpp"

namespace cdm {

class CDM_DLLAPI ModelConstructor {
public:
    ModelConstructor() = default;

    void addLink(const Joint& joint, const Body& body);
    void link(const std::string& bodyName, const std::string& nextJointName, const Transform& A_b_j);
    Model build(const std::string& rootName, const Transform& A_world_root = Transform::Identity()) const;

private:
    struct BodyBase {
        Body body;
        std::vector<std::pair<Index, Transform>> next;
    };
    struct JointBase {
        Joint joint;
        Index previous;
    };

private:
    std::vector<BodyBase> m_bodyBases;
    std::vector<JointBase> m_jointBases;
    std::unordered_map<std::string, Index> m_jointIndexFromName;
    std::unordered_map<std::string, Index> m_bodyIndexFromName;
};

} // namespace cdm
