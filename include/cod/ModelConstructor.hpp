#pragma once

#include "cod/API.hpp"
#include "cod/Model.hpp"
#include <unordered_map>
#include <vector>

namespace cod {

class COD_DLLAPI ModelConstructor {
public:
    ModelConstructor() = default;

    void addLink(const Joint& joint, const Body& body);
    void link(const std::string& bodyName, const std::string& nextJointName, const Transform& A_b_j);
    Model build(const std::string& rootName, const Transform& A_world_root = Transform::Identity()) const;

private:
    struct BodyBase {
        Body body;
        std::vector<std::pair<int, Transform>> next;
    };
    struct JointBase {
        Joint joint;
        int previous;
    };

private:
    std::vector<BodyBase> m_bodyBases;
    std::vector<JointBase> m_jointBases;
    std::unordered_map<std::string, int> m_jointIndexFromName;
    std::unordered_map<std::string, int> m_bodyIndexFromName;
};

} // namespace cod
