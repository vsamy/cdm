#pragma once

namespace cdm {

template <int Order>
struct ModelConfig {
    CMTM<Order> world;
    Eigen::VectorXd q;
    Eigen::MatrixXd dqs;
    std::vector<CMTM<Order>> jointMotions;
    std::vector<CMTM<Order>> bodyMotions;
    std::vector<ForceVectorX<Order>> jointMomentums;
    std::vector<ForceVectorX<Order>> bodyMomentums;
    std::vector<ForceVectorX<Order>> jointForces;
    std::vector<ForceVectorX<Order>> bodyForces;
    std::vector<Eigen::VectorXd> jointTorques;
};

} // namespace cdm
