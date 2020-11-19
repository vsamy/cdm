#pragma once

#include "cod/API.hpp"
#include "cod/typedefs.hpp"

namespace cod {

struct COD_DLLAPI ModelConfig {
    int order;
    CMTM world;
    Eigen::VectorXd q;
    Eigen::MatrixXd dqs;
    std::vector<CMTM> jointMotions;
    std::vector<CMTM> bodyMotions;
    std::vector<ForceVectorX> jointMomentums;
    std::vector<ForceVectorX> bodyMomentums;
    std::vector<ForceVectorX> jointForces;
    std::vector<ForceVectorX> bodyForces;
    std::vector<Eigen::VectorXd> jointTorques;
};

} // namespace cod
