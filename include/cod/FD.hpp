#pragma once

#include "cod/API.hpp"
#include "cod/Model.hpp"
#include "cod/ModelConfig.hpp"

namespace cod {

COD_DLLAPI Eigen::VectorXd FD(const Model& m, const ModelConfig& mc, const std::vector<Eigen::VectorXd>& Tau);
// COD_DLLAPI Eigen::VectorXd standardFD(const Model& m, ModelConfig& mc, const std::vector<Eigen::VectorXd>& Tau);

} // namespace cod
