#pragma once

#include "cdm/BasicJacobian.hpp"
#include "cdm/math_utility.hpp"

namespace cdm {

template <int Order>
Eigen::MatrixXd LinkMomentumJacobian(const Model& m, const ModelConfig<Order>& mc, const std::string& bodyName)
{
    auto I = makeDiag<Order>(m.body(m.bodyIndexByName(bodyName)).inertia().matrix()); // \tilde{I} = diag(I)
    return I * BasicJacobian<Order>(m, mc, bodyName);
}

template <int JacOrder, int Order>
Eigen::MatrixXd LinkMomentumJacobianOfOrder(const Model& m, const ModelConfig<Order>& mc, const std::string& bodyName)
{
    Index j = m.bodyIndexByName(bodyName);
    return m.body(j).inertia().matrix() * BasicJacobianOfOrder<JacOrder>(m, mc, bodyName);
}

} // namespace cdm
