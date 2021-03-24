#pragma once

#include "cdm/JointMomentumJacobian.hpp"

namespace cdm {

template <int Order>
Eigen::MatrixXd JointForceJacobian(const Model& m, const ModelConfig<Order>& mc, const std::string& bodyName)
{
    Index j = m.bodyIndexByName(bodyName);
    auto D = generateD(cdm::CrossN<Order>{ mc.bodyMotions[j].motion() });
    return D * JointMomentumJacobian(m, mc, bodyName);
}

template <int JacOrder, int Order>
Eigen::MatrixXd JointForceJacobianOfOrder(const Model& m, const ModelConfig<Order>& mc, const std::string& bodyName)
{
    if constexpr (JacOrder == 0) {
        return JointMomentumJacobianOfOrder<JacOrder>(m, mc, bodyName);
    } else {
        Index j = m.bodyIndexByName(bodyName);
        CrossN<JacOrder> cx{ mc.bodyMotions[j].motion() };
        Eigen::MatrixXd D = cx.matrix().middleRows<6 * JacOrder>(6 * (JacOrder - 1)) / static_cast<double>(JacOrder); // TODO: Should be Scalar
        D.middleCols<6>(6 * JacOrder).setIdentity();
        return D * JointMomentumJacobianOfOrder<JacOrder>(m, mc, bodyName);
    }
}

} // namespace cdm
