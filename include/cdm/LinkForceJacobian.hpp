#pragma once

#include "cdm/LinkMomentumJacobian.hpp"

namespace cdm {

template <int Order>
Eigen::MatrixXd LinkForceJacobian(const Model& m, const ModelConfig<Order>& mc, const std::string& bodyName)
{
    int j = m.bodyIndexByName(bodyName);
    auto D = generateD(CrossNd<Order>{ mc.bodyMotions[j].motion() });
    return D * LinkMomentumJacobian(m, mc, bodyName);
}

template <int Order, int JacOrder>
Eigen::MatrixXd LinkForceJacobianOfOrder(const Model& m, const ModelConfig<Order>& mc, const std::string& bodyName)
{
    if constexpr (JacOrder == 0) {
        return LinkMomentumJacobianOfOrder<JacOrder>(m, mc, bodyName);
    } else {
        int j = m.bodyIndexByName(bodyName);
        CrossNd<JacOrder> cx{ mc.bodyMotions[j].motions() };
        Eigen::MatrixXd D = cx.matrix().middleRows<6 * JacOrder>(6 * (JacOrder - 1)) / JacOrder;
        D.middleCols<6>(6 * JacOrder).setIdentity();
        return D * LinkMomentumJacobianOfOrder<JacOrder>(m, mc, bodyName);
    }
}

} // namespace cdm
