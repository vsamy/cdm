#pragma once

#include "LinkMomentumJacobian.hpp"

template <typename MI, typename Tree>
Eigen::MatrixXd LinkForceJacobian(const std::string& bodyName, MI& info, Tree& tree)
{
    constexpr int ord = Tree::order;
    const auto& mb = info.model.mb;
    int j = mb.bodyIndexByName(bodyName);
    auto D = generateD(coma::CrossNd<ord>{ tree.bodies[j].motion() });
    return D * LinkMomentumJacobian(bodyName, info, tree);
}

template <int JacOrder, typename MI, typename Tree>
Eigen::MatrixXd LinkForceJacobianOfOrder(const std::string& bodyName, MI& info, Tree& tree)
{
    constexpr int jac_ord = JacOrder;
    if constexpr (jac_ord == 0) {
        return LinkMomentumJacobianOfOrder<jac_ord>(bodyName, info, tree);
    } else {
        const auto& mb = info.model.mb;
        int j = mb.bodyIndexByName(bodyName);
        coma::CrossNd<jac_ord> cx{ tree.bodies[j].motions() };
        Eigen::MatrixXd D = cx.matrix().middleRows<6 * jac_ord>(6 * (jac_ord - 1)) / jac_ord;
        D.middleCols<6>(6 * jac_ord).setIdentity();
        return D * LinkMomentumJacobianOfOrder<jac_ord>(bodyName, info, tree);
    }
}