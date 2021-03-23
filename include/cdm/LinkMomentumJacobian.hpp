#pragma once

#include "MotionJacobian.hpp"
#include "utils.hpp"

template <typename MI, typename Tree>
Eigen::MatrixXd LinkMomentumJacobian(const std::string& bodyName, MI& info, Tree& tree)
{
    constexpr int ord = Tree::order;
    const auto& mb = info.model.mb;
    auto I = makeDiag<ord>(mb.body(mb.bodyIndexByName(bodyName)).inertia().matrix()); // \tilde{I} = diag(I)
    return I * MotionJacobian(bodyName, info, tree);
}

template <int JacOrder, typename MI, typename Tree>
Eigen::MatrixXd LinkMomentumJacobianOfOrder(const std::string& bodyName, MI& info, Tree& tree)
{
    constexpr int jac_ord = JacOrder;
    const auto& mb = info.model.mb;
    int j = mb.bodyIndexByName(bodyName);
    return mb.body(j).inertia().matrix() * MotionJacobianOfOrder<jac_ord>(bodyName, info, tree);
}