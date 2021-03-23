#pragma once

#include "cdm/Model.hpp"
#include "cdm/ModelConfig.hpp"

namespace cdm {

template <int Order>
Eigen::MatrixXd BasicJacobian(const Model& m, const ModelConfig<Order>& mc, const std::string& bodyName)
{
    static_assert(Order >= 0, "Not yet ready for Dynamic");
    const auto& parents = m.parents();
    const auto& posInDof = m.jointsPosInDof();
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(6 * order, m.nDof() * order);
    Index j = m.bodyIndexByName(bodyName);
    auto C_b_0 = mc.linkMotions[j].inverse();
    while (j != -1) {
        J.block(0, order * posInDof[j], 6 * order, m.joint(j).dof()) = ((C_b_0 * mc.linkMotions[j]) * DiMotionSupace{ m.joint(j).S() }).matrix();
        j = parents[j];
    }

    return J;
}

template <int Order, int JacOrder>
Eigen::MatrixXd BasicJacobianOfOrder(const Model& m, const ModelConfig<Order>& mc, const std::string& bodyName)
{
    return BasicJacobianOfOrder(m, mbc, bodyName, JacOrder);
}

template <int Order>
Eigen::MatrixXd BasicJacobianOfOrder(const Model& m, const ModelConfig<Order>& mc, const std::string& bodyName, int jacOrder)
{
    const auto& parents = m.parents();
    const auto& posInDof = m.jointsPosInDof();
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(6, m.nDof());
    Index j = m.bodyIndexByName(bodyName);
    auto C_b_0 = mc.linkMotions[j].inverse();
    while (j != -1) {
        // A faster code would recursively compute
        auto C_b_j = C_b_0 * mc.linkMotions[j];
        auto S = mb.joint(j).S();
        if (jacOrder == 0) {
            J.block(0, posInDof[j], 6, S.cols()) = C_b_j.transform() * S.matrix();
        } else {
            J.block(0, posInDof[j], 6, S.cols()) = (C_b_j[jacOrder - 1] * S).matrix();
        }
        j = parents[j];
    }

    return J;
}

} // namespace cdm
