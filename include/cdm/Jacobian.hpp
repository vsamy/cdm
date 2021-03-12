#pragma once

#include "cdm/Model.hpp"
#include "cdm/ModelConfig.hpp"

namespace cdm {

template <int Order>
Eigen::MatrixXd MotionJacobian(const Model& m, const ModelConfig<Order>& mc, const std::string& bodyName)
{
    const int order = mc.world.order();
    const auto& parents = m.parents();
    const auto& pos = m.jointsPosInDof();
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(6 * order, m.nDof() * order);
    Index j = m.bodyIndexByName(bodyName);
    auto C_b_0 = mc.linkMotions[j].inverse();
    while (j != -1) {
        J.block(0, order * pos[j], 6 * order, m.joint(j).dof()) = ((C_b_0 * mc.linkMotions[j]) * DiMotionSupace{ m.joint(j).S() }).matrix();
        j = parents[j];
    }

    return J;
}

template <int JacOrder, int Order>
Eigen::MatrixXd MotionJacobianOfOrder(const Model& m, const ModelConfig<Order>& mc, const std::string& bodyName)
{
    return MotionJacobianOfOrder(m, mbc, bodyName, JacOrder);
}

template <int Order>
Eigen::MatrixXd MotionJacobianOfOrder(const Model& m, const ModelConfig<Order>& mc, const std::string& bodyName, int jacOrder)
{
    const auto& parents = m.parents();
    const auto& pos = m.jointsPosInDof();
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(6, m.nDof());
    Index j = m.bodyIndexByName(bodyName);
    auto C_b_0 = mc.linkMotions[j].inverse();
    while (j != -1) {
        // A faster code would recursively compute
        auto C_b_j = C_b_0 * mc.linkMotions[j];
        auto S = mb.joint(j).motionSubspace();
        if (jacOrder == 0) {
            J.block(0, pos[j], 6, S.cols()) = C_b_j.transform() * S;
        } else {
            J.block(0, pos[j], 6, S.cols()) = C_b_j[jacOrder - 1] * S;
        }
        j = parents[j];
    }

    return J;
}

} // namespace cdm
