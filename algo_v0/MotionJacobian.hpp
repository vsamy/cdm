#pragma once

#include "utils.hpp"

template <typename MI, typename Tree>
Eigen::MatrixXd MotionJacobian(const std::string& bodyName, MI& info, Tree& tree)
{
    constexpr int ord = Tree::order;
    const auto& mb = info.model.mb;
    const auto& pred = mb.predecessors();
    const auto& pos = mb.jointsPosInDof();
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(6 * ord, mb.nrDof() * ord);
    int j = mb.bodyIndexByName(bodyName);
    auto C_b_0 = tree.links[j].inverse();
    while (j != -1) {
        // A faster code would recursively compute
        auto G = makeDiag<ord>(mb.joint(j).motionSubspace()); // G = diag(S)
        J.block(0, ord * pos[j], 6 * ord, G.cols()) = (C_b_0 * tree.links[j]).template matrix<ord>() * G;
        j = pred[j];
    }

    return J;
}

template <int JacOrder, typename MI, typename Tree>
Eigen::MatrixXd MotionJacobianOfOrder(const std::string& bodyName, MI& info, Tree& tree)
{
    constexpr int jac_ord = JacOrder;
    const auto& mb = info.model.mb;
    const auto& pred = mb.predecessors();
    const auto& pos = mb.jointsPosInDof();
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(6, mb.nrDof());
    int j = mb.bodyIndexByName(bodyName);
    auto C_b_0 = tree.links[j].inverse();
    while (j != -1) {
        // A faster code would recursively compute
        auto C_b_j = C_b_0 * tree.links[j];
        auto S = mb.joint(j).motionSubspace();
        if constexpr (jac_ord == 0) {
            J.block(0, pos[j], 6, S.cols()) = C_b_j.transform().matrix() * S;
        } else {
            J.block(0, pos[j], 6, S.cols()) = C_b_j[jac_ord - 1].matrix() * S;
        }
        j = pred[j];
    }

    return J;
}