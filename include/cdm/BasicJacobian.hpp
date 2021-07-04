/*
 * Copyright 2020-2021 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include "cdm/Model.hpp"
#include "cdm/ModelConfig.hpp"
#include "math_utility.hpp"

namespace cdm {

template <int Order>
Eigen::MatrixXd BasicJacobian(const Model& m, const ModelConfig<Order>& mc, const std::string& bodyName)
{
    static_assert(Order >= 0, "Not yet ready for Dynamic");
    // TODO: assert(m.hasBody(bodyName))
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(6 * Order, m.nDof() * Order);
    Index j = m.bodyIndexByName(bodyName);
    auto C_b_0 = mc.bodyMotions[static_cast<size_t>(j)].inverse();
    while (j != -1) {
        size_t uj = static_cast<size_t>(j);
        J.block(0, Order * m.jointPosInDof(j), 6 * Order, Order * m.joint(j).dof()) = (C_b_0 * mc.bodyMotions[uj]).template matrix<Order>() * makeDiag<Order>(m.joint(j).S().matrix());
        j = m.jointParent(j);
    }

    return J;
}

template <int JacOrder, int Order>
Eigen::MatrixXd BasicJacobianOfOrder(const Model& m, const ModelConfig<Order>& mc, const std::string& bodyName)
{
    return BasicJacobianOfOrder(m, mc, bodyName, JacOrder);
}

template <int Order>
Eigen::MatrixXd BasicJacobianOfOrder(const Model& m, const ModelConfig<Order>& mc, const std::string& bodyName, int jacOrder)
{
    static_assert(Order >= 0, "Not yet ready for Dynamic");
    // TODO: assert(m.hasBody(bodyName))
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(6, m.nDof());
    Index j = m.bodyIndexByName(bodyName);
    auto C_b_0 = mc.bodyMotions[static_cast<size_t>(j)].inverse();
    while (j != -1) {
        // A faster code would recursively compute
        size_t uj = static_cast<size_t>(j);
        auto C_b_j = C_b_0 * mc.bodyMotions[uj];
        auto S = m.joint(j).S();
        if (jacOrder == 0) {
            J.block(0, m.jointPosInDof(j), 6, S.cols()) = C_b_j.transform() * S;
        } else {
            J.block(0, m.jointPosInDof(j), 6, S.cols()) = C_b_j[jacOrder - 1] * S;
        }
        j = m.jointParent(j);
    }

    return J;
}

} // namespace cdm
