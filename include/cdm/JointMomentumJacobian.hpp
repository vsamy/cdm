/*
 * Copyright 2020-2021 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include "cdm/LinkMomentumJacobian.hpp"

namespace cdm {

template <int Order>
Eigen::MatrixXd JointMomentumJacobian(const Model& m, const ModelConfig<Order>& mc, const std::string& bodyName)
{
    static_assert(Order >= 0, "Not yet ready for Dynamic");
    // TODO: assert(m.hasBody(bodyName))
    Index bInd = m.bodyIndexByName(bodyName);
    size_t ubInd = static_cast<size_t>(bInd);
    auto C_b_0 = mc.bodyMotions[ubInd].inverse();
    // Get p part
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(6 * Order, m.nDof() * Order);
    auto M = getSubTreeInertia(m, mc);
    Index j = bInd;
    while (j != -1) {
        auto G = makeDiag<Order>(m.joint(j).S().matrix()); // G = diag(S)
        B.block(0, Order * m.jointPosInDof(j), 6 * Order, G.cols()) = M[ubInd] * (C_b_0 * mc.bodyMotions[static_cast<size_t>(j)]).template matrix<Order>() * G;
        j = m.jointParent(j);
    }

    auto findChildren = [&m](Index b) {
        std::vector<Index> children;
        for (Index i = b + 1; i < m.nLinks(); ++i) {
            if (m.jointParent(i) == b)
                children.push_back(i);
        }

        return children;
    };

    // Compute \sum Cd * p part
    std::function<void(const Model& m, const std::vector<Index>&)> computeChildMomentumsAndAdd;
    computeChildMomentumsAndAdd = [&, M](const Model& m, const std::vector<Index>& children) {
        for (auto c : children) {
            size_t uc = static_cast<size_t>(c);
            auto C_b_c = C_b_0 * mc.bodyMotions[uc];
            auto G = makeDiag<Order>(m.joint(c).S().matrix());
            B.block(0, Order * m.jointPosInDof(c), 6 * Order, G.cols()) = C_b_c.template dualMatrix<Order>() * M[uc] * G;

            computeChildMomentumsAndAdd(m, findChildren(c));
        }
    };

    computeChildMomentumsAndAdd(m, findChildren(bInd));
    return B;
}

template <int JacOrder, int Order>
Eigen::MatrixXd JointMomentumJacobianOfOrder(const Model& m, const ModelConfig<Order>& mc, const std::string& bodyName)
{
    static_assert(Order >= 0, "Not yet ready for Dynamic");

    std::vector<Eigen::Matrix6d> M(static_cast<size_t>(m.nLinks()), Eigen::Matrix6d::Zero());
    const auto& bodies = m.bodies();
    const auto& parents = m.jointParents();
    for (Index i = m.nLinks() - 1; i >= 0; --i) {
        size_t ui = static_cast<size_t>(i);
        M[ui] += bodies[ui].inertia().matrix();
        Index p = parents[ui]; // parent
        if (p != -1) {
            size_t up = static_cast<size_t>(p);
            auto C_p_b = mc.bodyMotions[up].inverse() * mc.bodyMotions[ui];
            if constexpr (JacOrder == 0) {
                M[up] += C_p_b.transform().dualMatrix() * M[ui] * C_p_b.transform().inverse().matrix();
            } else {
                constexpr size_t JOIndex = static_cast<size_t>(JacOrder - 1);
                M[up] += C_p_b[JOIndex].dualMatrix() * M[ui] * C_p_b.inverse()[JOIndex].matrix();
            }
        }
    }

    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(6, m.nDof());
    Index bInd = m.bodyIndexByName(bodyName);
    Index j = bInd;
    auto C_b_0 = mc.bodyMotions[bInd].inverse();
    while (j != -1) {
        auto S = m.joint(j).S();
        auto C_b_j = C_b_0 * mc.bodyMotions[j];
        if constexpr (JacOrder == 0) {
            B.block(0, m.jointPosInDof(j), 6, S.cols()) = M[bInd] * C_b_j.transform().matrix() * S.matrix();
        } else {
            constexpr size_t JOIndex = static_cast<size_t>(JacOrder - 1);
            B.block(0, m.jointPosInDof(j), 6, S.cols()) = M[bInd] * C_b_j[JOIndex].matrix() * S.matrix();
        }
        j = m.jointParent(j);
    }

    auto findChildren = [&m](Index b) {
        std::vector<Index> children;
        for (Index i = b + 1; i < m.nLinks(); ++i) {
            if (m.jointParent(i) == b) {
                children.push_back(i);
            }
        }

        return children;
    };

    // Compute \sum Cd * p part
    std::function<void(const std::vector<Index>&)> computeChildMomentumsAndAdd;
    computeChildMomentumsAndAdd = [&](const std::vector<Index>& children) {
        for (Index c : children) {
            size_t uc = static_cast<size_t>(c);
            auto C_b_c = C_b_0 * mc.bodyMotions[uc];
            auto S = m.joint(c).S();
            if constexpr (JacOrder == 0) {
                B.block(0, m.jointPosInDof(c), 6, S.cols()) = C_b_c.transform().dualMatrix() * M[uc] * S.matrix();
            } else {
                constexpr size_t JOIndex = static_cast<size_t>(JacOrder - 1);
                B.block(0, m.jointPosInDof(c), 6, S.cols()) = C_b_c[JOIndex].dualMatrix() * M[uc] * S.matrix();
            }

            computeChildMomentumsAndAdd(findChildren(c));
        }
    };

    computeChildMomentumsAndAdd(findChildren(bInd));
    return B;
}

} // namespace cdm
