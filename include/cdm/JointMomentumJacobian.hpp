#pragma once

#include "cdm/LinkMomentumJacobian.hpp"

namespace cdm {

template <int Order>
Eigen::MatrixXd JointMomentumJacobian(const Model& m, const ModelConfig<Order>& mc, const std::string& bodyName)
{
    const auto& posInDof = m.jointPosInDof();
    const auto& parents = m.jointParents();
    Index bInd = m.bodyIndexByName(bodyName);
    auto C_b_0 = mc.bodyMotions[bInd].inverse();
    // Get p part
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(6 * Order, m.nDof() * Order);
    auto M = getSubTreeInertia(m, mc);
    Index j = bInd;
    while (j != -1) {
        auto G = makeDiag<Order>(m.joint(j).S().matrix()); // G = diag(S)
        B.block(0, Order * posInDof[j], 6 * Order, G.cols()) = M[bInd] * (C_b_0 * mc.bodyMotions[j]).template matrix<Order>() * G;
        j = parents[j];
    }

    auto findChildren = [&parents, nLink = m.nLinks()](Index b) {
        std::vector<Index> children;
        for (Index i = b + 1; i < nLink; ++i) {
            if (parents[i] == b)
                children.push_back(i);
        }

        return children;
    };

    // Compute \sum Cd * p part
    std::function<void(const Model& m, const std::vector<Index>&)> computeChildMomentumsAndAdd;
    computeChildMomentumsAndAdd = [&, M](const Model& m, const std::vector<Index>& children) {
        const auto& posInDof = m.jointPosInDof();
        for (auto c : children) {
            auto C_b_c = C_b_0 * mc.bodyMotions[c];
            auto G = makeDiag<Order>(m.joint(c).S().matrix());
            B.block(0, Order * posInDof[c], 6 * Order, G.cols()) = C_b_c.template dualMatrix<Order>() * M[c] * G;

            computeChildMomentumsAndAdd(m, findChildren(c));
        }
    };

    computeChildMomentumsAndAdd(m, findChildren(bInd));
    return B;
}

template <int JacOrder, int Order>
Eigen::MatrixXd JointMomentumJacobianOfOrder(const Model& m, const ModelConfig<Order>& mc, const std::string& bodyName)
{
    std::vector<Eigen::Matrix6d> M(m.nLinks(), Eigen::Matrix6d::Zero());
    const auto& bodies = m.bodies();
    const auto& parents = m.parents();
    const auto& posInDof = m.jointsPosInDof();
    for (Index i = mb.nrBodies() - 1; i >= 0; --i) {
        M[i] += bodies[i].inertia().matrix();
        int p = pred[i]; // parent
        if (p != -1) {
            auto C_p_b = mc.bodyMotions[p].inverse() * mc.bodyMotions[i];
            if constexpr (JacOrder == 0) {
                M[p] += C_p_b.transform().dualMatrix() * M[i] * C_p_b.transform().inverse().matrix();
            } else {
                M[p] += C_p_b[JacOrder - 1].dualMatrix() * M[i] * C_p_b.inverse()[JacOrder - 1].matrix();
            }
        }
    }

    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(6, m.nDof());
    Index bInd = m.bodyIndexByName(bodyName);
    Index j = bInd;
    auto C_b_0 = mc.bodyMotions[bInd].inverse();
    while (j != -1) {
        auto S = mb.joint(j).S();
        auto C_b_j = C_b_0 * mc.bodyMotions[j];
        if constexpr (JacOrder == 0) {
            B.block(0, pos[j], 6, S.cols()) = M[bInd] * C_b_j.transform().matrix() * S.matrix();
        } else {
            B.block(0, pos[j], 6, S.cols()) = M[bInd] * C_b_j[JacOrder - 1].matrix() * S.matrix();
        }
        j = pred[j];
    }

    auto findChildren = [&m](Index b) {
        const auto& parents = m.parents();
        std::vector<Index> children;
        for (int i = b + 1; i < mb.nrBodies(); ++i) {
            if (parents[i] == b)
                children.push_back(i);
        }

        return children;
    };

    // Compute \sum Cd * p part
    std::function<void(const std::vector<Index>&)> computeChildMomentumsAndAdd;
    computeChildMomentumsAndAdd = [&](const std::vector<Index>& children) {
        const auto& posInDof = m.jointsPosInDof();
        for (int c : children) {
            auto C_b_c = C_b_0 * mc.bodyMotions[c];
            auto S = mb.joint(c).S();
            if constexpr (JacOrder == 0) {
                B.block(0, pos[c], 6, S.cols()) = C_b_c.transform().dualMatrix() * M[c] * S.matrix();
            } else {
                B.block(0, pos[c], 6, S.cols()) = C_b_c[JacOrder - 1].dualMatrix() * M[c] * S.matrix();
            }

            computeChildMomentumsAndAdd(findChildren(c));
        }
    };

    computeChildMomentumsAndAdd(findChildren(bInd));
    return B;
}

} // namespace cdm
