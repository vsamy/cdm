#pragma once

#include "LinkMomentumJacobian.hpp"
#include "utils.hpp"

template <typename MI, typename Tree>
Eigen::MatrixXd JointMomentumJacobian(const std::string& bodyName, MI& info, Tree& tree)
{
    constexpr int ord = Tree::order;
    const auto& mb = info.model.mb;
    const auto& pred = mb.predecessors();
    const auto& pos = mb.jointsPosInDof();
    int bInd = mb.bodyIndexByName(bodyName);
    auto C_b_0 = tree.bodies[bInd].inverse();
    // Get p part
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(6 * ord, mb.nrDof() * ord);
    auto M = getSubTreeInertia(mb, tree);
    int j = bInd;
    while (j != -1) {
        auto G = makeDiag<ord>(mb.joint(j).motionSubspace()); // G = diag(S)
        B.block(0, ord * pos[j], 6 * ord, G.cols()) = M[bInd] * (C_b_0 * tree.bodies[j]).template matrix<ord>() * G;
        j = pred[j];
    }

    auto findChildren = [&mb](int b) {
        const auto& pred = mb.predecessors();
        std::vector<int> children;
        for (int i = b + 1; i < mb.nrBodies(); ++i) {
            if (pred[i] == b)
                children.push_back(i);
        }

        return children;
    };

    // Compute \sum Cd * p part
    std::function<void(const rbd::MultiBody& mb, const std::vector<int>&)> computeChildMomentumsAndAdd;
    computeChildMomentumsAndAdd = [&, M](const rbd::MultiBody& mb, const std::vector<int>& children) {
        constexpr int ord = Tree::order;
        const auto& pos = mb.jointsPosInDof();
        for (int c : children) {
            auto C_b_c = C_b_0 * tree.bodies[c];
            auto G = makeDiag<ord>(mb.joint(c).motionSubspace());
            B.block(0, ord * pos[c], 6 * ord, G.cols()) = C_b_c.template dualMatrix<ord>() * M[c] * G;

            computeChildMomentumsAndAdd(mb, findChildren(c));
        }
    };

    computeChildMomentumsAndAdd(mb, findChildren(bInd));
    return B;
}

template <int JacOrder, typename MI, typename Tree>
Eigen::MatrixXd JointMomentumJacobianOfOrder(const std::string& bodyName, MI& info, Tree& tree)
{
    constexpr int jac_ord = JacOrder;
    const auto& mb = info.model.mb;
    std::vector<Eigen::Matrix6d> M(mb.nrBodies(), Eigen::Matrix6d::Zero());
    const auto& bodies = mb.bodies();
    const auto& pred = mb.predecessors();
    const auto& pos = mb.jointsPosInDof();
    for (int i = mb.nrBodies() - 1; i >= 0; --i) {
        M[i] += bodies[i].inertia().matrix();
        int p = pred[i]; // parent
        if (p != -1) {
            auto C_p_b = tree.bodies[p].inverse() * tree.bodies[i];
            if constexpr (jac_ord == 0) {
                M[p] += C_p_b.transform().dualMatrix() * M[i] * C_p_b.transform().inverse().matrix();
            } else {
                M[p] += C_p_b[jac_ord - 1].dualMatrix() * M[i] * C_p_b.inverse()[jac_ord - 1].matrix();
            }
        }
    }

    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(6, mb.nrDof());
    int bInd = mb.bodyIndexByName(bodyName);
    int j = bInd;
    auto C_b_0 = tree.bodies[bInd].inverse();
    while (j != -1) {
        auto S = mb.joint(j).motionSubspace();
        auto C_b_j = C_b_0 * tree.bodies[j];
        if constexpr (jac_ord == 0) {
            B.block(0, pos[j], 6, S.cols()) = M[bInd] * C_b_j.transform().matrix() * S;
        } else {
            B.block(0, pos[j], 6, S.cols()) = M[bInd] * C_b_j[jac_ord - 1].matrix() * S;
        }
        j = pred[j];
    }

    auto findChildren = [&mb](int b) {
        const auto& pred = mb.predecessors();
        std::vector<int> children;
        for (int i = b + 1; i < mb.nrBodies(); ++i) {
            if (pred[i] == b)
                children.push_back(i);
        }

        return children;
    };

    // Compute \sum Cd * p part
    std::function<void(const std::vector<int>&)> computeChildMomentumsAndAdd;
    computeChildMomentumsAndAdd = [&](const std::vector<int>& children) {
        constexpr int ord = Tree::order;
        const auto& pos = mb.jointsPosInDof();
        for (int c : children) {
            auto C_b_c = C_b_0 * tree.bodies[c];
            auto S = mb.joint(c).motionSubspace();
            if constexpr (jac_ord == 0) {
                B.block(0, pos[c], 6, S.cols()) = C_b_c.transform().dualMatrix() * M[c] * S;
            } else {
                B.block(0, pos[c], 6, S.cols()) = C_b_c[jac_ord - 1].dualMatrix() * M[c] * S;
            }

            computeChildMomentumsAndAdd(findChildren(c));
        }
    };

    computeChildMomentumsAndAdd(findChildren(bInd));
    return B;
}