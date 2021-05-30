#pragma once

#include "utils.hpp"
#include <Eigen/Core>

template <typename MI, typename Tree>
std::vector<Eigen::VectorXd> FD(const MI& info, Tree& tree, const std::vector<Eigen::VectorXd>& Tau)
{
    constexpr int ord = Tree::order;
    const auto& mb = info.model.mb;
    const auto& pred = mb.predecessors();
    const auto& succ = mb.successors();

    std::vector<Eigen::MatrixXd> C(mb.nrBodies());
    std::vector<Eigen::MatrixXd> CD(mb.nrBodies());
    std::vector<Eigen::VectorXd> PA(mb.nrBodies());
    std::vector<Eigen::MatrixXd> IA(mb.nrBodies());
    std::vector<Eigen::MatrixXd> G(mb.nrBodies());
    std::vector<Eigen::MatrixXd> U(mb.nrBodies());
    std::vector<Eigen::MatrixXd> UD(mb.nrBodies());
    std::vector<Eigen::MatrixXd> D(mb.nrBodies());
    std::vector<Eigen::VectorXd> T(mb.nrBodies());
    std::vector<Eigen::VectorXd> y(mb.nrBodies());

    for (int i = 0; i < mb.nrBodies(); ++i) {
        C[i] = tree.joints[i].inverse().template matrix<ord>();
        CD[i] = tree.joints[i].template dualMatrix<ord>();
        G[i] = makeDiag<ord>(mb.joint(i).motionSubspace());
        IA[i].setZero(6 * ord, 6 * ord);
        PA[i].setZero(6 * ord);
    }

    for (int i = mb.nrBodies() - 1; i >= 0; --i) {
        IA[i] += makeDiag<ord>(mb.body(i).inertia().matrix());
        U[i] = IA[i] * G[i];
        UD[i] = G[i].transpose() * IA[i];
        D[i] = G[i].transpose() * U[i];

        y[i] = D[i].inverse() * (Tau[i] - G[i].transpose() * PA[i]);
        if (pred[i] != -1) {
            auto tmp1 = IA[i] - U[i] * D[i].inverse() * UD[i];
            IA[pred[i]] += CD[i] * tmp1 * C[i];
            auto tmp2 = PA[i] + U[i] * y[i];
            PA[pred[i]] += CD[i] * tmp2;
        }
    }

    for (int i = 0; i < mb.nrJoints(); ++i) {
        int dof = mb.joint(i).dof();
        if (pred[i] != -1) {
            y[i] -= D[i].inverse() * UD[i] * C[i] * T[pred[i]];
        }
        T[i] = G[i] * y[i];
        if (pred[i] != -1) {
            T[i] += C[i] * T[pred[i]];
        }
    }

    return y;
}

template <typename MI, typename Tree>
Eigen::VectorXd standard_FD(const MI& info, Tree& tree, const Eigen::VectorXd& tau)
{
    constexpr int ord = Tree::order;
    const auto& mb = info.model.mb;
    const auto& pred = mb.predecessors();
    const auto& succ = mb.successors();
    const auto& jpd = mb.jointsPosInDof();

    std::vector<Eigen::Vector6d> PA(mb.nrBodies());
    std::vector<Eigen::Matrix6d> IA(mb.nrBodies());
    std::vector<Eigen::Matrix6d> X(mb.nrBodies());
    std::vector<Eigen::Vector6d> T(mb.nrBodies());
    std::vector<Eigen::MatrixXd> U(mb.nrBodies());
    std::vector<Eigen::MatrixXd> D(mb.nrBodies());
    Eigen::VectorXd y(tau.size());

    for (int i = 0; i < mb.nrBodies(); ++i) {
        IA[i].setZero();
        PA[i].setZero();
        X[i] = tree.joints[i].transform().inverse().matrix();
    }

    for (int i = mb.nrBodies() - 1; i >= 0; --i) {
        const auto& S = mb.joint(i).motionSubspace();
        int dof = mb.joint(i).dof();
        IA[i] += mb.body(i).inertia().matrix();
        U[i] = IA[i] * S;
        D[i] = S.transpose() * U[i];

        y.segment(jpd[i], dof) = D[i].inverse() * (tau.segment(jpd[i], dof) - S.transpose() * PA[i]);
        if (pred[i] != -1) {
            auto tmp1 = IA[i] - U[i] * D[i].inverse() * U[i].transpose();
            IA[pred[i]] += X[i].transpose() * tmp1 * X[i];
            auto tmp2 = PA[i] + U[i] * y.segment(jpd[i], dof);
            PA[pred[i]] += X[i].transpose() * tmp2;
        }
    }

    for (int i = 0; i < mb.nrJoints(); ++i) {
        int dof = mb.joint(i).dof();
        if (pred[i] != -1) {
            y.segment(jpd[i], dof) -= D[i].inverse() * U[i].transpose() * X[i] * T[pred[i]];
        }
        T[i] = mb.joint(i).motionSubspace() * y.segment(jpd[i], dof);
        if (pred[i] != -1) {
            T[i] += X[i] * T[pred[i]];
        }
    }

    return y;
}
