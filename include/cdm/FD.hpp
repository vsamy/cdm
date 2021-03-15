#pragma once

#include "cdm/Model.hpp"
#include "cdm/ModelConfig.hpp"
#include "cdm/typedefs.hpp"
#include "cdm/LowerBlockTriangularMatrix.hpp"

namespace cdm {

template <int Order>
std::vector<Eigen::VectorXd> FD(const Model& m, const ModelConfig<Order>& mc, const std::vector<Eigen::VectorXd>& tau)
{
    const auto& parents = m.jointParents();
    const int dynOrder = static_cast<int>(mc.world.order());

    std::vector<DiMotionSubspace<Order>> G(m.nLinks());
    std::vector<Eigen::Matrix<double, 6 * Order, Eigen::Dynamic>> PA(m.nLinks());
    std::vector<Eigen::MatrixXd> IA(m.nLinks());
    std::vector<LB66> U(m.nLinks());
    std::vector<LB66> UD(m.nLinks());
    std::vector<LB66> D(m.nLinks());
    std::vector<LB66> DInv(m.nLinks());
    std::vector<Eigen::VectorXd> y(dynOrder, Eigen::VectorXd(m.nDof()));
    std::vector<Eigen::Matrix<double, 6 * Order, Eigen::Dynamic>> T(m.nLinks());

    for (Index i = 0; i < m.nLinks(); ++i) {
        G[i] = DiMotionSubspace<Order>{ m.joint(i).S() };
        IA[i].setZero(6 * Order, 6 * Order);
        // UD[i].setOrder(Order).setZero(G[i].cols(), 6);
        // D[i].setOrder(Order).setZero(G[i].cols(), G[i].cols());
        PA[i].setZero();
    }

    for (Index i = m.nLinks() - 1; i >= 0; --i) {
        for (Index j = 0; j < Order; ++j) {
            IA[i].block<6, 6>(6 * j, 6 * j) += m.body(i).inertia().matrix();
        }
        U[i] = IA[i] * G[i];
        auto GT = G[i].transpose();
        UD[i] = GT * IA[i];
        D[i] = GT * U[i];

        DInv[i] = D[i].inverse();
        y[i] = DInv[i] * (tau[i] - (GT * PA[i]).vector());
        if (parents[i] != -1) {
            LB66 tmp1 = IA[i] - U[i] * DInv[i] * UD[i];
            IA[parents[i]] += DualMul(mc.jointMotions[i], tmp1) * mc.jointMotions[i];
            auto tmp2 = PA[i] + U[i] * y[i];
            PA[parents[i]] += DualMul(mc.jointMotions[i], tmp2);
        }
    }

    for (Index i = 0; i < m.nLinks(); ++i) {
        if (parents[i] != -1) {
            y[i] -= DInv[i] * UD[i] * (mc.jointMotions[i] * T[parents[i]]).vector();
        }
        T[i] = G[i] * y[i];
        if (parents[i] != -1) {
            T[i] += mc.jointMotions[i] * T[parents[i]];
        }
    }

    return y;
}

template <int Order>
std::vector<Eigen::VectorXd> FD(const Model& m, const ModelConfig<Order>& mc, const std::vector<Eigen::VectorXd>& tau)
{
    constexpr int ord = Tree::order;
    const auto& parent = m.jointParents();

    std::vector<Eigen::MatrixXd> C(m.nLinks());
    std::vector<Eigen::MatrixXd> CD(m.nLinks());
    std::vector<Eigen::VectorXd> PA(m.nLinks(), Eigen::VectorXd::Zero(6 * ord));
    std::vector<Eigen::MatrixXd> IA(m.nLinks(), Eigen::MatrixXd::Zero(6 * ord, 6 * ord));
    std::vector<Eigen::MatrixXd> G(m.nLinks());
    std::vector<Eigen::MatrixXd> U(m.nLinks());
    std::vector<Eigen::MatrixXd> UD(m.nLinks());
    std::vector<Eigen::MatrixXd> D(m.nLinks());
    std::vector<Eigen::VectorXd> T(m.nLinks());
    std::vector<Eigen::VectorXd> y(m.nLinks());

    for (Index i = 0; i < m.nLinks(); ++i) {
        C[i] = mc.jointMotions[i].inverse().template matrix<ord>();
        CD[i] = mc.jointMotions[i].template dualMatrix<ord>();
        Index dof = m.joint(i).dof();
        const auto& S = m.joint(i).S();
        G[i].setZero(6 * ord, 6 * dof);
        for (Index j = 0; j < ord; ++j) {
            G[i].block(6 * j, j * dof, 6, dof) = S;
        }
    }

    for (int i = mb.nrBodies() - 1; i >= 0; --i) {
        Eigen::Matrix6d I = m.body(i).inertia().matrix();
        Index dof = m.joint(i).dof();
        const auto& S = m.joint(i).S();
        UD[i].setZero(dof, 6 * ord);
        // U[i] = IA[i] * G[i];
        // UD[i] = G[i].transpose() * IA[i];
        // D[i] = G[i].transpose() * U[i];
        for (int j = 0; j < ord; ++j) {
            IA[i].block<6, 6>(6 * j, 6* j) += I;
            U[i].block(0, j * dof, 6 * ord, dof) = IA[i].block<6 * ord, 6>(0, 6 * j) * S;
            UD[i].block(j * dof, 0, dof, 6 * ord) = S.transpose() * IA[i].block<6, 6 * ord>(6 * j, 0);
            D[i].block(j * dof, 0, dof, dof) = S.transpose() * U[i].block(6 * j, 0, 6, 6 * dof);
        }

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

template <int Order>
Eigen::VectorXd standardFD(const Model& m, ModelConfig<Order>& mc, const Eigen::VectorXd& tau)
{
    const auto& parents = m.jointParents();
    const auto& jpd = m.jointPosInDof();

    std::vector<Eigen::Vector6d> PA(m.nLinks());
    std::vector<Eigen::Matrix6d> IA(m.nLinks());
    std::vector<Eigen::Matrix6d> A(m.nLinks());
    std::vector<Eigen::Vector6d> T(m.nLinks());
    std::vector<Eigen::MatrixXd> U(m.nLinks());
    std::vector<Eigen::MatrixXd> D(m.nLinks());
    Eigen::VectorXd y(tau.size());

    for (Index i = 0; i < m.nLinks(); ++i) {
        IA[i].setZero();
        PA[i].setZero();
        A[i] = mc.jointMotions[i].transform().inverse().matrix();
    }

    for (Index i = m.nLinks() - 1; i >= 0; --i) {
        auto S = m.joint(i).S().matrix();
        Index dof = m.joint(i).dof();
        IA[i] += m.body(i).inertia().matrix();
        U[i] = IA[i] * S;
        D[i] = S.transpose() * U[i];

        y.segment(jpd[i], dof) = D[i].inverse() * (tau.segment(jpd[i], dof) - S.transpose() * PA[i]);
        if (parents[i] != -1) {
            auto tmp1 = IA[i] - U[i] * D[i].inverse() * U[i].transpose();
            IA[parents[i]] += A[i].transpose() * tmp1 * A[i];
            auto tmp2 = PA[i] + U[i] * y.segment(jpd[i], dof);
            PA[parents[i]] += A[i].transpose() * tmp2;
        }
    }

    for (Index i = 0; i < m.nLinks(); ++i) {
        Index dof = m.joint(i).dof();
        if (parents[i] != -1) {
            y.segment(jpd[i], dof) -= D[i].inverse() * U[i].transpose() * A[i] * T[parents[i]];
        }
        T[i] = m.joint(i).S().matrix() * y.segment(jpd[i], dof);
        if (parents[i] != -1) {
            T[i] += A[i] * T[parents[i]];
        }
    }

    return y;
}

} // namespace cdm
