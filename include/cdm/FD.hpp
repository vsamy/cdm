#pragma once

#include "cdm/LowerBlockTriangularMatrix.hpp"
#include "cdm/Model.hpp"
#include "cdm/ModelConfig.hpp"
#include "cdm/math_utility.hpp"
#include "cdm/typedefs.hpp"

namespace cdm {

/*! \brief Comprehensive Forward Dynamics.
 * This computes all motion information given comprehensive generalized momentums.
 * \tparam Order Orde of the model.
 * \param m Model to compute the FD from.
 * \param mc Model information where to save the results.
 * \param tau Set of comprehensive generalized momentums.
 */
template <int Order>
std::vector<Eigen::VectorXd> FD(const Model& m, const ModelConfig<Order>& mc, const std::vector<Eigen::VectorXd>& tau)
{
    constexpr int ord = Order;

    size_t un = static_cast<size_t>(m.nLinks());
    std::vector<Eigen::MatrixXd> C(un);
    std::vector<Eigen::MatrixXd> CD(un);
    std::vector<Eigen::VectorXd> PA(un);
    std::vector<Eigen::MatrixXd> IA(un);
    std::vector<Eigen::MatrixXd> G(un);
    std::vector<Eigen::MatrixXd> U(un);
    std::vector<Eigen::MatrixXd> UD(un);
    std::vector<Eigen::MatrixXd> D(un);
    std::vector<Eigen::VectorXd> T(un);
    std::vector<Eigen::VectorXd> y(un);

    for (size_t i = 0; i < un; ++i) {
        C[i] = mc.jointMotions[i].inverse().template matrix<ord>();
        CD[i] = mc.jointMotions[i].template dualMatrix<ord>();
        G[i] = makeDiag<Order>(m.joints()[i].S().matrix());
        IA[i].setZero(6 * ord, 6 * ord);
        PA[i].setZero(6 * ord);
    }

    for (Index i = m.nLinks() - 1; i >= 0; --i) {
        size_t ui = static_cast<size_t>(i);
        IA[ui] += makeDiag<Order>(m.body(i).inertia().matrix());
        U[ui] = IA[ui] * G[ui];
        UD[ui] = G[ui].transpose() * IA[ui];
        D[ui] = G[ui].transpose() * U[ui];

        y[ui] = D[ui].inverse() * (tau[ui] - G[ui].transpose() * PA[ui]);
        Index p = m.jointParent(i);
        if (p != -1) {
            size_t up = static_cast<size_t>(p);
            auto tmp1 = IA[ui] - U[ui] * D[ui].inverse() * UD[ui];
            IA[up] += CD[ui] * tmp1 * C[ui];
            auto tmp2 = PA[ui] + U[ui] * y[ui];
            PA[up] += CD[ui] * tmp2;
        }
    }

    T[0] = G[0] * y[0];
    for (size_t i = 1; i < un; ++i) {
        size_t up = static_cast<size_t>(m.jointParents()[i]);
        y[i] -= D[i].inverse() * UD[i] * C[i] * T[up];
        T[i] = G[i] * y[i] + C[i] * T[up];
    }

    return y;
}

// template <int Order>
// std::vector<Eigen::VectorXd> FD(const Model& m, const ModelConfig<Order>& mc, const std::vector<Eigen::VectorXd>& tau)
// {
//     constexpr int ord = Order;
//     const auto& parents = m.jointParents();

//     std::vector<Eigen::MatrixXd> C(m.nLinks());
//     std::vector<Eigen::MatrixXd> CD(m.nLinks());
//     std::vector<Eigen::VectorXd> PA(m.nLinks(), Eigen::VectorXd::Zero(6 * ord));
//     std::vector<Eigen::MatrixXd> IA(m.nLinks(), Eigen::MatrixXd::Zero(6 * ord, 6 * ord));
//     std::vector<Eigen::MatrixXd> G(m.nLinks());
//     std::vector<Eigen::MatrixXd> U(m.nLinks());
//     std::vector<Eigen::MatrixXd> UD(m.nLinks());
//     std::vector<Eigen::MatrixXd> D(m.nLinks());
//     std::vector<Eigen::VectorXd> T(m.nLinks());
//     std::vector<Eigen::VectorXd> y(m.nLinks());

//     for (Index i = 0; i < m.nLinks(); ++i) {
//         C[i] = mc.jointMotions[i].inverse().template matrix<ord>();
//         CD[i] = mc.jointMotions[i].template dualMatrix<ord>();
//         Index dof = m.joint(i).dof();
//         const auto& S = m.joint(i).S();
//         G[i].setZero(6 * ord, 6 * dof);
//         for (Index j = 0; j < ord; ++j) {
//             G[i].block(6 * j, j * dof, 6, dof) = S;
//         }
//     }

//     for (int i = mb.nrBodies() - 1; i >= 0; --i) {
//         Eigen::Matrix6d I = m.body(i).inertia().matrix();
//         Index dof = m.joint(i).dof();
//         const auto& S = m.joint(i).S();
//         UD[i].setZero(dof, 6 * ord);
//         // U[i] = IA[i] * G[i];
//         // UD[i] = G[i].transpose() * IA[i];
//         // D[i] = G[i].transpose() * U[i];
//         for (int j = 0; j < ord; ++j) {
//             IA[i].block<6, 6>(6 * j, 6 * j) += I;
//             U[i].block(0, j * dof, 6 * ord, dof) = IA[i].block<6 * ord, 6>(0, 6 * j) * S;
//             UD[i].block(j * dof, 0, dof, 6 * ord) = S.transpose() * IA[i].block<6, 6 * ord>(6 * j, 0);
//             D[i].block(j * dof, 0, dof, dof) = S.transpose() * U[i].block(6 * j, 0, 6, 6 * dof);
//             T[i].segment<6>(6 * j) = Tau[i].segment<6>(6 * j) - S.transpose() * PA[i].segment<6>(6 * j);
//         }

//         // y[i] = D[i].inverse() * (Tau[i] - G[i].transpose() * PA[i]);
//         y[i] = D[i].inverse() * T[i];
//         if (pred[i] != -1) {
//             auto tmp1 = IA[i] - U[i] * D[i].inverse() * UD[i];
//             IA[pred[i]] += CD[i] * tmp1 * C[i];
//             auto tmp2 = PA[i] + U[i] * y[i];
//             PA[pred[i]] += CD[i] * tmp2;
//         }
//     }

//     int pos = 0;
//     for (int j = 0; j < ord; ++j) {
//         T[0].segment<6>(6 * j) = S * y[0].segment(pos, dof);
//         pos += dof;
//     }
//     for (int i = 0; i < mb.nrJoints(); ++i) {
//         int dof = mb.joint(i).dof();
//         const auto& S = mb.joint(i).S();
//         y[i] -= D[i].inverse() * UD[i] * C[i] * T[pred[i]];
//         pos = 0;
//         for (int j = 0; j < ord; ++j) {
//             T[i].segment<6>(6 * j) = S * y[i].segment(pos, dof);
//             pos += dof;
//         }
//         T[i] += C[i] * T[pred[i]];
//     }

//     return y;
// }

/*! \brief Standard Forward Dynamics.
 * This can be used to perform recursive FD when comprehensive torques are provided.
 * \tparam Order Orde of the model.
 * \param m Model to compute the FD from.
 * \param mc Model information where to save the results.
 * \param tau Set of comprehensive torques.
 */
template <int Order>
Eigen::VectorXd standardFD(const Model& m, ModelConfig<Order>& mc, const Eigen::VectorXd& tau)
{
    size_t un = static_cast<size_t>(m.nLinks());
    std::vector<Eigen::Vector6d> PA(un);
    std::vector<Eigen::Matrix6d> IA(un);
    std::vector<Eigen::Matrix6d> A(un);
    std::vector<Eigen::Vector6d> T(un);
    std::vector<Eigen::MatrixXd> U(un);
    std::vector<Eigen::MatrixXd> D(un);
    Eigen::VectorXd y(tau.size());

    for (size_t i = 0; i < un; ++i) {
        IA[i].setZero();
        PA[i].setZero();
        A[i] = mc.jointMotions[i].transform().inverse().matrix();
    }

    for (Index i = m.nLinks() - 1; i >= 0; --i) {
        size_t ui = static_cast<size_t>(i);
        const auto& S = m.joint(i).S().matrix();
        Index dof = m.joint(i).dof();
        IA[ui] += m.body(i).inertia().matrix();
        U[ui] = IA[ui] * S;
        D[ui] = S.transpose() * U[ui];

        Index posInDof = m.jointPosInDof(i);
        y.segment(posInDof, dof) = D[ui].inverse() * (tau.segment(posInDof, dof) - S.transpose() * PA[ui]);
        Index p = m.jointParent(i);
        if (p != -1) {
            size_t up = static_cast<size_t>(p);
            auto tmp1 = IA[ui] - U[ui] * D[ui].inverse() * U[ui].transpose();
            IA[up] += A[ui].transpose() * tmp1 * A[ui];
            auto tmp2 = PA[ui] + U[ui] * y.segment(posInDof, dof);
            PA[up] += A[ui].transpose() * tmp2;
        }
    }

    T[0] = m.joint(0).S().matrix() * y.segment(0, m.joint(0).dof());
    for (Index i = 1; i < m.nLinks(); ++i) {
        size_t ui = static_cast<size_t>(i);
        size_t up = static_cast<size_t>(m.jointParent(i));
        Index dof = m.joint(i).dof();
        Index posInDof = m.jointPosInDof(i);
        y.segment(posInDof, dof) -= D[ui].inverse() * U[ui].transpose() * A[ui] * T[up];
        T[ui] = m.joint(i).S().matrix() * y.segment(posInDof, dof) + A[ui] * T[up];
    }

    return y;
}

} // namespace cdm
