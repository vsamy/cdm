#pragma once

namespace cdm {

template <int Order>
Eigen::VectorXd FD(const Model& m, const ModelConfig<Order>& mc, const std::vector<Eigen::VectorXd>& Tau)
{
    const auto& parents = m.jointParents();

    std::vector<CMTM<Order>> C(m.nLinks());
    std::vector<DiMotionSubspace<Order>> G(m.nLinks());
    std::vector<ForceVectorX<Order>> PA(m.nLinks());
    std::vector<BTM66<Order>> IA(m.nLinks());
    std::vector<DBX<Order>> U(m.nLinks());
    std::vector<DBX<Order>> UD(m.nLinks());
    std::vector<DBX<Order>> D(m.nLinks());
    std::vector<DBX<Order>> DInv(m.nLinks());
    Eigen::VectorXd y(m.nDof());
    std::vector<MotionVectorX<Order>> T(m.nLinks());

    for (Index i = 0; i < m.nLinks(); ++i) {
        G[i] = DiMotionSubspace{ mb.joint(i).S() };
        IA[i].setZero(Order);
        PA[i].setZero(Order);
    }

    for (Index i = m.nLinks() - 1; i >= 0; --i) {
        IA[i] += DiInertia(mb.body(i).inertia());
        U[i] = IA[i] * G[i];
        auto GT = G[i].transpose();
        UD[i] = GT * IA[i];
        D[i] = GT * U[i];

        DInv[i] = D[i].inverse();
        y[i] = DInv[i] * (Tau[i] - GT * PA[i]);
        if (parents[i] != -1) {
            auto tmp1 = IA[i] - U[i] * DInv[i] * UD[i];
            IA[parents[i]] += DualMul(mc.jointMotions[i], tmp1) * mc.jointMotions[i];
            auto tmp2 = PA[i] + U[i] * y[i];
            PA[parents[i]] += DualMul(mc.jointMotions[i], tmp2);
        }
    }

    for (Index i = 0; i < m.nLinks(); ++i) {
        Index dof = m.joint(i).dof();
        if (parents[i] != -1) {
            y[i] -= DInv[i] * UD[i] * (mc.jointMotions[i] * T[parents[i]]);
        }
        T[i] = G[i] * y[i];
        if (parents[i] != -1) {
            T[i] += mc.jointMotions[i] * T[parents[i]];
        }
    }

    return y;
}

// template <int Order>
// Eigen::VectorXd standardFD(const Model& m, ModelConfig<Order>& mc, const std::vector<Eigen::VectorXd>& Tau)
// {
//     constexpr int ord = Tree::order;
//     const auto& mb = info.model.mb;
//     const auto& pred = mb.predecessors();
//     const auto& succ = mb.successors();
//     const auto& jpd = mb.jointsPosInDof();

//     std::vector<Eigen::Vector6d> PA(mb.nrBodies());
//     std::vector<Eigen::Matrix6d> IA(mb.nrBodies());
//     std::vector<Eigen::Matrix6d> X(mb.nrBodies());
//     std::vector<Eigen::Vector6d> T(mb.nrBodies());
//     std::vector<Eigen::MatrixXd> U(mb.nrBodies());
//     std::vector<Eigen::MatrixXd> D(mb.nrBodies());
//     Eigen::VectorXd y(tau.size());

//     for (Index i = 0; i < mb.nrBodies(); ++i) {
//         IA[i].setZero();
//         PA[i].setZero();
//         X[i] = tree.joints[i].transform().inverse().matrix();
//     }

//     for (Index i = mb.nrBodies() - 1; i >= 0; --i) {
//         const auto& S = mb.joint(i).motionSubspace();
//         Index dof = mb.joint(i).dof();
//         IA[i] += mb.body(i).inertia().matrix();
//         U[i] = IA[i] * S;
//         D[i] = S.transpose() * U[i];

//         y.segment(jpd[i], dof) = D[i].inverse() * (tau.segment(jpd[i], dof) - S.transpose() * PA[i]);
//         if (pred[i] != -1) {
//             auto tmp1 = IA[i] - U[i] * D[i].inverse() * U[i].transpose();
//             IA[pred[i]] += X[i].transpose() * tmp1 * X[i];
//             auto tmp2 = PA[i] + U[i] * y.segment(jpd[i], dof);
//             PA[pred[i]] += X[i].transpose() * tmp2;
//         }
//     }

//     for (Index i = 0; i < mb.nrJoints(); ++i) {
//         Index dof = mb.joint(i).dof();
//         if (pred[i] != -1) {
//             y.segment(jpd[i], dof) -= D[i].inverse() * U[i].transpose() * X[i] * T[pred[i]];
//         }
//         T[i] = mb.joint(i).motionSubspace() * y.segment(jpd[i], dof);
//         if (pred[i] != -1) {
//             T[i] += X[i] * T[pred[i]];
//         }
//     }

//     return y;
// }

} // namespace cdm
