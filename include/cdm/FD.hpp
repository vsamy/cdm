#pragma once

namespace cdm {

template <int Order>
std::vector<Eigen::VectorXd> FD(const Model& m, const ModelConfig<Order>& mc, const std::vector<Eigen::VectorXd>& tau)
{
    const auto& parents = m.jointParents();
    const int dynOrder = static_cast<int>(mc.world.order());

    std::vector<CMTM<Order>> C(m.nLinks());
    std::vector<DiMotionSubspace<Order>> G(m.nLinks());
    std::vector<ForceVectorX<Order>> PA(m.nLinks());
    std::vector<LBTM66<Order>> IA(m.nLinks());
    std::vector<LBTM66<Order>> U(m.nLinks());
    std::vector<LBTM66<Order>> UD(m.nLinks());
    std::vector<LBTM66<Order>> D(m.nLinks());
    std::vector<LBTM66<Order>> DInv(m.nLinks());
    std::vector<Eigen::VectorXd> y(dynOrder, Eigen::VectorXd(m.nDof()));
    std::vector<MotionVectorX<Order>> T(m.nLinks());

    for (Index i = 0; i < m.nLinks(); ++i) {
        G[i] = DiMotionSubspace<Order>{ m.joint(i).S() };
        IA[i].setZero(Order);
        PA[i].setZero(Order);
    }

    for (Index i = m.nLinks() - 1; i >= 0; --i) {
        IA[i] += DiInertia<Order>(m.body(i).inertia());
        U[i] = IA[i] * G[i];
        auto GT = G[i].transpose();
        UD[i] = GT * IA[i];
        D[i] = GT * U[i];

        DInv[i] = D[i].inverse();
        y[i] = DInv[i] * (tau[i] - GT * PA[i]);
        if (parents[i] != -1) {
            LBTM66<Order> tmp1 = IA[i] - U[i] * DInv[i] * UD[i];
            IA[parents[i]] += coma::DualMul(mc.jointMotions[i], tmp1) * mc.jointMotions[i];
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
