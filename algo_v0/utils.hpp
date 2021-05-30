#pragma once

#include "coma/Core"
#include <Eigen/Core>
#include <RBDyn/MultiBody.h>

template <int Order>
Eigen::MatrixXd makeDiag(const Eigen::MatrixXd& mat)
{
    constexpr int ord = Order;
    Eigen::MatrixXd out = Eigen::MatrixXd::Zero(ord * mat.rows(), ord * mat.cols());
    for (int i = 0; i < ord; ++i)
        out.block(i * mat.rows(), i * mat.cols(), mat.rows(), mat.cols()) = mat;

    return out;
}

template <typename CrossNOp>
Eigen::MatrixXd generateD(const CrossNOp& cx)
{
    constexpr int n_vec = coma::internal::traits<CrossNOp>::n_vec;
    using Scalar = typename coma::internal::traits<CrossNOp>::Scalar;
    using mat_t = Eigen::Matrix<Scalar, 6 * n_vec, 6 * n_vec>;
    using sub_mat_t = Eigen::Matrix<Scalar, 6, 6 * n_vec>;
    Eigen::MatrixXd D_N = Eigen::MatrixXd::Zero(6 * (n_vec - 1), 6 * (n_vec - 1));
    for (int i = 0; i < n_vec - 1; ++i)
        D_N.block<6, 6>(6 * i, 6 * i) = Eigen::Matrix6d::Identity() / (i + 1);
    return mat_t::Identity() + (mat_t() << sub_mat_t::Zero(), D_N * cx.dualMatrix().template topRows<6 * (n_vec - 1)>()).finished();
}

// M = I + Cd * I * C
template <typename Tree>
std::vector<Eigen::MatrixXd> getSubTreeInertia(const rbd::MultiBody& mb, const Tree& tree)
{
    constexpr int ord = Tree::order;
    std::vector<Eigen::MatrixXd> M(mb.nrBodies(), Eigen::MatrixXd::Zero(6 * ord, 6 * ord));
    const auto& bodies = mb.bodies();
    const auto& pred = mb.predecessors();
    for (int i = mb.nrBodies() - 1; i >= 0; --i) {
        M[i] += makeDiag<ord>(bodies[i].inertia().matrix());
        int p = pred[i]; // parent
        if (p != -1) {
            auto C_p_b = tree.links[p].inverse() * tree.links[i];
            M[p] += C_p_b.template dualMatrix<ord>() * M[i] * C_p_b.inverse().template matrix<ord>();
        }
    }

    return M;
}

namespace detail {

template <typename Transf, typename CMTM, size_t Order>
struct CMTMSetter {
    static void impl(const Transf& T, const Eigen::MatrixXd& S, const std::vector<Eigen::MatrixXd>& dqs, int t, int pos, int dof, CMTM& joints);
};

template <typename Transf, typename CMTM>
struct CMTMSetter<Transf, CMTM, 0> {
    static void impl(const Transf& T, const Eigen::MatrixXd& S, const std::vector<Eigen::MatrixXd>& dqs, int t, int pos, int dof, CMTM& joints)
    {
        joints.set(T);
    }
};

template <typename Transf, typename CMTM>
struct CMTMSetter<Transf, CMTM, 1> {
    static void impl(const Transf& T, const Eigen::MatrixXd& S, const std::vector<Eigen::MatrixXd>& dqs, int t, int pos, int dof, CMTM& joints)
    {
        joints.set(T, S * dqs[0].col(t).segment(pos, dof));
    }
};

template <typename Transf, typename CMTM>
struct CMTMSetter<Transf, CMTM, 2> {
    static void impl(const Transf& T, const Eigen::MatrixXd& S, const std::vector<Eigen::MatrixXd>& dqs, int t, int pos, int dof, CMTM& joints)
    {
        joints.set(T, S * dqs[0].col(t).segment(pos, dof), S * dqs[1].col(t).segment(pos, dof));
    }
};

template <typename Transf, typename CMTM>
struct CMTMSetter<Transf, CMTM, 3> {
    static void impl(const Transf& T, const Eigen::MatrixXd& S, const std::vector<Eigen::MatrixXd>& dqs, int t, int pos, int dof, CMTM& joints)
    {
        joints.set(T, S * dqs[0].col(t).segment(pos, dof), S * dqs[1].col(t).segment(pos, dof), S * dqs[2].col(t).segment(pos, dof));
    }
};

template <typename Transf, typename CMTM>
struct CMTMSetter<Transf, CMTM, 4> {
    static void impl(const Transf& T, const Eigen::MatrixXd& S, const std::vector<Eigen::MatrixXd>& dqs, int t, int pos, int dof, CMTM& joints)
    {
        joints.set(T, S * dqs[0].col(t).segment(pos, dof), S * dqs[1].col(t).segment(pos, dof), S * dqs[2].col(t).segment(pos, dof),
            S * dqs[3].col(t).segment(pos, dof));
    }
};

template <typename Transf, typename CMTM>
struct CMTMSetter<Transf, CMTM, 5> {
    static void impl(const Transf& T, const Eigen::MatrixXd& S, const std::vector<Eigen::MatrixXd>& dqs, int t, int pos, int dof, CMTM& joints)
    {
        joints.set(T, S * dqs[0].col(t).segment(pos, dof), S * dqs[1].col(t).segment(pos, dof), S * dqs[2].col(t).segment(pos, dof),
            S * dqs[3].col(t).segment(pos, dof), S * dqs[4].col(t).segment(pos, dof));
    }
};

template <typename Transf, typename CMTM>
struct CMTMSetter<Transf, CMTM, 6> {
    static void impl(const Transf& T, const Eigen::MatrixXd& S, const std::vector<Eigen::MatrixXd>& dqs, int t, int pos, int dof, CMTM& joints)
    {
        joints.set(T, S * dqs[0].col(t).segment(pos, dof), S * dqs[1].col(t).segment(pos, dof), S * dqs[2].col(t).segment(pos, dof),
            S * dqs[3].col(t).segment(pos, dof), S * dqs[4].col(t).segment(pos, dof), S * dqs[5].col(t).segment(pos, dof));
    }
};

template <typename Transf, typename CMTM>
struct CMTMSetter<Transf, CMTM, 7> {
    static void impl(const Transf& T, const Eigen::MatrixXd& S, const std::vector<Eigen::MatrixXd>& dqs, int t, int pos, int dof, CMTM& joints)
    {
        joints.set(T, S * dqs[0].col(t).segment(pos, dof), S * dqs[1].col(t).segment(pos, dof), S * dqs[2].col(t).segment(pos, dof),
            S * dqs[3].col(t).segment(pos, dof), S * dqs[4].col(t).segment(pos, dof), S * dqs[5].col(t).segment(pos, dof),
            S * dqs[6].col(t).segment(pos, dof));
    }
};

template <typename Transf, typename CMTM>
struct CMTMSetter<Transf, CMTM, 8> {
    static void impl(const Transf& T, const Eigen::MatrixXd& S, const std::vector<Eigen::MatrixXd>& dqs, int t, int pos, int dof, CMTM& joints)
    {
        joints.set(T, S * dqs[0].col(t).segment(pos, dof), S * dqs[1].col(t).segment(pos, dof), S * dqs[2].col(t).segment(pos, dof),
            S * dqs[3].col(t).segment(pos, dof), S * dqs[4].col(t).segment(pos, dof), S * dqs[5].col(t).segment(pos, dof),
            S * dqs[6].col(t).segment(pos, dof), S * dqs[7].col(t).segment(pos, dof));
    }
};

template <typename Transf, typename CMTM>
struct CMTMSetter<Transf, CMTM, 9> {
    static void impl(const Transf& T, const Eigen::MatrixXd& S, const std::vector<Eigen::MatrixXd>& dqs, int t, int pos, int dof, CMTM& joints)
    {
        joints.set(T, S * dqs[0].col(t).segment(pos, dof), S * dqs[1].col(t).segment(pos, dof), S * dqs[2].col(t).segment(pos, dof),
            S * dqs[3].col(t).segment(pos, dof), S * dqs[4].col(t).segment(pos, dof), S * dqs[5].col(t).segment(pos, dof),
            S * dqs[6].col(t).segment(pos, dof), S * dqs[7].col(t).segment(pos, dof), S * dqs[8].col(t).segment(pos, dof));
    }
};

template <typename Transf, typename CMTM>
struct CMTMSetter<Transf, CMTM, 10> {
    static void impl(const Transf& T, const Eigen::MatrixXd& S, const std::vector<Eigen::MatrixXd>& dqs, int t, int pos, int dof, CMTM& joints)
    {
        joints.set(T, S * dqs[0].col(t).segment(pos, dof), S * dqs[1].col(t).segment(pos, dof), S * dqs[2].col(t).segment(pos, dof),
            S * dqs[3].col(t).segment(pos, dof), S * dqs[4].col(t).segment(pos, dof), S * dqs[5].col(t).segment(pos, dof),
            S * dqs[6].col(t).segment(pos, dof), S * dqs[7].col(t).segment(pos, dof), S * dqs[8].col(t).segment(pos, dof),
            S * dqs[9].col(t).segment(pos, dof));
    }
};

template <typename Transf, typename CMTM>
struct CMTMSetter<Transf, CMTM, 11> {
    static void impl(const Transf& T, const Eigen::MatrixXd& S, const std::vector<Eigen::MatrixXd>& dqs, int t, int pos, int dof, CMTM& joints)
    {
        joints.set(T, S * dqs[0].col(t).segment(pos, dof), S * dqs[1].col(t).segment(pos, dof), S * dqs[2].col(t).segment(pos, dof),
            S * dqs[3].col(t).segment(pos, dof), S * dqs[4].col(t).segment(pos, dof), S * dqs[5].col(t).segment(pos, dof),
            S * dqs[6].col(t).segment(pos, dof), S * dqs[7].col(t).segment(pos, dof), S * dqs[8].col(t).segment(pos, dof),
            S * dqs[9].col(t).segment(pos, dof), S * dqs[10].col(t).segment(pos, dof));
    }
};

template <typename Transf, typename CMTM>
struct CMTMSetter<Transf, CMTM, 12> {
    static void impl(const Transf& T, const Eigen::MatrixXd& S, const std::vector<Eigen::MatrixXd>& dqs, int t, int pos, int dof, CMTM& joints)
    {
        joints.set(T, S * dqs[0].col(t).segment(pos, dof), S * dqs[1].col(t).segment(pos, dof), S * dqs[2].col(t).segment(pos, dof),
            S * dqs[3].col(t).segment(pos, dof), S * dqs[4].col(t).segment(pos, dof), S * dqs[5].col(t).segment(pos, dof),
            S * dqs[6].col(t).segment(pos, dof), S * dqs[7].col(t).segment(pos, dof), S * dqs[8].col(t).segment(pos, dof),
            S * dqs[9].col(t).segment(pos, dof), S * dqs[10].col(t).segment(pos, dof), S * dqs[11].col(t).segment(pos, dof));
    }
};

}

template<typename Transf, typename CMTM, size_t Order>
void CMTMSet(const Transf& T, const Eigen::MatrixXd& S, const std::vector<Eigen::MatrixXd>& dqs, int t, int pos, int dof, CMTM& joints)
{
    detail::CMTMSetter<Transf, CMTM, Order>::impl(T, S, dqs, t, pos, dof, joints);
}
