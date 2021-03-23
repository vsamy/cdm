#pragma once

#include <Eigen/Core>

namespace cdm {

/*! \brief Generate a block diagonal matrix.
 * \tparam NBlock Number of the block of the matrix.
 * \param blockMat Block matrix to repeat.
 * \return Block diagonal matrix.
 */
template <int NBlock>
Eigen::MatrixXd makeDiag(const Eigen::MatrixXd& blockMat)
{
    static_assert(N >= 0, "Not yet ready for dynamic");
    constexpr int N = NBlock;
    Eigen::MatrixXd out = Eigen::MatrixXd::Zero(N * blockMat.rows(), N * blockMat.cols());
    for (int i = 0; i < N; ++i)
        out.block(i * blockMat.rows(), i * blockMat.cols(), blockMat.rows(), blockMat.cols()) = blockMat;

    return out;
}

template <typename CrossNOp>
Eigen::MatrixXd generateD(const CrossNOp& cx)
{
    constexpr int n_vec = CrossNOp::n_vec;
    using Scalar = typename CrossNOp::Scalar;
    using mat_t = Eigen::Matrix<Scalar, 6 * n_vec, 6 * n_vec>;
    using sub_mat_t = Eigen::Matrix<Scalar, 6, 6 * n_vec>;
    Eigen::MatrixXd D_N = Eigen::MatrixXd::Zero(6 * (n_vec - 1), 6 * (n_vec - 1));
    for (int i = 0; i < n_vec - 1; ++i)
        D_N.block<6, 6>(6 * i, 6 * i) = Eigen::Matrix6d::Identity() / (i + 1);
    return mat_t::Identity() + (mat_t() << sub_mat_t::Zero(), D_N * cx.dualMatrix().template topRows<6 * (n_vec - 1)>()).finished();
}

// M = I + Cd * I * C
template <int Order>
std::vector<Eigen::MatrixXd> getSubTreeInertia(const Model& m, const ModelConfig<Order>& mc)
{
    static_assert(Order > 0, "Not yet ready for dynamic");
    std::vector<Eigen::MatrixXd> M(m.nLinks(), Eigen::MatrixXd::Zero(6 * Order, 6 * Order));
    const auto& bodies = m.bodies();
    const auto& parents = m.parents();
    for (Index i = m.nLinks() - 1; i >= 0; --i) {
        M[i] += makeDiag<ord>(bodies[i].inertia().matrix());
        Index p = parents[i]; // parent
        if (p != -1) {
            auto C_p_b = mc.bodyMotions[p].inverse() * mc.bodyMotions[i];
            M[p] += C_p_b.template dualMatrix<Order>() * M[i] * C_p_b.inverse().template matrix<Order>();
        }
    }

    return M;
}

} // namespace cdm
