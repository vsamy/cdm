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
    static_assert(NBlock >= 0, "Not yet ready for dynamic");
    constexpr int N = NBlock;
    Eigen::MatrixXd out = Eigen::MatrixXd::Zero(N * blockMat.rows(), N * blockMat.cols());
    for (int i = 0; i < N; ++i)
        out.block(i * blockMat.rows(), i * blockMat.cols(), blockMat.rows(), blockMat.cols()) = blockMat;

    return out;
}

template <int NVec>
Eigen::MatrixXd generateD(const CrossN<NVec>& cx)
{
    using mat_t = Eigen::Matrix<double, 6 * NVec, 6 * NVec>;
    using sub_mat_t = Eigen::Matrix<double, 6, 6 * NVec>;
    Eigen::MatrixXd D_N = Eigen::MatrixXd::Zero(6 * (NVec - 1), 6 * (NVec - 1));
    for (int i = 0; i < NVec - 1; ++i)
        D_N.block<6, 6>(6 * i, 6 * i) = Eigen::Matrix6d::Identity() / (i + 1);
    return mat_t::Identity() + (mat_t() << sub_mat_t::Zero(), D_N * cx.dualMatrix().template topRows<6 * (NVec - 1)>()).finished();
}

// M = I + Cd * I * C
template <int Order>
std::vector<Eigen::MatrixXd> getSubTreeInertia(const Model& m, const ModelConfig<Order>& mc)
{
    static_assert(Order > 0, "Not yet ready for dynamic");
    std::vector<Eigen::MatrixXd> M(static_cast<size_t>(m.nLinks()), Eigen::MatrixXd::Zero(6 * Order, 6 * Order));
    for (Index i = m.nLinks() - 1; i >= 0; --i) {
        size_t ui = static_cast<size_t>(i);
        M[ui] += makeDiag<Order>(m.body(i).inertia().matrix());
        Index p = m.jointParent(i);
        if (p != -1) {
            size_t up = static_cast<size_t>(p);
            auto C_p_b = mc.bodyMotions[up].inverse() * mc.bodyMotions[ui];
            M[up] += C_p_b.template dualMatrix<Order>() * M[ui] * C_p_b.inverse().template matrix<Order>();
        }
    }

    return M;
}

} // namespace cdm
