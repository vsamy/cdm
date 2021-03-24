#pragma once

#include "cdm/Model.hpp"
#include "cdm/typedefs.hpp"

namespace cdm {

/*! \brief Contains all model information.
 * \tparam Order Order of the model or coma::Dynamic.
 */
template <int Order>
struct ModelConfig {
    ModelConfig() = default;
    ModelConfig(const Model& m, Index order = Order);

    void setZero(const Model& m);
    Eigen::VectorXd getAleph() const;

    CMTM<Order> world; /*!< World transformation \f$C_0\f$. */
    Eigen::VectorXd q; /*!< Vector of generalized coordinates */
    Eigen::MatrixXd dqs; /*!< Concatenation of vector of generalized motions */
    std::vector<CMTM<Order>> jointMotions; /*!< Comprehensive joint motions */
    std::vector<CMTM<Order>> bodyMotions; /*!< Comprehensive body motions */
    std::vector<ForceVectorX<Order>> jointMomentums; /*!< Comprehensive joint momentums */
    std::vector<ForceVectorX<Order>> bodyMomentums; /*!< Comprehensive body momentums */
    std::vector<ForceVectorX<Order>> jointForces; /*!< Comprehensive joint forces */
    std::vector<ForceVectorX<Order>> bodyForces; /*!< Comprehensive body forces */
    std::vector<Eigen::VectorXd> jointTorques; /*!< Comprehensive joint torques */
};

template <int Order>
ModelConfig<Order>::ModelConfig(const Model& m, Index order)
    : q(m.nParam())
    , dqs(m.nDof(), order * m.nLinks())
    , jointMotions(m.nLinks())
    , bodyMotions(m.nLinks())
    , jointMomentums(m.nLinks())
    , bodyMomentums(m.nLinks())
    , jointForces(m.nLinks())
    , bodyForces(m.nLinks())
    , jointTorques(m.nLinks())
{
    for (Index i = 0; i < m.nLinks(); ++i) {
        jointTorques[i].resize(m.nDof());
    }
}

template <int Order>
void ModelConfig<Order>::setZero(const Model& m)
{
    q.setZero();
    dqs.setZero();
    for (Index i = 0; i < m.nLinks(); ++i) {
        jointMotions[i].setZero();
        bodyMotions[i].setZero();
        jointMomentums[i].setZero();
        bodyMomentums[i].setZero();
        jointForces[i].setZero();
        bodyForces[i].setZero();
        jointTorques[i].setZero();
    }
}

template <int Order>
Eigen::VectorXd ModelConfig<Order>::getAleph() const
{
    const Index order = static_cast<Index>(dqs.cols());
    const auto& factors = coma::factorial_factors<double, Order>; // TODO: Use inverse_factorial_factors instead
    Eigen::VectorXd f(order);
    for (Index i = 0; i < order; ++i) {
        f(i) = 1. / factors[i];
    }

    Eigen::MatrixXd alephMat = (dqs * f.asDiagonal()).transpose();
    return Eigen::Map<Eigen::VectorXd>(alephMat.data(), alephMat.size());
}

} // namespace cdm
