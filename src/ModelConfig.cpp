#include "cdm/ModelConfig.hpp"

namespace cdm {

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

    Eigen::MatrixXd alephMat = (mc.dqs * f.asDiagonal()).transpose();
    return Eigen::Map<Eigen::VectorXd>(alephMat.data(), alephMat.size());
}

} // namespace cdm
