/*
 * Copyright 2020-2021 CNRS-UM LIRMM, CNRS-AIST JRL
 */

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
    Eigen::VectorXd getAleph(const Model& m) const;

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
    , dqs(m.nDof(), order)
    , jointMotions(static_cast<size_t>(m.nLinks()))
    , bodyMotions(static_cast<size_t>(m.nLinks()))
    , jointMomentums(static_cast<size_t>(m.nLinks()))
    , bodyMomentums(static_cast<size_t>(m.nLinks()))
    , jointForces(static_cast<size_t>(m.nLinks()))
    , bodyForces(static_cast<size_t>(m.nLinks()))
    , jointTorques(static_cast<size_t>(m.nLinks()))
{
    size_t nLinks = static_cast<size_t>(m.nLinks());
    for (size_t i = 0; i < nLinks; ++i) {
        jointTorques[i].resize(m.nDof());
    }
}

template <int Order>
void ModelConfig<Order>::setZero(const Model& m)
{
    q.setZero();
    dqs.setZero();
    size_t nLinks = static_cast<size_t>(m.nLinks());
    for (size_t i = 0; i < nLinks; ++i) {
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
Eigen::VectorXd ModelConfig<Order>::getAleph(const Model& m) const
{
    auto order = dqs.cols();
    const auto& jointPosInDof = m.jointPosInDof();
    const auto& factors = coma::factorial_factors<double, Order>;
    Eigen::VectorXd v(order * m.nDof());
    Index curOrderPos = 0;
    for (Index i = 0; i < m.nLinks(); ++i) {
        Index dof = m.joint(i).dof();
        for (Index k = 0; k < order; ++k) {
            v.segment(curOrderPos + k * dof, dof) = dqs.col(k).segment(jointPosInDof[static_cast<size_t>(i)], dof) / factors[static_cast<size_t>(k)];
        }
        curOrderPos += order * dof;
    }

    return v;
}

} // namespace cdm
