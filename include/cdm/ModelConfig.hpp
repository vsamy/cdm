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

    void setZero(const Model& m, Index order = Order);
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

} // namespace cdm
