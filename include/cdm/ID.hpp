#pragma once

#include "cdm/Model.hpp"
#include "cdm/ModelConfig.hpp"

namespace cdm {

/*! \brief Comprehensive Inverse Dynamics.
 * Compute comprehensive joint forces and torques from motions and model.
 * \tparam Order Order of the model.
 * \param m Model to compute the ID from.
 * \param mc Model information where to save the results.
 * \warning Index 0 of comprehensive joint force is set to the joint momentum.
 * Thus its derivative is not the index 1.
 */
template <int Order>
void ID(const Model& m, ModelConfig<Order>& mc)
{
    const Index order = mc.world.order();
    const auto& bodies = m.bodies();
    const auto& joints = m.joints();

    // Forward
    for (Index i = 0; i < m.nLinks(); ++i) {
        size_t ui = static_cast<size_t>(i);
        size_t uc = static_cast<size_t>(m.jointChild(i));
        const auto& motion = mc.bodyMotions[uc].motion();
        mc.bodyMomentums[uc] = DiInertia<Order>{ bodies[uc].inertia() } * motion;
        auto cnx = CrossN<Order>{ motion };
        mc.bodyForces[uc] = mc.bodyMomentums[uc];
        auto df = cnx.dualMul(mc.bodyMomentums[uc]);
        for (Index j = 0; j < order - 1; ++j) {
            mc.bodyForces[uc][j + 1] += df[j];
        }
        mc.jointMomentums[ui].setZero();
        mc.jointForces[ui].setZero();
    }

    // Backward
    for (Index i = m.nLinks() - 1; i >= 0; --i) {
        size_t ui = static_cast<size_t>(i);
        size_t uc = static_cast<size_t>(m.jointChild(i));
        mc.jointMomentums[ui] += mc.bodyMomentums[uc];
        mc.jointForces[ui] += mc.bodyForces[uc];
        mc.jointForces[ui][0] = mc.jointMomentums[ui][0];
        auto GT = DiMotionSubspace<Order>{ joints[ui].S().transpose() };
        mc.jointTorques[ui] = GT * mc.jointForces[ui];
        Index parent = m.jointParent(i);
        if (parent < 0) {
            continue;
        }

        size_t uParent = static_cast<size_t>(parent);
        mc.jointMomentums[uParent] += mc.jointMotions[ui].dualMul(mc.jointMomentums[ui]);
        if constexpr (Order == coma::Dynamic) {
            auto df = mc.jointMotions[ui].dualMul(mc.jointForces[ui].subVector(1, order - 1));
            for (Index j = 0; j < order - 1; ++j)
                mc.jointForces[uParent][j + 1] += df[j];
        } else {
            auto df = mc.jointMotions[ui].dualMul(mc.jointForces[ui].template subVector<Order - 1>(1));
            for (Index j = 0; j < order - 1; ++j)
                mc.jointForces[uParent][j + 1] += df[j];
        }
    }
}

} // namespace cdm
