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
    const auto& parents = m.jointParents();
    const auto& childs = m.jointChilds();
    const auto& links = m.links();
    const auto& joints = m.joints();

    // Forward
    for (Index i = 0; i < m.nLinks(); ++i) {
        const auto& motion = mc.LinkMotions[childs[i]].motion();
        mc.linkMomentums[childs[i]] = DiInertia<Order>{ links[childs[i]].inertia() } * motion;
        auto cnx = CrossN<Order>{ motion };
        mc.linkForces[childs[i]] = mc.linkMomentums[childs[i]];
        auto df = cnx.dualMul(mc.linkMomentums[childs[i]]);
        for (Index j = 0; j < order - 1; ++j)
            mc.linkForces[childs[i]][j + 1] += df[j];
        mc.jointMomentums[i].setZero();
        mc.jointForces[i].setZero();
    }

    // Backward
    for (Index i = m.nLinks() - 1; i >= 0; --i) {
        mc.jointMomentums[i] += mc.linkMomentums[childs[i]];
        mc.jointForces[i] += mc.linkForces[childs[i]];
        mc.jointForces[i][0] = mc.jointMomentums[i][0];
        auto GT = DiMotionSubspace<Order>{ joints[i].S().transpose() };
        mc.jointTorques[i] = GT * mc.jointForces[i];
        if (parents[i] < 0)
            continue;

        mc.jointMomentums[parents[i]] += mc.jointMotions[i].dualMul(mc.jointMomentums[i]);
        if constexpr (Order == coma::Dynamic) {
            auto df = mc.jointMotions[i].dualMul(mc.jointForces[i].subVector(1, order - 1));
            for (Index j = 0; j < order - 1; ++j)
                mc.jointForces[parents[i]][j + 1] += df[j];
        } else {
            auto df = mc.jointMotions[i].dualMul(mc.jointForces[i].template subVector<Order - 1>(1));
            for (Index j = 0; j < order - 1; ++j)
                mc.jointForces[parents[i]][j + 1] += df[j];
        }
    }
}

} // namespace cdm
