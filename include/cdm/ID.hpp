#pragma once

#include "cdm/Model.hpp"
#include "cdm/ModelConfig.hpp"

namespace cdm {

template <int Order>
void ID(const Model& m, ModelConfig<Order>& mc)
{
    const Index order = mc.world.order();
    const auto& parents = m.jointParents();
    const auto& childs = m.jointChilds();
    const auto& bodies = m.bodies();
    const auto& joints = m.joints();

    // Forward
    for (Index i = 0; i < m.nLinks(); ++i) {
        const auto& motion = mc.bodyMotions[childs[i]].motion();
        mc.bodyMomentums[childs[i]] = DiInertia<Order>{ bodies[childs[i]].inertia() } * motion;
        auto cnx = CrossN<Order>{ motion };
        mc.bodyForces[childs[i]] = mc.bodyMomentums[childs[i]];
        auto df = cnx.dualMul(mc.bodyMomentums[childs[i]]);
        for (Index j = 0; j < order - 1; ++j)
            mc.bodyForces[childs[i]][j + 1] += df[j];
        mc.jointMomentums[i].setZero();
        mc.jointForces[i].setZero();
    }

    // Backward
    for (Index i = m.nLinks() - 1; i >= 0; --i) {
        mc.jointMomentums[i] += mc.bodyMomentums[childs[i]];
        mc.jointForces[i] += mc.bodyForces[childs[i]];
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
