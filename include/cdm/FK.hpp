#pragma once

#include "cdm/Model.hpp"
#include "cdm/ModelConfig.hpp"

namespace cdm {

/*! \brief Comprehensive Forward Kinematics.
 * This compute the FK accross the whole model.
 * \param m Model to compute the FK from.
 * \param mc Model information where to save the results.
 */
template <int Order>
void FK(const Model& m, ModelConfig<Order>& mc)
{
    const auto& parents = m.jointParents();
    // TODO: Bench with inverse
    for (Index i = 0; i < m.nLinks(); ++i) {
        if (parents[i] != -1) {
            mc.LinkMotions[i] = mc.LinkMotions[parents[i]] * mc.jointMotions[i];
        } else {
            mc.LinkMotions[i] = mc.world * mc.jointMotions[i];
        }
    }
}

} // namespace cdm
