/*
 * Copyright 2020-2021 CNRS-UM LIRMM, CNRS-AIST JRL
 */

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
    // TODO: Bench with inverse
    for (Index i = 0; i < m.nLinks(); ++i) {
        size_t ui = static_cast<size_t>(i);
        Index parent = m.jointParent(i);
        if (parent != -1) {
            mc.bodyMotions[ui] = mc.bodyMotions[static_cast<size_t>(parent)] * mc.jointMotions[ui];
        } else {
            mc.bodyMotions[ui] = mc.world * mc.jointMotions[ui];
        }
    }
}

} // namespace cdm
