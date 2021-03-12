#pragma once

#include "cdm/Model.hpp"
#include "cdm/ModelConfig.hpp"

namespace cdm {

template <int Order>
void FK(const Model& m, ModelConfig<Order>& mc)
{
    const auto& parents = m.jointParents();
    const auto& childs = m.jointChilds();
    // TODO: Bench with inverse
    for (Index i = 0; i < m.nLinks(); ++i) {
        if (parents[i] != -1)
            mc.bodyMotions[childs[i]] = mc.bodyMotions[parents[i]] * mc.jointMotions[i];
        else
            mc.bodyMotions[childs[i]] = mc.world * mc.jointMotions[i];
    }
}

} // namespace cdm
