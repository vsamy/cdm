#include "cod/FK.hpp"

namespace cod {

void FK(const Model& m, ModelConfig& mc)
{
    const auto& parents = m.jointParents();
    const auto& childs = m.jointChilds();
    for (int i = 0; i < m.nLinks(); ++i) {
        if (parents[i] != -1)
            mc.bodyMotions[childs[i]] = mc.bodyMotions[parents[i]] * mc.jointMotions[i];
        else
            mc.bodyMotions[childs[i]] = mc.world * mc.jointMotions[i];
    }
}

} // namespace cod
