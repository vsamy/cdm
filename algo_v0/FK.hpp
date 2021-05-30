#pragma once

template <typename MI, typename Tree>
void FK(const MI& info, Tree& tree)
{
    // Forward
    const auto& mb = info.model.mb;
    const auto& pred = mb.predecessors();
    const auto& succ = mb.successors();
    for (int i = 0; i < mb.nrJoints(); ++i) {
        if (pred[i] != -1)
            tree.links[succ[i]] = tree.links[pred[i]] * tree.joints[i];
        else
            tree.links[succ[i]] = tree.linkWorld * tree.joints[i];
    }
}
