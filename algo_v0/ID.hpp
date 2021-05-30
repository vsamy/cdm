#pragma once

#include "FK.hpp"
#include "utils.hpp"

template <typename MI, typename Tree>
void ID(MI& info, Tree& tree)
{
    const auto& mb = info.model.mb;
    FK(info, tree);

    constexpr int ord = Tree::order;
    const auto& pred = mb.predecessors();
    const auto& succ = mb.successors();
    // Forward
    for (int i = 0; i < mb.nrJoints(); ++i) {
        const auto& motions = tree.links[succ[i]].motion();
        tree.linkMomentums[succ[i]] = coma::DiInertiaNd<ord>{ tree.linkInertias[succ[i]] } * motions;
        auto cnx = coma::CrossNd<ord> { motions };
        tree.linkForces[succ[i]] = tree.linkMomentums[succ[i]];
        auto df = cnx.dualMul(tree.linkMomentums[succ[i]]);
        for (int j = 0; j < ord - 1; ++j)
            tree.linkForces[succ[i]][j + 1] += df[j];
        tree.jointMomentums[i].setZero();
        tree.jointForces[i].setZero();
    }

    // Backward
    for (int i = mb.nrJoints() - 1; i >= 0; --i) {
        tree.jointMomentums[i] += tree.linkMomentums[succ[i]];
        tree.jointForces[i] += tree.linkForces[succ[i]];
        tree.jointForces[i][0] = tree.jointMomentums[i][0];
        auto G = makeDiag<ord>(mb.joint(i).motionSubspace());
        tree.jointTorques[i] = G.transpose() * tree.jointForces[i].vector();
        if (pred[i] < 0)
            continue;

        tree.jointMomentums[pred[i]] += tree.joints[i].dualMul(tree.jointMomentums[i]);
        auto df = tree.joints[i].dualMul(tree.jointForces[i].template subVector<ord - 1>(1));
        for (int j = 0; j < ord - 1; ++j)
            tree.jointForces[pred[i]][j + 1] += df[j];
    }
}