#pragma once

#include "utils.hpp"

template <int Space, int Order>
struct Tree {
    static constexpr int space = Space;
    static constexpr int order = Order;
    using Transf = coma::Transform<double>;
    using CMTM = coma::CMTM<double, Space, Order>;
    CMTM linkWorld;
    std::vector<CMTM> joints;
    std::vector<CMTM> links;
    std::vector<coma::SpatialInertia<double>> linkInertias;
    std::vector<coma::ForceVectorX<double, Order>> jointMomentums;
    std::vector<coma::ForceVectorX<double, Order>> linkMomentums;
    std::vector<coma::ForceVectorX<double, Order>> jointForces;
    std::vector<coma::ForceVectorX<double, Order>> linkForces;
    std::vector<Eigen::VectorXd> jointTorques;

    template <typename MI>
    void init(int t, MI& info)
    {
        const auto& mb = info.model.mb;
        auto& mbc = info.model.mbc;

        links.resize(mb.nrBodies());
        joints.resize(mb.nrJoints());
        linkInertias.resize(mb.nrBodies());
        jointMomentums.resize(mb.nrJoints());
        jointForces.resize(mb.nrJoints());
        linkMomentums.resize(mb.nrBodies());
        linkForces.resize(mb.nrBodies());
        jointTorques.resize(mb.nrBodies());

        mbc.q = rbd::vectorToParam(mb, info.q.col(t));
        mbc.alpha = rbd::vectorToDof(mb, info.dqs[0].col(t));
        mbc.alphaD = rbd::vectorToDof(mb, info.dqs[1].col(t));

        const auto& Xt = mb.transforms();
        int pos = 0;
        for (int i = 0; i < mb.nrJoints(); ++i) {
            const auto& j = mb.joint(i);
            const auto& S = j.motionSubspace();
            const auto& inertia = mb.body(i).inertia();
            linkInertias[i] = coma::SpatialInertia<double>{ inertia.mass(), inertia.momentum(), inertia.inertia() };
            auto X = j.pose(mbc.q[i]) * Xt[i];
            CMTMSet<Transf, CMTM, order>(Transf { X.rotation().transpose(), X.translation() }, S, info.dqs, t, pos, j.dof(), joints[i]);
            pos += j.dof();
        }

        auto v = coma::MotionVectorNd<order>::Zero();
        if constexpr (order >= 2) {
            v[1] = coma::MotionVectord(Eigen::Vector3d::Zero(), mbc.gravity);
        }
        linkWorld.set(coma::Transform<double>::Identity(), v);
    }
};

template <size_t Order>
using Tree4d = Tree<4, Order>;
template <size_t Order>
using Tree6d = Tree<6, Order>;
