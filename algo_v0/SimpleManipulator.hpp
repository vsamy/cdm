#pragma once

// includes
// std
#include <tuple>

// RBDyn
#include <RBDyn/Body.h>
#include <RBDyn/Joint.h>
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>
#include <RBDyn/MultiBodyGraph.h>

/// @return An simple XYZ spherical arm with Y as up axis.
std::tuple<rbd::MultiBody, rbd::MultiBodyConfig, rbd::MultiBodyGraph> makeManipulator(size_t nrJoints = 5)
{
    rbd::MultiBodyGraph mbg;

    double mass = 1.;
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    Eigen::Vector3d h = Eigen::Vector3d::Zero();

    sva::RBInertiad rbi(mass, h, I);

    auto genJType = [](size_t j) {
        switch (j % 3) {
        case 0:
            return rbd::Joint::RevX;
        case 1:
            return rbd::Joint::RevY;
        default:
        case 2:
            return rbd::Joint::RevZ;
        }
    };

    sva::PTransformd to(Eigen::Vector3d(0.5, 0.5, -0.5));
    sva::PTransformd from(Eigen::Vector3d(0.5, -0.5, 0.5));
    std::string bpName = "b0";
    mbg.addBody(rbd::Body(rbi, bpName));
    for (size_t i = 0; i < nrJoints; ++i) {
        std::string bName = "b" + std::to_string(i + 1);
        std::string jName = "j" + std::to_string(i + 1);
        rbd::Body b(rbi, bName);
        rbd::Joint j(genJType(i), true, jName);
        mbg.addBody(b);
        mbg.addJoint(j);
        mbg.linkBodies(bpName, to, bName, from, jName);
        bpName = bName;
    }

    rbd::MultiBody mb = mbg.makeMultiBody("b0", true);

    rbd::MultiBodyConfig mbc(mb);
    mbc.zero(mb);

    return std::make_tuple(mb, mbc, mbg);
}
