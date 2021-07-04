/*
 * Copyright 2020-2021 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <cdm/ModelConstructor.hpp>

// includes
// std
#include <tuple>

// RBDyn
#include <RBDyn/Body.h>
#include <RBDyn/Joint.h>
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>
#include <RBDyn/MultiBodyGraph.h>

namespace rbd {

void createLeg(MultiBodyGraph& mbg, const Eigen::Vector3d& direction,
    const std::string& parentName, const std::string& prefix)
{
    using namespace Eigen;
    using namespace sva;
    using namespace rbd;

    double mass = 1.;
    Matrix3d I = Matrix3d::Identity();
    Vector3d h = Vector3d::Zero();

    RBInertiad rbi(mass, h, I);

    mbg.addBody({ rbi, prefix + "LEG0" });
    mbg.addBody({ rbi, prefix + "LEG1" });
    mbg.addBody({ rbi, prefix + "LEG2" });

    mbg.addJoint({ Joint::Spherical, true, prefix + "LEGBASE_0" });
    mbg.addJoint({ Joint::RevX, true, prefix + "LEG0_1" });
    mbg.addJoint({ Joint::Spherical, true, prefix + "LEG1_2" });

    PTransformd to(direction);
    PTransformd from(PTransformd::Identity());

    mbg.linkBodies(parentName, to, prefix + "LEG0", from, prefix + "LEGBASE_0");
    mbg.linkBodies(prefix + "LEG0", to, prefix + "LEG1", from, prefix + "LEG0_1");
    mbg.linkBodies(prefix + "LEG1", to, prefix + "LEG2", from, prefix + "LEG1_2");
}

void createArm(MultiBodyGraph& mbg, const Eigen::Vector3d& direction,
    const std::string& parentName, const std::string& prefix)
{
    using namespace Eigen;
    using namespace sva;

    double mass = 1.;
    Matrix3d I = Matrix3d::Identity();
    Vector3d h = Vector3d::Zero();

    RBInertiad rbi(mass, h, I);

    mbg.addBody({ rbi, prefix + "ARM0" });
    mbg.addBody({ rbi, prefix + "ARM1" });
    mbg.addBody({ rbi, prefix + "ARM2" });

    mbg.addJoint({ Joint::Spherical, true, prefix + "ARMBASE_0" });
    mbg.addJoint({ Joint::RevX, true, prefix + "ARM0_1" });
    mbg.addJoint({ Joint::Spherical, true, prefix + "ARM1_2" });

    PTransformd to(direction);
    PTransformd from(PTransformd::Identity());

    mbg.linkBodies(parentName, to, prefix + "ARM0", from, prefix + "ARMBASE_0");
    mbg.linkBodies(prefix + "ARM0", to, prefix + "ARM1", from, prefix + "ARM0_1");
    mbg.linkBodies(prefix + "ARM1", to, prefix + "ARM2", from, prefix + "ARM1_2");
}

std::tuple<MultiBody, MultiBodyConfig, MultiBodyGraph> makeHumanBody()
{
    using namespace Eigen;
    using namespace sva;
    using namespace rbd;

    MultiBodyGraph mbg;

    double mass = 1.;
    Matrix3d I = Matrix3d::Identity();
    Vector3d h = Vector3d::Zero();

    RBInertiad rbi(mass, h, I);

    mbg.addBody({ rbi, "BODY0" });
    mbg.addBody({ rbi, "BODY1" });
    mbg.addBody({ rbi, "TORSO" });
    mbg.addBody({ rbi, "HEAD0" });
    mbg.addBody({ rbi, "HEAD1" });

    mbg.addJoint({ Joint::RevY, true, "BODY0_BODY1" });
    mbg.addJoint({ Joint::RevX, true, "BODY1_TORSO" });
    mbg.addJoint({ Joint::RevY, true, "TORSO_HEAD0" });
    mbg.addJoint({ Joint::RevX, true, "HEAD0_HEAD1" });

    PTransformd to(Vector3d(0., 0.1, 0.));
    PTransformd from(PTransformd::Identity());

    mbg.linkBodies("BODY0", to, "BODY1", from, "BODY0_BODY1");
    mbg.linkBodies("BODY1", to, "TORSO", from, "BODY1_TORSO");
    mbg.linkBodies("TORSO", to, "HEAD0", from, "TORSO_HEAD0");
    mbg.linkBodies("HEAD0", to, "HEAD1", from, "HEAD0_HEAD1");

    // left arm
    createArm(mbg, Vector3d(-0.1, 0.05, 0.), "TORSO", "L");
    // right arm
    createArm(mbg, Vector3d(0.1, 0.05, 0.), "TORSO", "R");
    // left leg
    createLeg(mbg, Vector3d(-0.1, -0.05, 0.), "BODY0", "L");
    // right leg
    createLeg(mbg, Vector3d(0.1, -0.05, 0.), "BODY0", "R");

    MultiBody mb = mbg.makeMultiBody("BODY0", false);

    MultiBodyConfig mbc(mb);
    mbc.zero(mb);

    return std::make_tuple(mb, mbc, mbg);
}

} // namespace rbd

namespace cdm {

void createLeg(ModelConstructor& mc, const Eigen::Vector3d& direction,
    const std::string& parentName, const std::string& prefix)
{
    using m3_t = Eigen::Matrix3d;
    using v3_t = Eigen::Vector3d;

    double mass = 1.;
    m3_t I = m3_t::Identity();
    v3_t h = v3_t::Zero();

    Inertia si(mass, h, I);
    mc.addLink({ prefix + "LEGBASE_0", Joint::Type::Spherical }, { prefix + "LEG0", si });
    mc.addLink({ prefix + "LEG0_1", Joint::Type::Revolute, v3_t::UnitX() }, { prefix + "LEG1", si });
    mc.addLink({ prefix + "LEG1_2", Joint::Type::Spherical }, { prefix + "LEG2", si });

    Transform to{ I, direction };
    mc.connectLink(parentName, prefix + "LEGBASE_0", to);
    mc.connectLink(prefix + "LEG0", prefix + "LEG0_1", to);
    mc.connectLink(prefix + "LEG1", prefix + "LEG1_2", to);
}

void createArm(ModelConstructor& mc, const Eigen::Vector3d& direction,
    const std::string& parentName, const std::string& prefix)
{
    using m3_t = Eigen::Matrix3d;
    using v3_t = Eigen::Vector3d;

    double mass = 1.;
    m3_t I = m3_t::Identity();
    v3_t h = v3_t::Zero();

    Inertia si(mass, h, I);
    mc.addLink({ prefix + "ARMBASE_0", Joint::Type::Spherical }, { prefix + "ARM0", si });
    mc.addLink({ prefix + "ARM0_1", Joint::Type::Revolute, v3_t::UnitX() }, { prefix + "ARM1", si });
    mc.addLink({ prefix + "ARM1_2", Joint::Type::Spherical }, { prefix + "ARM2", si });

    Transform to{ I, direction };
    mc.connectLink(parentName, prefix + "ARMBASE_0", to);
    mc.connectLink(prefix + "ARM0", prefix + "ARM0_1", to);
    mc.connectLink(prefix + "ARM1", prefix + "ARM1_2", to);
}

Model makeHumanBody()
{
    using m3_t = Eigen::Matrix3d;
    using v3_t = Eigen::Vector3d;

    ModelConstructor mc;

    double mass = 1.;
    m3_t I = m3_t::Identity();
    v3_t h = v3_t::Zero();

    Inertia si(mass, h, I);
    mc.addLink({ "Root", Joint::Type::Free }, { "BODY0", si });
    mc.addLink({ "BODY0_BODY1", Joint::Type::Revolute, v3_t::UnitY() }, { "BODY1", si });
    mc.addLink({ "BODY1_TORSO", Joint::Type::Revolute, v3_t::UnitX() }, { "TORSO", si });
    mc.addLink({ "TORSO_HEAD0", Joint::Type::Revolute, v3_t::UnitY() }, { "HEAD0", si });
    mc.addLink({ "HEAD0_HEAD1", Joint::Type::Revolute, v3_t::UnitX() }, { "HEAD1", si });

    Transform to(I, v3_t(0., 0.1, 0.));
    mc.connectLink("BODY0", "BODY0_BODY1", to);
    mc.connectLink("BODY1", "BODY1_TORSO", to);
    mc.connectLink("TORSO", "TORSO_HEAD0", to);
    mc.connectLink("HEAD0", "HEAD0_HEAD1", to);

    // left arm
    createArm(mc, v3_t(-0.1, 0.05, 0.), "TORSO", "L");
    // right arm
    createArm(mc, v3_t(0.1, 0.05, 0.), "TORSO", "R");
    // left leg
    createLeg(mc, v3_t(-0.1, -0.05, 0.), "BODY0", "L");
    // right leg
    createLeg(mc, v3_t(0.1, -0.05, 0.), "BODY0", "R");

    return mc.build("Root");
}

} // namespace cdm