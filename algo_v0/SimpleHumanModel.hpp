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

void createLeg(rbd::MultiBodyGraph& mbg,
    const Eigen::Vector3d& direction,
    const std::string& parentName,
    const std::string& prefix)
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

void createArm(rbd::MultiBodyGraph& mbg,
    const Eigen::Vector3d& direction,
    const std::string& parentName,
    const std::string& prefix)
{
    using namespace Eigen;
    using namespace sva;
    using namespace rbd;

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

std::tuple<rbd::MultiBody, rbd::MultiBodyConfig, rbd::MultiBodyGraph> makeHumanBody()
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
    mbg.linkBodies("HEAD1", to, "HEAD1", from, "HEAD0_HEAD1");

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
