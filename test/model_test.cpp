/*
 * Copyright 2020-2021 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include "SimpleHumanModel.hpp"
#include "doctest/doctest.h"
#include "macros.hpp"
#include "model_generation.hpp"
#include <RBDyn/EulerIntegration.h>

TEST_CASE("body")
{
    using namespace cdm;
    using v3_t = Eigen::Vector3d;
    using m3_t = Eigen::Matrix3d;
    v3_t vec = v3_t::Random();
    m3_t mat = m3_t::Random();

    {
        Inertia I{ vec(0), vec, mat };
        Body b{ "name", I };
        REQUIRE(b.name() == "name");
        REQUIRE(b.inertia() == I);
    }
    {
        Body b{ "name", vec(0), vec, mat };
        REQUIRE(b.name() == "name");
        REQUIRE(b.inertia() == Inertia{ vec(0), vec, mat });
    }
    {
        Body b1{ "name", vec(0), vec, mat };
        Body b2{ "name", vec(2), vec, mat };
        REQUIRE(b1 == b2);
        REQUIRE(!(b1 != b2));
    }
}

TEST_CASE("joint")
{
    using namespace cdm;
    using v3_t = Eigen::Vector3d;

    v3_t vec = v3_t::Random();
    vec.normalize();

    {
        // Type Free
        Joint j{ "name", Joint::Type::Free };
        REQUIRE(j.S().matrix() == Eigen::Matrix6d::Identity());
    }
    {
        // Type Spherical
        Joint j{ "name", Joint::Type::Spherical };
        REQUIRE(j.S().matrix() == Eigen::Matrix<double, 6, 3>::Identity());
    }
    {
        // Type Revolute
        Joint j{ "name", Joint::Type::Revolute, vec };
        REQUIRE(j.S().matrix() == (Eigen::Vector6d() << vec, v3_t::Zero()).finished());
    }
    {
        // Type Prismatic
        Joint j{ "name", Joint::Type::Prismatic, vec };
        REQUIRE(j.S().matrix() == (Eigen::Vector6d() << v3_t::Zero(), vec).finished());
    }
    {
        // Type Fixed
        Joint j{ "name", Joint::Type::Fixed };
        REQUIRE(j.S().matrix() == Eigen::Matrix<double, 6, 0>());
    }
    {
        Joint j1{ "name", Joint::Type::Fixed };
        Joint j2{ "name", Joint::Type::Fixed };
        REQUIRE(j1 == j2);
        REQUIRE(!(j1 != j2));
    }
}

TEST_CASE("model")
{
    using namespace cdm;
    using v3_t = Eigen::Vector3d;
    using m3_t = Eigen::Matrix3d;
    v3_t vec = v3_t::Random();
    m3_t mat = m3_t::Random();
    Inertia I{ vec(0), vec, mat };

    /////////////////////////////////////////////////
    //                   /- j2|b2 --- j3b3         //
    // j0|b0 --- j1|b1 -|                          //
    //                  \- j4|b4                   //
    /////////////////////////////////////////////////
    std::vector<Body> bs;
    std::vector<Joint> js;
    std::vector<Transform> T0;
    std::vector<int> jointParents = { -1, 0, 1, 2, 1 };
    std::vector<int> jointChildren = { 0, 1, 2, 3, 4 };
    bs.emplace_back("b0", I);
    bs.emplace_back("b1", I);
    bs.emplace_back("b2", I);
    bs.emplace_back("b3", I);
    bs.emplace_back("b4", I);
    js.emplace_back("j0", Joint::Type::Free);
    js.emplace_back("j1", Joint::Type::Revolute, vec);
    js.emplace_back("j2", Joint::Type::Prismatic, vec);
    js.emplace_back("j3", Joint::Type::Fixed);
    js.emplace_back("j4", Joint::Type::Spherical);
    DISABLE_CONVERSION_WARNING_BEGIN
    T0.emplace_back(Eigen::Quaterniond::UnitRandom(), v3_t::Random());
    T0.emplace_back(Eigen::Quaterniond::UnitRandom(), v3_t::Random());
    T0.emplace_back(Eigen::Quaterniond::UnitRandom(), v3_t::Random());
    T0.emplace_back(Eigen::Quaterniond::UnitRandom(), v3_t::Random());
    Transform TRoot{ Eigen::Quaterniond::UnitRandom(), v3_t::Random() };
    DISABLE_CONVERSION_WARNING_END

    // Add in random order
    ModelConstructor mc;
    mc.addLink(js[1], bs[1]);
    mc.addLink(js[3], bs[3]);
    mc.addLink(js[2], bs[2]);
    mc.addLink(js[4], bs[4]);
    mc.addLink(js[0], bs[0]);
    mc.connectLink("b1", "j2", T0[1]);
    mc.connectLink("b2", "j3", T0[2]);
    mc.connectLink("b0", "j1", T0[0]);
    mc.connectLink("b1", "j4", T0[3]);
    // REQUIRE_THROWS_AS(mc.body("b0", "j7", TRoot), std::runtime_error);
    // REQUIRE_THROWS_AS(mc.body("b7", "j0", TRoot), std::runtime_error);
    // REQUIRE_THROWS_AS(mc.build("j7"), std::runtime_error);
    auto model = mc.build("j0", TRoot);

    // checks
    REQUIRE(model.nLinks() == 5);
    for (int i = 0; i < model.nLinks(); ++i) {
        size_t ui = static_cast<size_t>(i);
        REQUIRE(model.joint(i) == js[ui]);
        REQUIRE(model.body(i) == bs[ui]);
        if (i != 0) {
            REQUIRE(model.T0(i) == T0[ui - 1]);
        } else {
            REQUIRE(model.T0(i) == TRoot);
        }
        REQUIRE(model.jointParent(i) == jointParents[ui]);
        REQUIRE(model.jointChild(i) == jointChildren[ui]);
    }

    REQUIRE_THROWS_AS(model.jointParentAt(7), std::out_of_range);
    REQUIRE_THROWS_AS(model.jointChildAt(7), std::out_of_range);
    REQUIRE_THROWS_AS(model.jointAt(7), std::out_of_range);
    REQUIRE_THROWS_AS(model.bodyAt(7), std::out_of_range);
    REQUIRE_THROWS_AS(model.T0At(7), std::out_of_range);
}

TEST_CASE("rbd model")
{
    rbd::MultiBody mb;
    rbd::MultiBodyConfig mbc;
    rbd::MultiBodyGraph mbg;
    std::tie(mb, mbc, mbg) = rbd::makeHumanBody();
    mbc.zero(mb);

    cdm::Model model = cdm::makeHumanBody();

    // Checks
    REQUIRE(model.nLinks() == mb.nrBodies());
    REQUIRE(model.nDof() == mb.nrDof());
    for (int i = 0; i < model.nLinks(); ++i) {
        REQUIRE(model.joint(i).name() == mb.joint(i).name());
        REQUIRE(model.joint(i).dof() == mb.joint(i).dof());
        REQUIRE(model.joint(i).S() == mb.joint(i).motionSubspace());
        REQUIRE(model.body(i).name() == mb.body(i).name());
        REQUIRE(model.body(i).inertia().mass() == mb.body(i).inertia().mass());
        REQUIRE(model.body(i).inertia().momentum() == mb.body(i).inertia().momentum());
        REQUIRE(model.body(i).inertia().inertia() == mb.body(i).inertia().inertia());
        REQUIRE(model.jointPosInParam(i) == mb.jointPosInParam(i));
        REQUIRE(model.jointPosInDof(i) == mb.jointPosInDof(i));
        REQUIRE(model.T0(i).translation() == mb.transform(i).translation());
        REQUIRE(model.T0(i).rotation().isApprox(mb.transform(i).rotation().transpose()));
    }
}
