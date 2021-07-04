/*
 * Copyright 2020-2021 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include "SimpleHumanModel.hpp"
#include "doctest/doctest.h"
#include "macros.hpp"
#include "model_generation.hpp"
#include <cdm/FD.hpp>
#include <cdm/FK.hpp>
#include <cdm/ID.hpp>
#include <RBDyn/FA.h>
#include <RBDyn/FD.h>
#include <RBDyn/FK.h>
#include <RBDyn/FV.h>
#include <RBDyn/ID.h>
#include <tuple>

struct FixedOrder {
    static constexpr int order = 5;
};

struct DynamicOrder {
    static constexpr int order = coma::Dynamic;
};

// TODO: this does not test dynamic currently
TEST_CASE_TEMPLATE("FD", T, FixedOrder) //, DynamicOrder)
{
    using Index = cdm::Index;
    constexpr int order = T::order;

    rbd::MultiBody mb;
    rbd::MultiBodyConfig mbc;
    rbd::MultiBodyGraph mbg;
    std::tie(mb, mbc, mbg) = rbd::makeHumanBody();

    cdm::Model model = cdm::makeHumanBody();
    cdm::ModelConfig<order> mcNoGrav;
    cdm::ModelConfig<order> mc;

    int nt = 21;
    double dt = 1e-8;
    int t = nt / 2;
    auto data = GenerateData<order>(mb, mbc, nt, dt);
    data.setCurData(t);

    // No gravity
    data.gravity.setZero();
    Init(data, model, mcNoGrav);
    cdm::FK(model, mcNoGrav);
    cdm::ID(model, mcNoGrav);

    // With gravity
    data.gravity = Eigen::Vector3d(0, 0, 9.81); // gravity acceleration applied on robot
    Init(data, model, mc);
    cdm::FK(model, mc);
    cdm::ID(model, mc);

    // RBDyn
    Init(data, mb, mbc);
    rbd::InverseDynamics id{ mb };
    rbd::forwardKinematics(mb, mbc);
    rbd::forwardVelocity(mb, mbc);
    rbd::forwardAcceleration(mb, mbc);
    id.inverseDynamics(mb, mbc);

    size_t un = static_cast<size_t>(model.nLinks());
    std::vector<Eigen::VectorXd> tauP(un);
    std::vector<Eigen::VectorXd> tauF(un);
    const auto& factors = coma::factorial_factors<double, order>;
    for (Index i = 0; i < model.nLinks(); ++i) {
        size_t ui = static_cast<size_t>(i);
        auto GT = cdm::DiMotionSubspace<order>{ model.joint(i).S().transpose() };
        tauP[ui] = GT * mcNoGrav.jointMomentums[ui];
        tauF[ui] = mc.jointTorques[ui];
        Index dof = model.joint(i).dof();
        for (int n = 0; n < order; ++n) {
            tauP[ui].segment(n * dof, dof) /= factors[static_cast<size_t>(n)];
        }
    }

    // RBDyn FD
    Eigen::VectorXd tauF1{ mb.nrDof() };
    for (int i = 0; i < mb.nrJoints(); ++i) {
        int dof = mb.joint(i).dof();
        tauF1.segment(mb.jointPosInDof(i), dof) = tauF[static_cast<size_t>(i)].segment(dof, dof);
    }
    mbc.jointTorque = rbd::vectorToDof(mb, tauF1);
    rbd::ForwardDynamics fd{ mb };
    fd.forwardDynamics(mb, mbc);

    // Prepare recursive
    Eigen::MatrixXd dqs = data.dqs[static_cast<size_t>(t)];
    Eigen::VectorXd tau(model.nDof());
    Eigen::VectorXd f = Eigen::VectorXd::Zero(model.nDof());
    data.dqs[static_cast<size_t>(t)].setZero();

    // FD recursion
    Init(data, model, mc);
    cdm::FK(model, mc);
    cdm::ID(model, mc);
    for (int n = 0; n < order; ++n) {
        for (Index j = 0; j < model.nLinks(); ++j) {
            size_t uj = static_cast<size_t>(j);
            Index dof = model.joint(j).dof();
            tau.segment(model.jointPosInDof(j), dof) = tauF[uj].segment(n * dof, dof);
            f.segment(model.jointPosInDof(j), dof) = mc.jointTorques[uj].segment(n * dof, dof);
        }
        data.dqs[static_cast<size_t>(t)].col(n) = cdm::standardFD<order>(model, mc, tau - f);
        Init(data, model, mc);
        cdm::FK(model, mc);
        cdm::ID(model, mc);
    }

    REQUIRE(data.dqs[static_cast<size_t>(t)].isApprox(dqs));

    // FD
    std::vector<Eigen::VectorXd> y = FD(model, mc, tauP);
    Eigen::VectorXd ddq{ model.nDof() };
    for (Index i = 0; i < model.nLinks(); ++i) {
        Index pos = model.jointPosInDof(i);
        Index dof = model.joint(i).dof();
        ddq.segment(pos, dof) = y[static_cast<size_t>(i)].segment(dof, dof);
    }

    // Check with RBDyn
    Eigen::VectorXd alphaD = rbd::dofToVector(mb, mbc.alphaD);
    REQUIRE(alphaD.isApprox(ddq));

    // Check motion
    for (Index i = 0; i < model.nLinks(); ++i) {
        Index dof = model.joint(i).dof();
        for (int n = 0; n < order; ++n) {
            Eigen::VectorXd mot = dqs.col(n).segment(model.jointPosInDof(i), dof);
            Eigen::VectorXd yMot = y[static_cast<size_t>(i)].segment(n * dof, dof) * factors[static_cast<size_t>(n)];
            REQUIRE(mot.size() == yMot.size());
            for (Index k = 0; k < mot.size(); ++k) {
                REQUIRE(std::abs(mot(k) - yMot(k)) < coma::dummy_precision<double>());
            }
            // REQUIRE(mot.isApprox(yMot)); // isApprox does not output good rsults for some reasons
        }
    }
}