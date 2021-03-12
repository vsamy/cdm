#include "SimpleHumanModel.hpp"
#include "macros.hpp"
#include "model_generation.hpp"
#include <cdm/FK.hpp>
#include <cdm/ID.hpp>
#include <cdm/FD.hpp>
#include <catch2/catch.hpp>
#include <rbdyn/FA.h>
#include <rbdyn/FD.h>
#include <rbdyn/FK.h>
#include <rbdyn/FV.h>
#include <rbdyn/ID.h>
#include <tuple>

struct FixedOrder {
    static constexpr int order = 5;
};

struct DynamicOrder {
    static constexpr int order = coma::Dynamic;
};

// TODO: this does not test dynamic currently
TEMPLATE_TEST_CASE("FD", "[FD]", FixedOrder, DynamicOrder)
{
    using Index = cdm::Index;
    constexpr int order = FixedOrder::order;

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

    std::vector<Eigen::VectorXd> tauP(model.nLinks());
    std::vector<Eigen::VectorXd> tauF(model.nLinks());
    const auto& factors = coma::factorial_factors<double, order>;
    for (Index i = 0; i < model.nLinks(); ++i) {
        auto GT = cdm::DiMotionSubspace<order>{ model.joint(i).S().transpose() };
        tauP[i] = GT * mcNoGrav.jointMomentums[i];
        tauF[i] = mc.jointTorques[i];
        Index dof = model.joint(i).dof();
        for (int n = 0; n < order; ++n) {
            tauP[i].segment(n * dof, dof) /= factors[n];
        }
    }

    // RBDyn FD
    Eigen::VectorXd tauF1{ mb.nrDof() };
    for (int i = 0; i < mb.nrJoints(); ++i) {
        int dof = mb.joint(i).dof();
        tauF1.segment(mb.jointPosInDof(i), dof) = tauF[i].segment(dof, dof);
    }
    mbc.jointTorque = rbd::vectorToDof(mb, tauF1);
    rbd::ForwardDynamics fd{ mb };
    fd.forwardDynamics(mb, mbc);

    // Prepare recursive
    Eigen::MatrixXd dqs = data.dqs[t];
    Eigen::VectorXd tau(model.nDof());
    Eigen::VectorXd f = Eigen::VectorXd::Zero(model.nDof());
    data.dqs[t].setZero();

    // FD recursion
    Init(data, model, mc);
    cdm::FK(model, mc);
    cdm::ID(model, mc);
    for (int n = 0; n < order; ++n) {
        for (Index j = 0; j < model.nLinks(); ++j) {
            Index dof = model.joint(j).dof();
            tau.segment(model.jointPosInDof(j), dof) = tauF[j].segment(n * dof, dof);
            f.segment(model.jointPosInDof(j), dof) = mc.jointTorques[j].segment(n * dof, dof);
        }
        data.dqs[t].col(n) = cdm::standardFD<order>(model, mc, tau - f);
        Init(data, model, mc);
        cdm::FK(model, mc);
        cdm::ID(model, mc);
    }

    REQUIRE(data.dqs[t].isApprox(dqs));

    // FD
    std::vector<Eigen::VectorXd> y = FD(model, mc, tauP);
    Eigen::VectorXd ddq{ model.nDof() };
    for (Index i = 0; i < model.nLinks(); ++i) {
        Index pos = model.jointPosInDof(i);
        Index dof = model.joint(i).dof();
        ddq.segment(pos, dof) = y[i].segment(dof, dof);
    }

    // Check with RBDyn
    Eigen::VectorXd alphaD = rbd::dofToVector(mb, mbc.alphaD);
    REQUIRE(alphaD.isApprox(ddq));

    // Check motion
    for (Index i = 0; i < model.nLinks(); ++i) {
        Index dof = model.joint(i).dof();
        for (int n = 0; n < order; ++n) {
            REQUIRE(dqs.col(n).segment(model.jointPosInDof(i), dof).isApprox(y[i].segment(n * dof, dof) * factors[n]));
        }
    }
}