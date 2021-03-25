#include "SimpleHumanModel.hpp"
#include "doctest/doctest.h"
#include "macros.hpp"
#include "model_generation.hpp"
#include <cdm/FK.hpp>
#include <RBDyn/FA.h>
#include <RBDyn/FK.h>
#include <RBDyn/FV.h>
#include <tuple>

struct FixedOrder {
    static constexpr int order = 5;
};

struct DynamicOrder {
    static constexpr int order = coma::Dynamic;
};

TEST_CASE_TEMPLATE("FK", T, FixedOrder, DynamicOrder)
{
    constexpr int order = T::order;

    rbd::MultiBody mb;
    rbd::MultiBodyConfig mbc;
    rbd::MultiBodyGraph mbg;
    std::tie(mb, mbc, mbg) = rbd::makeHumanBody();

    cdm::Model model = cdm::makeHumanBody();
    cdm::ModelConfig<order> mc1;
    cdm::ModelConfig<order> mc2;

    int nt = 21;
    double dt = 1e-8;
    int t1 = nt / 2;
    int t2 = t1 + 1;
    auto data = GenerateData<order>(mb, mbc, nt, dt);
    data.setCurData(t1);

    // First FK used for numerical comparison
    Init(data, model, mc1);
    cdm::FK(model, mc1);

    // RBDyn FK
    Init(data, mb, mbc);
    rbd::forwardKinematics(mb, mbc);
    rbd::forwardVelocity(mb, mbc);
    rbd::forwardAcceleration(mb, mbc);

    // Check joint motion first
    for (cdm::Index i = 0; i < mb.nrJoints(); ++i) {
        size_t ui = static_cast<size_t>(i);
        REQUIRE(mbc.parentToSon[ui].translation().isApprox(mc1.jointMotions[ui].transform().translation()));
        REQUIRE(mbc.parentToSon[ui].rotation().transpose().isApprox(mc1.jointMotions[ui].transform().rotation()));
        REQUIRE(mbc.jointVelocity[ui].vector() == mc1.jointMotions[ui].motion()[0].vector());
    }

    // Check FK with RBDyn
    auto linkWorldInv = mc1.world.inverse();
    double prec = coma::dummy_precision<double>();
    for (cdm::Index i = 0; i < mb.nrJoints(); ++i) {
        size_t ui = static_cast<size_t>(i);
        auto linkMotions = (linkWorldInv * mc1.bodyMotions[ui]).motion();
        REQUIRE((mbc.bodyPosW[ui].translation() - mc1.bodyMotions[ui].transform().translation()).norm() < prec);
        REQUIRE((mbc.bodyPosW[ui].rotation().transpose() - mc1.bodyMotions[ui].transform().rotation()).norm() < prec);
        REQUIRE((mbc.bodyVelB[ui].vector() - linkMotions[0].vector()).norm() < prec);
        REQUIRE((mbc.bodyAccB[ui].vector() - linkMotions[1].vector()).norm() < prec);
    }

    // Second FK
    data.setCurData(t2);
    Init(data, model, mc2);
    cdm::FK(model, mc2);

    // Numerical check
    for (cdm::Index i = 0; i < mb.nrBodies(); ++i) {
        size_t ui = static_cast<size_t>(i);
        auto m1 = (linkWorldInv * mc1.bodyMotions[ui]).motion();
        auto m2 = (linkWorldInv * mc2.bodyMotions[ui]).motion();
        auto dV = (m2 - m1) / dt;
        for (cdm::Index n = 0; n < order - 1; ++n) {
            REQUIRE((dV[n] - m2[n + 1]).vector().norm() < dt * 1000.);
        }
    }
}