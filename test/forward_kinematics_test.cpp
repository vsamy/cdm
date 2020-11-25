#include "SimpleHumanModel.hpp"
#include "macros.hpp"
#include "model_generation.hpp"
#include <catch2/catch.hpp>
#include <cdm/Core>
#include <rbdyn/FA.h>
#include <rbdyn/FK.h>
#include <rbdyn/FV.h>
#include <tuple>

struct FixedOrder {
    static constexpr int order = 5;
};

struct DynamicOrder {
    static constexpr int order = coma::Dynamic;
};

TEMPLATE_TEST_CASE("FK", "[FK]", FixedOrder, DynamicOrder)
{
    constexpr int order = TestType::order;

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
    FK(model, mc1);

    // RBDyn FK
    Init(data, mb, mbc);
    rbd::forwardKinematics(mb, mbc);
    rbd::forwardVelocity(mb, mbc);
    rbd::forwardAcceleration(mb, mbc);

    // Check joint motion first
    for (cdm::Index i = 0; i < mb.nrJoints(); ++i) {
        REQUIRE(mbc.parentToSon[i].translation().isApprox(mc1.jointMotions[i].transform().translation()));
        REQUIRE(mbc.parentToSon[i].rotation().transpose().isApprox(mc1.jointMotions[i].transform().rotation()));
        REQUIRE(mbc.jointVelocity[i].vector() == mc1.jointMotions[i].motion()[0].vector());
    }

    // Check FK with RBDyn
    auto linkWorldInv = mc1.world.inverse();
    double prec = coma::dummy_precision<double>();
    for (cdm::Index i = 0; i < mb.nrJoints(); ++i) {
        auto linkMotions = (linkWorldInv * mc1.bodyMotions[i]).motion();
        REQUIRE((mbc.bodyPosW[i].translation() - mc1.bodyMotions[i].transform().translation()).norm() < prec);
        REQUIRE((mbc.bodyPosW[i].rotation().transpose() - mc1.bodyMotions[i].transform().rotation()).norm() < prec);
        REQUIRE((mbc.bodyVelB[i].vector() - linkMotions[0].vector()).norm() < prec);
        REQUIRE((mbc.bodyAccB[i].vector() - linkMotions[1].vector()).norm() < prec);
    }

    // Second FK
    data.setCurData(t2);
    Init(data, model, mc2);
    FK(model, mc2);

    // Numerical check
    for (cdm::Index i = 0; i < mb.nrBodies(); ++i) {
        auto m1 = (linkWorldInv * mc1.bodyMotions[i]).motion();
        auto m2 = (linkWorldInv * mc2.bodyMotions[i]).motion();
        auto dV = (m2 - m1) / dt;
        for (cdm::Index n = 0; n < order - 1; ++n) {
            REQUIRE((dV[n] - m2[n + 1]).vector().norm() < dt * 1000.);
        }
    }
}