#include "SimpleHumanModel.hpp"
#include "macros.hpp"
#include "model_generation.hpp"
#include <catch2/catch.hpp>
#include <cdm/Core>
#include <rbdyn/FA.h>
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

TEMPLATE_TEST_CASE("ID", "[ID]", FixedOrder, DynamicOrder)
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

    // First ID used for numerical comparison
    Init(data, model, mc1);
    FK(model, mc1);
    ID(model, mc1);

    // RBDyn ID
    Init(data, mb, mbc);
    rbd::InverseDynamics id{ mb };
    rbd::forwardKinematics(mb, mbc);
    rbd::forwardVelocity(mb, mbc);
    rbd::forwardAcceleration(mb, mbc);
    id.inverseDynamics(mb, mbc);

    // Check with RBDyn
    Eigen::VectorXd tau = rbd::dofToVector(mb, mbc.jointTorque);
    double prec = coma::dummy_precision<double>();
    for (int i = 0; i < mb.nrBodies(); ++i) {
        int posInDof = mb.jointPosInDof(i);
        int dof = mb.joint(i).dof();
        // First momentum and force (=momentum)
        auto f1 = mb.body(i).inertia() * mbc.bodyVelB[i];
        REQUIRE((f1.vector() - mc1.bodyMomentums[i][0].vector()).norm() < prec);
        REQUIRE((f1.vector() - mc1.bodyForces[i][0].vector()).norm() < prec);
        // Momentum derivative
        auto f2 = mb.body(i).inertia() * mbc.bodyAccB[i];
        REQUIRE((f2.vector() - mc1.bodyMomentums[i][1].vector()).norm() < prec);
        // force
        f2 += mbc.bodyVelB[i].crossDual(f1);
        REQUIRE((f2.vector() - mc1.bodyForces[i][1].vector()).norm() < prec);
        REQUIRE((id.f()[i].vector() - mc1.jointForces[i][1].vector()).norm() < prec);
        // Torque
        REQUIRE((tau.segment(posInDof, dof) - mc1.jointTorques[i].segment(dof, dof)).norm() < prec);
    }

    // Second ID
    data.setCurData(t2);
    Init(data, model, mc2);
    FK(model, mc2);
    ID(model, mc2);

    // Numerical check
    // Please remember that here the 0-order force is the momentum and f[1] corresponds to the classical force so f[0] == p[0]
    // And we have df[0]/dt == dp[0]/dt == p[1], thus df[0]/dt != f[1] (the numerical derivative of f[0] is p[1] and not f[1])
    for (int i = 0; i < mb.nrBodies(); ++i) {
        auto dlP = (mc2.bodyMomentums[i] - mc1.bodyMomentums[i]) / dt;
        auto dlF = (mc2.bodyForces[i] - mc1.bodyForces[i]) / dt;
        auto djP = (mc2.jointMomentums[i] - mc1.jointMomentums[i]) / dt;
        auto djF = (mc2.jointForces[i] - mc1.jointForces[i]) / dt;
        for (int n = 0; n < order - 1; ++n) {
            REQUIRE((dlP[n] - mc2.bodyMomentums[i][n + 1]).vector().norm() < dt * 1000.);
            REQUIRE((djP[n] - mc2.jointMomentums[i][n + 1]).vector().norm() < dt * 1000.);
            if (n != 0) {
                REQUIRE((dlF[n] - mc2.bodyForces[i][n + 1]).vector().norm() < dt * 1000.);
                REQUIRE((djF[n] - mc2.jointForces[i][n + 1]).vector().norm() < dt * 1000.);
            } else {
                REQUIRE((dlF[n] - mc2.bodyMomentums[i][n + 1]).vector().norm() < dt * 1000.);
                REQUIRE((djF[n] - mc2.jointMomentums[i][n + 1]).vector().norm() < dt * 1000.);
            }
        }
    }
}