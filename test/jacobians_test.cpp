#include "SimpleHumanModel.hpp"
#include "macros.hpp"
#include "model_generation.hpp"
#include <catch2/catch.hpp>
#include <cdm/Core>
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

template <int Order>
Eigen::VectorXd GetAleph(const Model& m, const ModelConfig<Order>& mc, int t)
{
    const int ord = mc.world.order();
    const auto& pos = m.jointPosInDof();
    const auto& factors = coma::factorial_factors<double, Order>;
    Eigen::VectorXd v(ord * model.mb.nrDof());
    int curOrderPos = 0;
    for (Index i = 0; i < m.nLinks(); ++i) {
        Index dof = m.joint(i).dof();
        for (size_t k = 0; k < ord; ++k) {
            v.segment(curOrderPos + k * dof, dof) = mc.dqs[t].col(k).segment(jointPosInDof[i], dof) / factors[k];
        }
        curOrderPos += ord * dof;
    }

    return v;
}

TEMPLATE_TEST_CASE("FK", "[FK]", FixedOrder, DynamicOrder)
{
    constexpr int order = TestType::order;

    rbd::MultiBody mb;
    rbd::MultiBodyConfig mbc;
    rbd::MultiBodyGraph mbg;
    std::tie(mb, mbc, mbg) = rbd::makeHumanBody();

    cdm::Model model = cdm::makeHumanBody();
    cdm::ModelConfig<order> mc;

    int nt = 21;
    double dt = 1e-8;
    int t = nt / 2;
    auto data = GenerateData<order>(mb, mbc, n, dt);
    data.setCurData(t);

    // First FK used for numerical comparison
    Init(data, model, mc);
    FK(model, mc);
    Eigen::MatrixXd BigJ = MotionJacobian(model, mc, bodyName); // dx = BigJ * \aleph

    // Extract J and dJ
    const auto& pos = m.jointPosInDof();
    for (Index i = 0; i < m.nLinks(); ++i) {
        Index dof = m.joint(i).dof();
        J.block(0, pos[i], 6, dof) = BigJ.block(0, pos[i] * order, 6, dof);
        dJ.block(0, pos[i], 6, dof) = BigJ.block(6, pos[i] * order, 6, dof);
    }

    // RBDyn FK
    Init(data, mb, mbc);
    rbd::forwardKinematics(mb, mbc);
    rbd::forwardVelocity(mb, mbc);
    rbd::forwardAcceleration(mb, mbc);

    // Get RBDyn J and dJ
    Eigen::MatrixXd J(6, mb.nrDof()); // v = J * \psi
    Eigen::MatrixXd dJ(6, mb.nrDof()); // dv = dJ * \psi + J * d\psi
    auto rbdJac = rbd::Jacobian(mb, bodyName);
    auto rbdShortJ = rbdJac.bodyJacobian(mb, mbc);
    auto rbdShortJDot = rbdJac.bodyJacobianDot(mb, mbc);
    Eigen::MatrixXd rbdJ(6, mb.nrDof());
    Eigen::MatrixXd rbdJDot(6, mb.nrDof());
    rbdJac.fullJacobian(mb, rbdShortJ, rbdJ);
    rbdJac.fullJacobian(mb, rbdShortJDot, rbdJDot);

    // Checks
    REQUIRE(rbdJ.isApprox(J));
    REQUIRE(rbdJDot.isApprox(dJ));
    REQUIRE(rbdJ.isApprox(cdm::MotionJacobianOfOrder<0>(bodyName, info, tree)));
    REQUIRE(rbdJDot.isApprox(cdm::MotionJacobianOfOrder<1>(bodyName, info, tree)));

    // Get \aleph
    Eigen::VectorXd aleph = GetAleph(model, mc, t);
    // Get dx
    Eigen::VectorXd dx_j = BigJ * aleph;

    // Check dx
    const auto& factors = coma::factorial_factors<double, order>;
    auto linkMotions = (mc.world.inverse() * mc.linkMotions[m.bodyIndexByName(bodyName)]).motion();
    for (int i = 0; i < mc.world.order(); ++i) {
        REQUIRE(linkMotions[i].vector().isApprox(dx_j.segment<6>(6 * i) * factors[i]));
    }
}