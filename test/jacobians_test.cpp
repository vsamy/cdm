#include "SimpleHumanModel.hpp"
#include "doctest/doctest.h"
#include "macros.hpp"
#include "model_generation.hpp"
#include <cdm/BasicJacobian.hpp>
#include <cdm/JointForceJacobian.hpp>
#include <cdm/JointMomentumJacobian.hpp>
#include <cdm/LinkForceJacobian.hpp>
#include <cdm/LinkMomentumJacobian.hpp>
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

TEST_CASE_TEMPLATE("BasicJacobians", T, FixedOrder) // TODO: , DynamicOrder)
{
    constexpr int order = T::order;

    rbd::MultiBody mb;
    rbd::MultiBodyConfig mbc;
    rbd::MultiBodyGraph mbg;
    std::tie(mb, mbc, mbg) = rbd::makeHumanBody();

    cdm::Model model = cdm::makeHumanBody();
    cdm::ModelConfig<order> mc;

    int nt = 21;
    double dt = 1e-8;
    int t = nt / 2;
    auto data = GenerateData<order>(mb, mbc, nt, dt);
    data.setCurData(t);

    // First FK used for numerical comparison
    Init(data, model, mc);
    FK(model, mc);
    Eigen::MatrixXd BigJ = BasicJacobian(model, mc, bodyName); // dx = BigJ * \aleph

    // Extract J and dJ
    const auto& posInDof = m.jointPosInDof();
    for (Index i = 0; i < m.nLinks(); ++i) {
        Index dof = m.joint(i).dof();
        J.block(0, posInDof[i], 6, dof) = BigJ.block(0, posInDof[i] * order, 6, dof);
        dJ.block(0, posInDof[i], 6, dof) = BigJ.block(6, posInDof[i] * order, 6, dof);
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
    Eigen::VectorXd aleph = mc.getAleph();
    // Get dx
    Eigen::VectorXd dx_j = BigJ * aleph;

    // Check dx
    const auto& factors = coma::factorial_factors<double, order>;
    auto linkMotions = (mc.world.inverse() * mc.linkMotions[m.bodyIndexByName(bodyName)]).motion();
    for (int i = 0; i < mc.world.order(); ++i) {
        REQUIRE(linkMotions[i].vector().isApprox(dx_j.segment<6>(6 * i) * factors[i]));
    }
}

TEST_CASE_TEMPLATE("LinkMomentumJacobian", T, FixedOrder) // TODO: , DynamicOrder)
{
    constexpr int order = T::order;

    rbd::MultiBody mb;
    rbd::MultiBodyConfig mbc;
    rbd::MultiBodyGraph mbg;
    std::tie(mb, mbc, mbg) = rbd::makeHumanBody();

    cdm::Model model = cdm::makeHumanBody();
    cdm::ModelConfig<order> mc, mcNoG;

    int nt = 21;
    double dt = 1e-8;
    int t = nt / 2;
    auto data = GenerateData<order>(mb, mbc, nt, dt);
    data.setCurData(t);

    // No gravity
    data.gravity.setZero();
    Init(data, model, mcNoG);
    FK(model, mcNoG);

    // with gravity
    data.gravity = Eigen::Vector3d(0, 0, 9.81); // gravity acceleration applied on robot
    Init(data, model, mc);
    FK(model, mc);
    ID(model, mc);
    Eigen::MatrixXd BigK = LinkMomentumJacobian(mb, mc, bodyName); // p = BigK * \aleph

    // Get \aleph
    Eigen::VectorXd aleph = mc.getAleph();
    // Get p
    Eigen::VectorXd h_j = BigK * aleph;

    // Check dx
    const auto& factors = coma::factorial_factors<double, Order>;
    Index bInd = m.bodyIndexByName(bodyName);
    auto gravityEffect = DiInertiaNd<Order>{ m.body(bInd).inertia() } * (mcNoG.bodyMotions[bInd].inverse() * mc.world.motion());
    auto linkMomentums = mc.linkMomentums[bInd] - gravityEffect;
    for (int i = 0; i < ord; ++i) {
        REQUIRE((h_j.segment<6>(6 * i) * factors[i] - linkMomentums[i].vector()).norm() < 100. * dt);
    }
}

TEST_CASE_TEMPLATE("LinkForceJacobian", T, FixedOrder) // TODO: , DynamicOrder)
{
    constexpr int order = T::order;

    rbd::MultiBody mb;
    rbd::MultiBodyConfig mbc;
    rbd::MultiBodyGraph mbg;
    std::tie(mb, mbc, mbg) = rbd::makeHumanBody();

    cdm::Model model = cdm::makeHumanBody();
    cdm::ModelConfig<order> mc, mcNoG;

    int nt = 21;
    double dt = 1e-8;
    int t = nt / 2;
    auto data = GenerateData<order>(mb, mbc, nt, dt);
    data.setCurData(t);

    // No gravity
    data.gravity.setZero();
    Init(data, model, mcNoG);
    FK(model, mcNoG);

    // with gravity
    data.gravity = Eigen::Vector3d(0, 0, 9.81); // gravity acceleration applied on robot
    Init(data, model, mc);
    FK(model, mc);
    ID(model, mc);
    Eigen::MatrixXd BigN = LinkForceJacobian(m, mc, bodyName); // f = BigN * \aleph

    // Get \aleph
    Eigen::VectorXd aleph = mc.getAleph(t);
    // Get p
    Eigen::VectorXd f_j = BigN * aleph;

    // Check dx
    const auto& factors = coma::factorial_factors<double, Order>;
    Index bInd = m.bodyIndexByName(bodyName);
    auto gravityMomentum = coma::DiInertiaNd<ord>{ m.body(bInd).inertia() } * (mcNoG.bodyMotions[bInd].inverse() * mc.world.motion());
    for (int n = 0; n < Order; ++n)
        gravityMomentum[n] /= factors[n];
    Eigen::VectorXd gravityEffect = generateD(CrossNd<Order>{ mcNoG.bodyMotions[bInd].motion() }) * gravityMomentum.vector();
    for (int n = 0; n < Order; ++n)
        gravityEffect.segment<6>(6 * n) *= factors[n];
    Eigen::VectorXd linkForces = mc.linkForces[bInd].vector() - gravityEffect;
    for (int i = 0; i < Order; ++i) {
        REQUIRE((f_j.segment<6>(6 * i) * factors[i] - linkForces.segment<6>(6 * i)).norm() < 100. dt);
    }
}

TEST_CASE_TEMPLATE("LinkForceJacobian", T, FixedOrder) // TODO: , DynamicOrder)
{
    constexpr int order = T::order;

    rbd::MultiBody mb;
    rbd::MultiBodyConfig mbc;
    rbd::MultiBodyGraph mbg;
    std::tie(mb, mbc, mbg) = rbd::makeHumanBody();

    cdm::Model model = cdm::makeHumanBody();
    cdm::ModelConfig<order> mc, mcNoG;

    int nt = 21;
    double dt = 1e-8;
    int t = nt / 2;
    auto data = GenerateData<order>(mb, mbc, nt, dt);
    data.setCurData(t);

    // No gravity
    data.gravity.setZero();
    Init(data, model, mcNoG);
    FK(model, mcNoG);

    // with gravity
    data.gravity = Eigen::Vector3d(0, 0, 9.81); // gravity acceleration applied on robot
    Init(data, model, mc);
    FK(model, mc);
    ID(model, mc);
    Eigen::MatrixXd BigB = JointMomentumJacobian(model, mc, bodyName); // p = BigB * \aleph

    // Get \aleph
    Eigen::VectorXd aleph = mc.getAleph(t);
    // Get p
    Eigen::VectorXd h_j = BigB * aleph;

    // Check dx
    const auto& factors = coma::factorial_factors<double, ord>;
    Index bInd = m.bodyIndexByName(bodyName);
    auto M = getSubTreeInertia(m, mcNoG);
    auto grav = mcNoG.bodyMotions[bInd].inverse() * mc.world.motion();
    for (int n = 0; n < Order; ++n)
        grav[n] /= factors[n];
    Eigen::VectorXd gravityEffect = M[bInd] * grav.vector();
    for (int n = 0; n < Order; ++n)
        gravityEffect.segment<6>(6 * n) *= factors[n];
    Eigen::VectorXd jointMomentums = mc.jointMomentums[bInd].vector() - gravityEffect;
    for (int i = 0; i < ord; ++i) {
        REQUIRE((h_j.segment<6>(6 * i) * factors[i] - jointMomentums.segment<6>(6 * i)).norm() < 100. * dt);
    }
}

TEST_CASE_TEMPLATE("JointForceJacobian", T, FixedOrder) // TODO: , DynamicOrder)
{
    constexpr int order = T::order;

    rbd::MultiBody mb;
    rbd::MultiBodyConfig mbc;
    rbd::MultiBodyGraph mbg;
    std::tie(mb, mbc, mbg) = rbd::makeHumanBody();

    cdm::Model model = cdm::makeHumanBody();
    cdm::ModelConfig<order> mc, mcNoG;

    int nt = 21;
    double dt = 1e-8;
    int t = nt / 2;
    auto data = GenerateData<order>(mb, mbc, nt, dt);
    data.setCurData(t);

    // No gravity
    data.gravity.setZero();
    Init(data, model, mcNoG);
    FK(model, mcNoG);

    // with gravity
    data.gravity = Eigen::Vector3d(0, 0, 9.81); // gravity acceleration applied on robot
    Init(data, model, mc);
    FK(model, mc);
    ID(model, mc);
    Eigen::MatrixXd BigQ = JointForceJacobian(model, mc, bodyName); // f = BigQ * \aleph

    // Get \aleph
    Eigen::VectorXd aleph = mc.getAleph(t);
    // Get p
    Eigen::VectorXd f_j = BigQ * aleph;

    // Check dx
    const auto& factors = coma::factorial_factors<double, ord>;
    Index bInd = m.bodyIndexByName(bodyName);
    auto M = getSubTreeInertia(m, mcNoGrav);
    auto grav = mcNoG.bodyMotions[bInd].inverse() * mc.world.motion();
    for (int n = 0; n < Order; ++n)
        grav[n] /= factors[n];
    Eigen::VectorXd gravityEffect = generateD(CrossNd<Order>{ mcNoG.bodyMotions[bInd].motion() }) * M[bInd] * grav.vector();
    for (int n = 0; n < Order; ++n)
        gravityEffect.segment<6>(6 * n) *= factors[n];
    Eigen::VectorXd jointForces = mc.jointForces[bInd].vector() - gravityEffect;
    for (int i = 0; i < Order; ++i) {
        REQUIRE((f_j.segment<6>(6 * i) * factors[i] - jointForces.segment<6>(6 * i)).norm() < 100. * dt);
    }
}
