#include "SimpleHumanModel.hpp"
#include "doctest/doctest.h"
#include "macros.hpp"
#include "model_generation.hpp"
#include <RBDyn/FA.h>
#include <RBDyn/FD.h>
#include <RBDyn/FK.h>
#include <RBDyn/FV.h>
#include <RBDyn/ID.h>
#include <RBDyn/Jacobian.h>
#include <cdm/BasicJacobian.hpp>
#include <cdm/FK.hpp>
#include <cdm/ID.hpp>
#include <cdm/JointForceJacobian.hpp>
#include <cdm/JointMomentumJacobian.hpp>
#include <cdm/LinkForceJacobian.hpp>
#include <cdm/LinkMomentumJacobian.hpp>
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
    std::string bodyName{ "RARM1" };

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
    Eigen::MatrixXd J(6, model.nDof());
    Eigen::MatrixXd dJ(6, model.nDof());
    for (auto i = 0; i < model.nLinks(); ++i) {
        auto dof = model.joint(i).dof();
        auto posInDof = model.jointPosInDof(i);
        J.block(0, posInDof, 6, dof) = BigJ.block(0, posInDof * order, 6, dof);
        dJ.block(0, posInDof, 6, dof) = BigJ.block(6, posInDof * order, 6, dof);
    }

    // RBDyn FK
    Init(data, mb, mbc);
    rbd::forwardKinematics(mb, mbc);
    rbd::forwardVelocity(mb, mbc);
    rbd::forwardAcceleration(mb, mbc);

    // Get RBDyn J and dJ
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
    REQUIRE(rbdJ.isApprox(cdm::BasicJacobianOfOrder<0>(model, mc, bodyName)));
    REQUIRE(rbdJDot.isApprox(cdm::BasicJacobianOfOrder<1>(model, mc, bodyName)));

    // Get \aleph
    Eigen::VectorXd aleph = mc.getAleph();
    // Get dx
    Eigen::VectorXd dx_j = BigJ * aleph;

    // Check dx
    const auto& factors = coma::factorial_factors<double, order>;
    size_t bodyIndex = static_cast<size_t>(model.bodyIndexByName(bodyName));
    auto linkMotions = (mc.world.inverse() * mc.bodyMotions[bodyIndex]).motion();
    for (auto i = 0; i < mc.world.order(); ++i) {
        REQUIRE(linkMotions[i].vector().isApprox(dx_j.segment<6>(6 * i) * factors[static_cast<size_t>(i)]));
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
    std::string bodyName{ "RARM1" };

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
    Eigen::MatrixXd BigK = cdm::LinkMomentumJacobian(model, mc, bodyName); // p = BigK * \aleph

    // Get \aleph
    Eigen::VectorXd aleph = mc.getAleph();
    // Get p
    Eigen::VectorXd h_j = BigK * aleph;

    // Check dx
    const auto& factors = coma::factorial_factors<double, order>;
    auto bInd = model.bodyIndexByName(bodyName);
    auto ubInd = static_cast<size_t>(bInd);
    auto gravityEffect = cdm::DiInertia<order>{ model.body(bInd).inertia() } * (mcNoG.bodyMotions[ubInd].inverse() * mc.world.motion());
    auto linkMomentums = mc.bodyMomentums[ubInd] - gravityEffect;
    for (int i = 0; i < order; ++i) {
        REQUIRE((h_j.segment<6>(6 * i) * factors[static_cast<size_t>(i)] - linkMomentums[i].vector()).norm() < 100. * dt);
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
    std::string bodyName{ "RARM1" };

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
    Eigen::MatrixXd BigN = LinkForceJacobian(model, mc, bodyName); // f = BigN * \aleph

    // Get \aleph
    Eigen::VectorXd aleph = mc.getAleph();
    // Get p
    Eigen::VectorXd f_j = BigN * aleph;

    // Check dx
    const auto& factors = coma::factorial_factors<double, order>;
    auto bInd = model.bodyIndexByName(bodyName);
    auto ubInd = static_cast<size_t>(bInd);
    auto gravityMomentum = coma::DiInertiaNd<order>{ model.body(bInd).inertia() } * (mcNoG.bodyMotions[ubInd].inverse() * mc.world.motion());
    for (int n = 0; n < order; ++n) {
        gravityMomentum[n] /= factors[static_cast<size_t>(n)];
    }
    Eigen::VectorXd gravityEffect = cdm::generateD(cdm::CrossN<order>{ mcNoG.bodyMotions[ubInd].motion() }) * gravityMomentum.vector();
    for (int n = 0; n < order; ++n) {
        gravityEffect.segment<6>(6 * n) *= factors[static_cast<size_t>(n)];
    }
    Eigen::VectorXd linkForces = mc.bodyForces[ubInd].vector() - gravityEffect;
    for (int i = 0; i < order; ++i) {
        REQUIRE((f_j.segment<6>(6 * i) * factors[static_cast<size_t>(i)] - linkForces.segment<6>(6 * i)).norm() < 100. * dt);
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
    std::string bodyName{ "RARM1" };

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
    Eigen::VectorXd aleph = mc.getAleph();
    // Get p
    Eigen::VectorXd h_j = BigB * aleph;

    // Check dx
    const auto& factors = coma::factorial_factors<double, order>;
    auto bInd = model.bodyIndexByName(bodyName);
    auto ubInd = static_cast<size_t>(bInd);
    auto M = getSubTreeInertia(model, mcNoG);
    auto grav = mcNoG.bodyMotions[ubInd].inverse() * mc.world.motion();
    for (int n = 0; n < order; ++n) {
        grav[n] /= factors[static_cast<size_t>(n)];
    }
    Eigen::VectorXd gravityEffect = M[ubInd] * grav.vector();
    for (int n = 0; n < order; ++n) {
        gravityEffect.segment<6>(6 * n) *= factors[static_cast<size_t>(n)];
    }
    Eigen::VectorXd jointMomentums = mc.jointMomentums[ubInd].vector() - gravityEffect;
    for (int i = 0; i < order; ++i) {
        REQUIRE((h_j.segment<6>(6 * i) * factors[static_cast<size_t>(i)] - jointMomentums.segment<6>(6 * i)).norm() < 100. * dt);
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
    std::string bodyName{ "RARM1" };

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
    Eigen::VectorXd aleph = mc.getAleph();
    // Get p
    Eigen::VectorXd f_j = BigQ * aleph;

    // Check dx
    const auto& factors = coma::factorial_factors<double, order>;
    cdm::Index bInd = model.bodyIndexByName(bodyName);
    auto ubInd = static_cast<size_t>(bInd);
    auto M = getSubTreeInertia(model, mcNoG);
    auto grav = mcNoG.bodyMotions[ubInd].inverse() * mc.world.motion();
    for (int n = 0; n < order; ++n) {
        grav[n] /= factors[static_cast<size_t>(n)];
    }
    Eigen::VectorXd gravityEffect = cdm::generateD(cdm::CrossN<order>{ mcNoG.bodyMotions[ubInd].motion() }) * M[ubInd] * grav.vector();
    for (int n = 0; n < order; ++n) {
        gravityEffect.segment<6>(6 * n) *= factors[static_cast<size_t>(n)];
    }
    Eigen::VectorXd jointForces = mc.jointForces[ubInd].vector() - gravityEffect;
    for (int i = 0; i < order; ++i) {
        REQUIRE((f_j.segment<6>(6 * i) * factors[static_cast<size_t>(i)] - jointForces.segment<6>(6 * i)).norm() < 100. * dt);
    }
}
