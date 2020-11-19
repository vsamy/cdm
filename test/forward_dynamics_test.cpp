#include "SimpleHumanModel.hpp"
#include "macros.hpp"
#include "model_generation.hpp"
#include <catch2/catch.hpp>
#include <cod/FD.hpp>
#include <rbdyn/FA.h>
#include <rbdyn/FD.h>
#include <rbdyn/FK.h>
#include <rbdyn/FV.h>
#include <tuple>

TEST_CASE("FD", "[FD]")
{
    constexpr int order = 5;

    rbd::MultiBody mb;
    rbd::MultiBodyConfig mbc;
    rbd::MultiBodyGraph mbg;
    std::tie(mb, mbc, mbg) = rbd::makeHumanBody();

    coma::Model<double> model = coma::makeHumanBody<double>();
    coma::ModelConfig<double, 6, order> mcNoGrav;
    coma::ModelConfig<double, 6, order> mc;
    // coma::ModelConfig<double, 6, order> mc2;

    int nt = 21;
    double dt = 1e-8;
    int t1 = nt / 2;
    int t2 = t1 + 1;
    auto data = GenerateData<order>(mb, mbc, nt, dt);
    data.setCurData(t1);

    // No gravity
    data.gravity.setZero();
    Init(data, model, mcNoGrav);
    FK(model, mcNoGrav);
    ID(model, mcNoGrav);

    // With gravity
    data.gravity = Eigen::Vector3d(0, 0, 9.81); // gravity acceleration applied on robot
    Init(data, model, mc);
    FK(model, mc);
    ID(model, mc);

    // RBDyn
    Init(data, mb, mbc);
    rbd::InverseDynamics id{ mb };
    rbd::forwardKinematics(mb, mbc);
    rbd::forwardVelocity(mb, mbc);
    rbd::forwardAcceleration(mb, mbc);
    id.inverseDynamics(mb, mbc);

    std::vector<Eigen::VectorXd> tauP(mb.nrJoints());
    std::vector<Eigen::VectorXd> tauF(mb.nrJoints());
    const auto& factors = coma::factorial_factors<double, ord>;
    for (int i = 0; i < mb.nrJoints(); ++i) {
        auto G = makeDiag<ord>(mbc.motionSubspace[i]);
        tauP[i] = G.transpose() * treeNoGrav.jointMomentums[i].vector();
        tauF[i] = tree.jointTorques[i];
        int dof = mb.joint(i).dof();
        for (int n = 0; n < ord; ++n) {
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
    Eigen::MatrixXd dqs(mb.nrDof(), 5);
    Eigen::VectorXd tau(mb.nrDof());
    Eigen::VectorXd f = Eigen::VectorXd::Zero(mb.nrDof());
    for (int i = 0; i < 5; ++i) {
        dqs.col(i) = info.dqs[i].col(t);
        info.dqs[i].col(t).setZero();
    }

    // Check recursion
    tree.init(t, info);
    ID(info, tree);
    for (int n = 0; n < ord; ++n) {
        for (int j = 0; j < mb.nrJoints(); ++j) {
            int dof = mb.joint(j).dof();
            tau.segment(mb.jointPosInDof(j), dof) = tauF[j].segment(n * dof, dof);
            f.segment(mb.jointPosInDof(j), dof) = tree.jointTorques[j].segment(n * dof, dof);
        }
        info.dqs[n].col(t) = standard_FD(info, tree, tau - f);
        tree.init(t, info);
        ID(info, tree);
    }
    Eigen::MatrixXd yErr2{ mb.nrJoints(), ord };
    for (int i = 0; i < mb.nrJoints(); ++i) {
        int pos = mb.jointPosInDof(i);
        int dof = mb.joint(i).dof();
        for (int n = 0; n < ord; ++n) {
            yErr2.col(n)(i) = (info.dqs[n].col(t).segment(pos, dof) - dqs.col(n).segment(pos, dof)).norm();
        }
    }

    // FD
    std::vector<Eigen::VectorXd> y = FD(info, tree, tauP);
    Eigen::VectorXd ddq{ mb.nrDof() };
    for (int i = 0; i < mb.nrJoints(); ++i) {
        int pos = mb.jointPosInDof(i);
        int dof = mb.joint(i).dof();
        ddq.segment(pos, dof) = y[i].segment(dof, dof);
    }

    // Check with RBDyn
    Eigen::VectorXd alphaD = rbd::dofToVector(mb, mbc.alphaD);
    double rbdErr = (alphaD - ddq).norm();

    // Check motion
    Eigen::MatrixXd yErr{ mb.nrJoints(), ord };
    for (int i = 0; i < mb.nrJoints(); ++i) {
        int dof = mb.joint(i).dof();
        for (int n = 0; n < ord; ++n) {
            yErr.col(n)(i) = (y[i].segment(n * dof, dof) * factors[n] - dqs.col(n).segment(mb.jointPosInDof(i), dof)).norm();
        }
    }

    std::cout << "\nComparison with RBDyn output: " << rbdErr
              << "\nComparison with ID input: ";
    if (VERBOSITY) {
        std::cout << "\nMax dq norm error: " << yErr.col(0).maxCoeff()
                  << "\nMax d2q norm error: " << yErr.col(1).maxCoeff()
                  << "\nMax d3q norm error: " << yErr.col(2).maxCoeff()
                  << "\nMax d4q norm error: " << yErr.col(3).maxCoeff()
                  << "\nMax d5q norm error: " << yErr.col(4).maxCoeff();
    } else {
        std::cout << "Max norm error = " << yErr.maxCoeff();
    }
    std::cout << "\nComparison for recursive FD:";
    if (VERBOSITY) {
        std::cout << "\nMax dq norm error: " << yErr2.col(0).maxCoeff()
                  << "\nMax d2q norm error: " << yErr2.col(1).maxCoeff()
                  << "\nMax d3q norm error: " << yErr2.col(2).maxCoeff()
                  << "\nMax d4q norm error: " << yErr2.col(3).maxCoeff()
                  << "\nMax d5q norm error: " << yErr2.col(4).maxCoeff();
    } else {
        std::cout << "Max norm error = " << yErr2.maxCoeff();
    }
}