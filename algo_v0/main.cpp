#include "FD.hpp"
#include "FK.hpp"
#include "FKAD.hpp"
#include "ID.hpp"
#include "IDAD.hpp"
#include "EigenAD.hpp"
#include "JointForceJacobian.hpp"
#include "JointMomentumJacobian.hpp"
#include "LinkForceJacobian.hpp"
#include "LinkMomentumJacobian.hpp"
#include "ModelInfo.hpp"
#include "MotionJacobian.hpp"
#include "Tree.hpp"
#include "utilityAD.hpp"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <RBDyn/FA.h>
#include <RBDyn/FD.h>
#include <RBDyn/FK.h>
#include <RBDyn/FV.h>
#include <RBDyn/ID.h>
#include <RBDyn/Jacobian.h>
#include <benchmark/benchmark.h>

static bool VERBOSITY = false;

template <typename MI, typename Tree>
void checkFK(MI& info, Tree& tree)
{
    constexpr int ord = Tree::order;
    const auto& mb = info.model.mb;
    auto& mbc = info.model.mbc;
    int t1 = info.nt / 2;
    int t2 = t1 + 1;
    mbc.gravity = Eigen::Vector3d(0, 0, 9.81);

    // First FK used for numerical comparison
    tree.init(t1, info);
    FK(info, tree);

    // RBDyn FK
    rbd::forwardKinematics(mb, mbc);
    rbd::forwardVelocity(mb, mbc);
    rbd::forwardAcceleration(mb, mbc);

    // Second FK
    auto tree2 = tree;
    tree2.init(t2, info);
    FK(info, tree2);

    // Check FK with RBDyn
    auto linkWorldInv = tree.linkWorld.inverse();
    Eigen::MatrixXd rbdErr { mb.nrBodies(), 4 };
    for (int i = 0; i < mb.nrJoints(); ++i) {
        auto linkMotions = (linkWorldInv * tree.links[i]).motion();
        rbdErr(i, 0) = (mbc.bodyPosW[i].translation() - tree.links[i].transform().translation()).norm();
        rbdErr(i, 1) = (mbc.bodyPosW[i].rotation().transpose() - tree.links[i].transform().rotation()).norm();
        rbdErr(i, 2) = (mbc.bodyVelB[i].vector() - linkMotions[0].vector()).norm();
        rbdErr(i, 3) = (mbc.bodyAccB[i].vector() - linkMotions[1].vector()).norm();
    }

    // Numerical check
    Eigen::MatrixXd numErr { mb.nrBodies(), ord - 1 };
    for (int i = 0; i < mb.nrBodies(); ++i) {
        auto m1 = (linkWorldInv * tree.links[i]).motion();
        auto m2 = (linkWorldInv * tree2.links[i]).motion();
        auto dV = (m2 - m1) / info.dt;
        for (int n = 0; n < ord - 1; ++n) {
            numErr(i, n) = (dV[n] - m2[n + 1]).vector().norm();
        }
    }

    std::cout << "\nComparison with RBDyn output: ";
    if (VERBOSITY) {
        std::cout << "\nMax translation norm error: " << rbdErr.col(0).maxCoeff()
                  << "\nMax rotation norm error: " << rbdErr.col(1).maxCoeff()
                  << "\nMax velocity norm error: " << rbdErr.col(2).maxCoeff()
                  << "\nMax acceleration norm error: " << rbdErr.col(3).maxCoeff();
    } else {
        std::cout << "Max norm error = " << rbdErr.maxCoeff();
    }
    std::cout << "\nComparison with numerical diff: ";
    if (VERBOSITY) {
        std::cout << "\nMax acceleration norm error: " << numErr.col(0).maxCoeff()
                  << "\nMax jerk norm error: " << numErr.col(1).maxCoeff()
                  << "\nMax snap norm error: " << numErr.col(2).maxCoeff()
                  << "\nMax crackle norm error: " << numErr.col(3).maxCoeff();
    } else {
        std::cout << "Max norm error = " << numErr.maxCoeff();
    }
}

template <typename MI, typename Tree>
void checkID(MI& info, Tree& tree, bool withGravity)
{
    constexpr int ord = Tree::order;
    const auto& mb = info.model.mb;
    auto& mbc = info.model.mbc;
    int t1 = info.nt / 2;
    int t2 = t1 + 1;
    mbc.gravity = withGravity ? Eigen::Vector3d(0, 0, 9.81) : Eigen::Vector3d::Zero(); // gravity acceleration applied on robot

    // First ID used for numerical comparison
    tree.init(t1, info);
    ID(info, tree);

    // RBDyn ID
    rbd::InverseDynamics id { mb };
    rbd::forwardKinematics(mb, mbc);
    rbd::forwardVelocity(mb, mbc);
    rbd::forwardAcceleration(mb, mbc);
    id.inverseDynamics(mb, mbc);
    const auto& rbdFJoint = id.f();
    std::vector<Eigen::VectorXd> rbdTau(mb.nrBodies());
    std::vector<sva::ForceVecd> rbdLinkP(mb.nrBodies());
    std::vector<sva::ForceVecd> rbdLinkFD(mb.nrBodies());
    std::vector<sva::ForceVecd> rbdLinkF(mb.nrBodies());
    std::vector<sva::ForceVecd> rbdLinkDP(mb.nrBodies());
    for (int i = 0; i < mb.nrBodies(); ++i) {
        rbdLinkP[i] = mb.body(i).inertia() * mbc.bodyVelB[i];
        rbdLinkDP[i] = mb.body(i).inertia() * mbc.bodyAccB[i];
        rbdLinkFD[i] = rbdLinkP[i];
        rbdLinkF[i] = rbdLinkDP[i] + mbc.bodyVelB[i].crossDual(rbdLinkP[i]);
        rbdTau[i] = Eigen::Map<Eigen::VectorXd>(mbc.jointTorque[i].data(), mbc.jointTorque[i].size());
    }

    // Second ID
    auto tree2 = tree;
    tree2.init(t2, info);
    ID(info, tree2);

    // Check with RBDyn
    Eigen::MatrixXd rbdErr { mb.nrBodies(), 6 };
    for (int i = 0; i < mb.nrBodies(); ++i) {
        int dof = mb.joint(i).dof();
        rbdErr(i, 0) = (rbdLinkP[i].vector() - tree.linkMomentums[i][0].vector()).norm();
        rbdErr(i, 1) = (rbdLinkDP[i].vector() - tree.linkMomentums[i][1].vector()).norm();
        rbdErr(i, 2) = (rbdLinkFD[i].vector() - tree.linkForces[i][0].vector()).norm();
        rbdErr(i, 3) = (rbdLinkF[i].vector() - tree.linkForces[i][1].vector()).norm();
        rbdErr(i, 4) = (rbdFJoint[i].vector() - tree.jointForces[i][1].vector()).norm();
        rbdErr(i, 5) = (rbdTau[i] - tree.jointTorques[i].segment(dof, dof)).norm();
    }

    // Numerical check
    Eigen::MatrixXd numLinkMomentumErr { mb.nrBodies(), ord - 1 };
    Eigen::MatrixXd numLinkForceErr { mb.nrBodies(), ord - 1 };
    Eigen::MatrixXd numJointMomentumErr { mb.nrBodies(), ord - 1 };
    Eigen::MatrixXd numJointForceErr { mb.nrBodies(), ord - 1 };
    for (int i = 0; i < mb.nrBodies(); ++i) {
        auto dlP = (tree2.linkMomentums[i] - tree.linkMomentums[i]) / info.dt;
        auto dlF = (tree2.linkForces[i] - tree.linkForces[i]) / info.dt;
        auto djP = (tree2.jointMomentums[i] - tree.jointMomentums[i]) / info.dt;
        auto djF = (tree2.jointForces[i] - tree.jointForces[i]) / info.dt;
        for (int n = 0; n < ord - 1; ++n) {
            numLinkMomentumErr(i, n) = (dlP[n] - tree2.linkMomentums[i][n + 1]).vector().norm();
            numLinkForceErr(i, n) = (dlF[n] - tree2.linkForces[i][n + 1]).vector().norm();
            numJointMomentumErr(i, n) = (djP[n] - tree2.jointMomentums[i][n + 1]).vector().norm();
            numJointForceErr(i, n) = (djF[n] - tree2.jointForces[i][n + 1]).vector().norm();
        }
    }

    // AD check
    using CppAD::AD;

    // domain space vector
    std::vector<AD<AD<double>>> X(2 * mb.nrDof());
    std::vector<AD<AD<double>>> Y(mb.nrDof());
    int pos = 0;
    for (int i = 0; i < mb.nrJoints(); ++i) {
        int dof = mb.joint(i).dof();
        for (int j = 0; j < dof; ++j) {
            X[pos + j] = info.dqs[0](pos + j, t2);
            X[pos + j + dof] = info.dqs[1](pos + j, t2);
        }
        pos += dof;
    }

    // // declare independent variables and start recording operation sequence
    // CppAD::Independent(X);

    // // range space vector
    // Y = IDAD(X, info.model); // value during recording of operations

    // // store operation sequence in f: X -> Y and stop recording
    // CppAD::ADFun<AD<double>> f(X, Y);

    // // compute derivative using operation sequence stored in f
    // std::vector<AD<double>> jac(mb.nrDof() * 2 * mb.nrDof()); // Jacobian of f (m by n matrix)
    // std::vector<AD<double>> X2(2 * mb.nrDof());
    // pos = 0;
    // for (int i = 0; i < mb.nrJoints(); ++i) {
    //     int dof = mb.joint(i).dof();
    //     for (int j = 0; j < dof; ++j) {
    //         X2[pos + j] = info.dqs[0](pos + j, t2);
    //         X2[pos + j + dof] = info.dqs[1](pos + j, t2);
    //     }
    //     pos += dof;
    // }

    // // declare independent variables and start recording operation sequence
    // CppAD::Independent(X2);

    // jac = f.Jacobian(X2); // Jacobian for operation sequence
    // auto jacT = TransposeAD(jac, 6);

    // // store operation sequence in f: X -> Y and stop recording
    // CppAD::ADFun<double> f2(X2, jacT);

    // // compute 2nd-order derivative using operation sequence stored in f
    // std::vector<double> H(mb.nrDof() * 2 * mb.nrDof() * 2 * mb.nrDof());
    // std::vector<double> x(2 * mb.nrDof()); // domain space vector
    // std::vector<double> dx(2 * mb.nrDof()); // domain space vector
    // std::vector<double> ddx(2 * mb.nrDof()); // domain space vector
    // pos = 0;
    // for (int i = 0; i < mb.nrJoints(); ++i) {
    //     int dof = mb.joint(i).dof();
    //     for (int j = 0; j < dof; ++j) {
    //         x[pos + j] = info.dqs[0](pos + j, t2);
    //         x[pos + j + dof] = info.dqs[1](pos + j, t2);
    //         dx[pos + j] = info.dqs[1](pos + j, t2);
    //         dx[pos + j + dof] = info.dqs[2](pos + j, t2);
    //         ddx[pos + j] = info.dqs[2](pos + j, t2);
    //         ddx[pos + j + dof] = info.dqs[3](pos + j, t2);
    //     }
    //     pos += dof;
    // }
    // H = f2.Jacobian(x);

    // auto jacTNoAD = CastOutAD(jacT);
    // auto HT = TransposeAD(H, 2 * mb.nrDof() * 2 * mb.nrDof());
    // auto Hdq = ProductAD(HT, dx, 2 * mb.nrDof());
    // auto Hdqdq = ProductAD(Hdq, dx, 2 * mb.nrDof());
    // auto Jddq = ProductAD(jacTNoAD, ddx, 2 * mb.nrDof());

    // // DX = J * dQ
    // auto DX = ProductAD(jacTNoAD, dx, 2 * mb.nrDof());
    // // D2X = J * d2Q + dJ * DQ
    // auto D2X = SumAD(Jddq, Hdqdq);

    std::cout << "\nComparison with RBDyn output: ";
    if (VERBOSITY) {
        std::cout << "\nMax link momentum norm error: " << rbdErr.col(0).maxCoeff()
                  << "\nMax link momentum derivative norm error: " << rbdErr.col(1).maxCoeff()
                  << "\nMax link force (integral) norm error: " << rbdErr.col(2).maxCoeff()
                  << "\nMax link force norm error: " << rbdErr.col(3).maxCoeff()
                  << "\nMax joint force norm error: " << rbdErr.col(4).maxCoeff()
                  << "\nMax joint torque norm error: " << rbdErr.col(5).maxCoeff();
    } else {
        std::cout << "Max norm error = " << rbdErr.maxCoeff();
    }

    // Please remember that here the 0-order force is the momentum and f[1] corresponds to the classical force so f[0] == p[0]
    // And we have df[0]/dt == dp[0]/dt == p[1], thus df[0]/dt != f[1] (the numerical derivative of f[0] is p[1] and not f[1])
    std::cout << "\nComparison with numerical diff: ";
    if (VERBOSITY) {
        if (!withGravity) {
            for (int n = 0; n < 4; ++n) {
                std::cout << "\nMax link momentum " << n + 1 << "-derivative norm error: " << numLinkMomentumErr.col(n).maxCoeff()
                          << "\nMax joint momentum " << n + 1 << "-derivative norm error: " << numJointMomentumErr.col(n).maxCoeff();
            }
        }
        for (int n = 1; n < 4; ++n) {
            std::cout << "\nMax link force " << n << "-derivative norm error: " << numLinkForceErr.col(n).maxCoeff();
            std::cout << "\nMax joint force " << n << "-derivative norm error: " << numJointForceErr.col(n).maxCoeff();
        }
    } else {
        if (withGravity) {
            std::cout << "Max norm error = " << std::max(numLinkForceErr.rightCols<ord - 2>().maxCoeff(), numJointForceErr.rightCols<ord - 2>().maxCoeff());
        } else {
            std::cout << "Max norm error = " << std::max({ numLinkMomentumErr.maxCoeff(), numLinkForceErr.rightCols<ord - 2>().maxCoeff(), numJointMomentumErr.maxCoeff(), numJointForceErr.rightCols<ord - 2>().maxCoeff() });
        }
    }
}

template <typename MI, typename Tree>
void checkMotionJacobian(const std::string& bodyName, MI& info, Tree& tree)
{
    const auto& mb = info.model.mb;
    auto& mbc = info.model.mbc;
    int t = info.nt / 2;
    mbc.gravity = Eigen::Vector3d(0, 0, 9.81); // gravity acceleration applied on robot

    tree.init(t, info);
    FK(info, tree);
    Eigen::MatrixXd BigJ = MotionJacobian(bodyName, info, tree); // dx = BigJ * \aleph
    Eigen::MatrixXd J(6, mb.nrDof()); // v = J * \psi
    Eigen::MatrixXd dJ(6, mb.nrDof()); // dv = dJ * \psi + J * d\psi

    // Extract J and dJ
    const auto& pos = mb.jointsPosInDof();
    for (int i = 0; i < mb.nrJoints(); ++i) {
        int dof = mb.joint(i).dof();
        J.block(0, pos[i], 6, dof) = BigJ.block(0, pos[i] * Tree::order, 6, dof);
        dJ.block(0, pos[i], 6, dof) = BigJ.block(6, pos[i] * Tree::order, 6, dof);
    }

    // Get RBDyn J and dJ
    rbd::forwardKinematics(mb, mbc);
    rbd::forwardVelocity(mb, mbc);
    rbd::forwardAcceleration(mb, mbc);
    auto rbdJac = rbd::Jacobian(mb, bodyName);
    auto rbdShortJ = rbdJac.bodyJacobian(mb, mbc);
    auto rbdShortJDot = rbdJac.bodyJacobianDot(mb, mbc);
    Eigen::MatrixXd rbdJ(6, mb.nrDof());
    Eigen::MatrixXd rbdJDot(6, mb.nrDof());
    rbdJac.fullJacobian(mb, rbdShortJ, rbdJ);
    rbdJac.fullJacobian(mb, rbdShortJDot, rbdJDot);

    // Checks
    double JErr = (rbdJ - J).norm();
    double dJErr = (rbdJDot - dJ).norm();
    double JpErr = (rbdJ - MotionJacobianOfOrder<0>(bodyName, info, tree)).norm();
    double dJpErr = (rbdJDot - MotionJacobianOfOrder<1>(bodyName, info, tree)).norm();
    std::cout << "\nComparison with RBDyn output: ";
    if (VERBOSITY) {
        std::cout << "\nJacobian (CMTM full Jac) Froebius norm error: " << JErr
                  << "\nJacobian dot (CMTM full Jac) Froebius norm error: " << dJErr
                  << "\nJacobian (CMTM direct computation Jac) Froebius norm error: " << JpErr
                  << "\nJacobian dot (CMTM direct computation Jac) Froebius norm error: " << dJpErr;
    } else {
        std::cout << "Max norm error = " << std::max({ JErr, dJErr, JpErr, dJpErr });
    }

    // Get \aleph
    Eigen::VectorXd aleph = info.getAleph(t);
    // Get dx
    Eigen::VectorXd dx_j = BigJ * aleph;

    // Check dx
    std::vector<double> dxErr;
    constexpr int dx_size = 6 * Tree::order;
    const auto& factors = coma::factorial_factors<double, Tree::order>;
    auto linkMotions = (tree.linkWorld.inverse() * tree.links[mb.bodyIndexByName(bodyName)]).motion();
    for (int i = 0; i < Tree::order; ++i) {
        dxErr.push_back((dx_j.segment<6>(6 * i) * factors[i] - linkMotions[i].vector()).norm());
    }

    std::cout << "\nComparison of dx = J * \\aleph with FK output: ";
    if (VERBOSITY) {
        std::cout << "\nMax velocity norm error: " << dxErr[0]
                  << "\nMax acceleration norm error: " << dxErr[1]
                  << "\nMax jerk norm error: " << dxErr[2]
                  << "\nMax snap norm error: " << dxErr[3]
                  << "\nMax crackle norm error: " << dxErr[4];
    } else {
        std::cout << "Max norm error = " << *std::max_element(dxErr.begin(), dxErr.end());
    }
}

template <typename MI, typename Tree>
void checkLinkMomentumJacobian(const std::string& bodyName, MI& info, Tree& tree)
{
    constexpr int ord = Tree::order;
    const auto& mb = info.model.mb;
    auto& mbc = info.model.mbc;
    int t = info.nt / 2;

    // No gravity
    mbc.gravity.setZero();
    auto treeNoGrav = tree;
    treeNoGrav.init(t, info);
    FK(info, treeNoGrav);

    // with gravity
    mbc.gravity = Eigen::Vector3d(0, 0, 9.81); // gravity acceleration applied on robot
    tree.init(t, info);
    ID(info, tree);
    Eigen::MatrixXd BigK = LinkMomentumJacobian(bodyName, info, tree); // p = BigK * \aleph

    // Get \aleph
    Eigen::VectorXd aleph = info.getAleph(t);
    // Get p
    Eigen::VectorXd p_j = BigK * aleph;

    // Check dx
    std::vector<double> dpErr;
    const auto& factors = coma::factorial_factors<double, ord>;
    int bInd = mb.bodyIndexByName(bodyName);
    auto gravityEffect = coma::DiInertiaNd<ord>{ tree.linkInertias[bInd] } * (treeNoGrav.links[bInd].inverse() * tree.linkWorld.motion());
    auto linkMomentums = tree.linkMomentums[bInd] - gravityEffect;
    for (int i = 0; i < ord; ++i) {
        dpErr.push_back((p_j.segment<6>(6 * i) * factors[i] - linkMomentums[i].vector()).norm());
    }

    std::cout << "\nComparison of p = K * \\aleph with ID output: ";
    if (VERBOSITY) {
        for (int n = 0; n < ord; ++n) {
            std::cout << "\nMax " << n << "-derivative link momentum norm error: " << dpErr[n];
        }
    } else {
        std::cout << "Max norm error = " << *std::max_element(dpErr.begin(), dpErr.end());
    }
}

template <typename MI, typename Tree>
void checkLinkForceJacobian(const std::string& bodyName, MI& info, Tree& tree)
{
    constexpr int ord = Tree::order;
    const auto& mb = info.model.mb;
    auto& mbc = info.model.mbc;
    int t = info.nt / 2;

    // No gravity
    mbc.gravity.setZero();
    auto treeNoGrav = tree;
    treeNoGrav.init(t, info);
    FK(info, treeNoGrav);

    // with gravity
    mbc.gravity = Eigen::Vector3d(0, 0, 9.81); // gravity acceleration applied on robot
    tree.init(t, info);
    ID(info, tree);
    Eigen::MatrixXd BigN = LinkForceJacobian(bodyName, info, tree); // f = BigN * \aleph

    // Get \aleph
    Eigen::VectorXd aleph = info.getAleph(t);
    // Get p
    Eigen::VectorXd f_j = BigN * aleph;

    // Check dx
    std::vector<double> fErr;
    const auto& factors = coma::factorial_factors<double, ord>;
    int bInd = mb.bodyIndexByName(bodyName);
    auto gravityMomentum = coma::DiInertiaNd<ord>{ tree.linkInertias[bInd] } * (treeNoGrav.links[bInd].inverse() * tree.linkWorld.motion());
    for (int n = 0; n < ord; ++n)
        gravityMomentum[n] /= factors[n];
    Eigen::VectorXd gravityEffect = generateD(coma::CrossNd<ord> { treeNoGrav.links[bInd].motion() }) * gravityMomentum.vector();
    for (int n = 0; n < ord; ++n)
        gravityEffect.segment<6>(6 * n) *= factors[n];
    Eigen::VectorXd linkForces = tree.linkForces[bInd].vector() - gravityEffect;
    for (int i = 0; i < ord; ++i) {
        fErr.push_back((f_j.segment<6>(6 * i) * factors[i] - linkForces.segment<6>(6 * i)).norm());
    }

    std::cout << "\nComparison of f = N * \\aleph with ID output: ";
    if (VERBOSITY) {
        for (int n = 0; n < ord; ++n) {
            std::cout << "\nMax " << n << "-derivative link force norm error: " << fErr[n];
        }
    } else {
        std::cout << "Max norm error = " << *std::max_element(fErr.begin(), fErr.end());
    }
}

template <typename MI, typename Tree>
void checkJointMomentumJacobian(const std::string& bodyName, MI& info, Tree& tree)
{
    constexpr int ord = Tree::order;
    const auto& mb = info.model.mb;
    auto& mbc = info.model.mbc;
    int t = info.nt / 2;

    // No gravity
    mbc.gravity.setZero();
    auto treeNoGrav = tree;
    treeNoGrav.init(t, info);
    FK(info, treeNoGrav);

    // with gravity
    mbc.gravity = Eigen::Vector3d(0, 0, 9.81); // gravity acceleration applied on robot
    tree.init(t, info);
    ID(info, tree);
    Eigen::MatrixXd BigB = JointMomentumJacobian(bodyName, info, tree); // p = BigB * \aleph

    // Get \aleph
    Eigen::VectorXd aleph = info.getAleph(t);
    // Get p
    Eigen::VectorXd p_j = BigB * aleph;

    // Check dx
    std::vector<double> pErr;
    const auto& factors = coma::factorial_factors<double, ord>;
    int bInd = mb.bodyIndexByName(bodyName);
    auto M = getSubTreeInertia(mb, treeNoGrav);
    auto grav = treeNoGrav.links[bInd].inverse() * tree.linkWorld.motion();
    for (int n = 0; n < ord; ++n)
        grav[n] /= factors[n];
    Eigen::VectorXd gravityEffect = M[bInd] * grav.vector();
    for (int n = 0; n < ord; ++n)
        gravityEffect.segment<6>(6 * n) *= factors[n];
    Eigen::VectorXd jointMomentums = tree.jointMomentums[bInd].vector() - gravityEffect;
    for (int i = 0; i < ord; ++i) {
        pErr.push_back((p_j.segment<6>(6 * i) * factors[i] - jointMomentums.segment<6>(6 * i)).norm());
    }

    std::cout << "\nComparison of p = B * \\aleph with ID output: ";
    if (VERBOSITY) {
        for (int n = 0; n < ord; ++n) {
            std::cout << "\nMax " << n << "-derivative joint momentum norm error: " << pErr[n];
        }
    } else {
        std::cout << "Max norm error = " << *std::max_element(pErr.begin(), pErr.end());
    }
}

template <typename MI, typename Tree>
void checkJointForceJacobian(const std::string& bodyName, MI& info, Tree& tree)
{
    constexpr int ord = Tree::order;
    const auto& mb = info.model.mb;
    auto& mbc = info.model.mbc;
    int t = info.nt / 2;

    // No gravity
    mbc.gravity.setZero();
    auto treeNoGrav = tree;
    treeNoGrav.init(t, info);
    FK(info, treeNoGrav);

    // with gravity
    mbc.gravity = Eigen::Vector3d(0, 0, 9.81); // gravity acceleration applied on robot
    tree.init(t, info);
    ID(info, tree);
    Eigen::MatrixXd BigQ = JointForceJacobian(bodyName, info, tree); // f = BigQ * \aleph

    // Get \aleph
    Eigen::VectorXd aleph = info.getAleph(t);
    // Get p
    Eigen::VectorXd f_j = BigQ * aleph;

    // Check dx
    std::vector<double> fErr;
    const auto& factors = coma::factorial_factors<double, ord>;
    int bInd = mb.bodyIndexByName(bodyName);
    auto M = getSubTreeInertia(mb, treeNoGrav);
    auto grav = treeNoGrav.links[bInd].inverse() * tree.linkWorld.motion();
    for (int n = 0; n < ord; ++n)
        grav[n] /= factors[n];
    Eigen::VectorXd gravityEffect = generateD(coma::CrossNd<ord> { treeNoGrav.links[bInd].motion() }) * M[bInd] * grav.vector();
    for (int n = 0; n < ord; ++n)
        gravityEffect.segment<6>(6 * n) *= factors[n];
    Eigen::VectorXd jointForces = tree.jointForces[bInd].vector() - gravityEffect;
    for (int i = 0; i < ord; ++i) {
        fErr.push_back((f_j.segment<6>(6 * i) * factors[i] - jointForces.segment<6>(6 * i)).norm());
    }

    std::cout << "\nComparison of f = Q * \\aleph with ID output: ";
    if (VERBOSITY) {
        for (int n = 0; n < ord; ++n) {
            std::cout << "\nMax " << n << "-derivative joint force norm error: " << fErr[n];
        }
    } else {
        std::cout << "Max norm error = " << *std::max_element(fErr.begin(), fErr.end());
    }
}

template <typename MI, typename Tree>
void checkFD(MI& info, Tree& tree)
{
    constexpr int ord = Tree::order;
    const auto& mb = info.model.mb;
    auto& mbc = info.model.mbc;
    int t = info.nt / 2;

    // No gravity
    mbc.gravity.setZero();
    auto treeNoGrav = tree;
    treeNoGrav.init(t, info);
    ID(info, treeNoGrav);

    // With gravity
    mbc.gravity = Eigen::Vector3d(0, 0, 9.81); // gravity acceleration applied on robot
    tree.init(t, info);
    ID(info, tree);

    rbd::InverseDynamics id { mb };
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
    Eigen::VectorXd tauF1 { mb.nrDof() };
    for (int i = 0; i < mb.nrJoints(); ++i) {
        int dof = mb.joint(i).dof();
        tauF1.segment(mb.jointPosInDof(i), dof) = tauF[i].segment(dof, dof);
    }
    mbc.jointTorque = rbd::vectorToDof(mb, tauF1);
    rbd::ForwardDynamics fd { mb };
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
    Eigen::MatrixXd yErr2 { mb.nrJoints(), ord };
    for (int i = 0; i < mb.nrJoints(); ++i) {
        int pos = mb.jointPosInDof(i);
        int dof = mb.joint(i).dof();
        for (int n = 0; n < ord; ++n) {
            yErr2.col(n)(i) = (info.dqs[n].col(t).segment(pos, dof) - dqs.col(n).segment(pos, dof)).norm();
        }
    }

    // FD
    std::vector<Eigen::VectorXd> y = FD(info, tree, tauP);
    Eigen::VectorXd ddq { mb.nrDof() };
    for (int i = 0; i < mb.nrJoints(); ++i) {
        int pos = mb.jointPosInDof(i);
        int dof = mb.joint(i).dof();
        ddq.segment(pos, dof) = y[i].segment(dof, dof);
    }

    // Check with RBDyn
    Eigen::VectorXd alphaD = rbd::dofToVector(mb, mbc.alphaD);
    double rbdErr = (alphaD - ddq).norm();

    // Check motion
    Eigen::MatrixXd yErr { mb.nrJoints(), ord };
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

template <typename MI, typename Tree>
void checkNum2Diff(MI& info, Tree& tree)
{
    // On FK
    constexpr int ord = Tree::order;
    const auto& mb = info.model.mb;
    auto& mbc = info.model.mbc;
    int t = info.nt / 2;

    auto getErr = [&mb, dt = info.dt](const Tree& treeCur, const Tree& treeOld, Eigen::MatrixXd& numDiff) {
        constexpr int ord = Tree::order;
        auto linkWorldInv = treeCur.linkWorld.inverse();
        for (int i = 0; i < mb.nrBodies(); ++i) {
            auto mOld = (linkWorldInv * treeOld.links[i]).motion();
            auto mCur = (linkWorldInv * treeCur.links[i]).motion();
            auto dV = (mCur - mOld) / dt;
            for (int n = 0; n < ord - 1; ++n) {
                numDiff.col(n).segment<6>(6 * i) = dV[n].vector();
            }
        }
    };

    // FK at t
    tree.init(t, info);
    FK(info, tree);

    // FK at t + 1
    auto tree2 = tree;
    tree2.init(t + 1, info);
    FK(info, tree2);

    // Get numerical diff
    Eigen::MatrixXd numDiff1 { 6 * mb.nrBodies(), ord - 1 };
    getErr(tree2, tree, numDiff1);

    // FK at t + 2
    tree.init(t + 2, info);
    FK(info, tree);

    // Get numerical diff
    Eigen::MatrixXd numDiff2 { 6 * mb.nrBodies(), ord - 1 };
    getErr(tree, tree2, numDiff2);

    // Get 2nd order numerical diff
    Eigen::MatrixXd num2Diff = (numDiff2 - numDiff1) / info.dt;

    // Check error
    Eigen::MatrixXd numErr = Eigen::MatrixXd::Zero(mb.nrBodies(), ord - 2);
    auto linkWorldInv = tree.linkWorld.inverse();
    for (int i = 0; i < mb.nrBodies(); ++i) {
        auto mOld = (linkWorldInv * tree2.links[i]).motion();
        auto mCur = (linkWorldInv * tree.links[i]).motion();
        for (int n = 0; n < ord - 2; ++n) {
            numErr(i, n) = (num2Diff.col(n).segment<6>(6 * i) - mCur[n + 2].vector()).norm();
        }
    }

    std::cout << "\nComparison with CMTM computation: ";
    if (VERBOSITY) {
        std::cout << "\nMax jerk norm error: " << numErr.col(0).maxCoeff()
                  << "\nMax snap norm error: " << numErr.col(1).maxCoeff()
                  << "\nMax crackle norm error: " << numErr.col(2).maxCoeff();
    } else {
        std::cout << "Max norm error = " << numErr.maxCoeff();
    }
}

template <typename MI, typename Tree>
void checkFKVecAD(MI& info, Tree& tree)
{
    const auto& mb = info.model.mb;
    auto& mbc = info.model.mbc;
    int t = info.nt / 2;
    mbc.gravity.setZero(); // DO NOT USE GRAVITY HERE

    // Compute CMTM FK
    tree.init(t, info);
    FK(info, tree);
    // RBDyn fk
    mbc.q = rbd::vectorToParam(mb, info.q.col(t));
    mbc.alpha = rbd::vectorToDof(mb, info.dqs[0].col(t));
    mbc.alphaD = rbd::vectorToDof(mb, info.dqs[1].col(t));
    rbd::forwardKinematics(mb, mbc);
    rbd::forwardVelocity(mb, mbc);
    rbd::forwardAcceleration(mb, mbc);

    using CppAD::AD; // use AD as abbreviation for CppAD::AD
    using std::vector; // use vector as abbreviation for std::vector

    // domain space vector
    vector<AD<AD<double>>> X(2 * mb.nrDof()); // vector of domain space variables
    vector<AD<AD<double>>> Y(6);
    std::vector<JointAD<AD<AD<double>>>> jAD;
    for (int i = 0; i < mb.nrDof(); ++i) {
        X[i] = info.q(i, t);
        X[i + mb.nrDof()] = info.dqs[0](i, t);
    }

    // declare independent variables and start recording operation sequence
    CppAD::Independent(X);

    // range space vector
    Y = FKVecAD(X, info.model, jAD); // value during recording of operations

    // store operation sequence in f: X -> Y and stop recording
    CppAD::ADFun<AD<double>> f(X, Y);

    // compute derivative using operation sequence stored in f
    vector<AD<double>> jac(6 * 2 * mb.nrDof()); // Jacobian of f (m by n matrix)
    vector<AD<double>> X2(2 * mb.nrDof()); // domain space vector
    for (int i = 0; i < mb.nrDof(); ++i) {
        X2[i] = info.q(i, t);
        X2[i + mb.nrDof()] = info.dqs[0](i, t);
    }

    // declare independent variables and start recording operation sequence
    CppAD::Independent(X2);

    jac = f.Jacobian(X2); // Jacobian for operation sequence
    auto jacT = TransposeAD(jac, 6);

    // store operation sequence in f: X -> Y and stop recording
    CppAD::ADFun<double> f2(X2, jacT);

    // compute 2nd-order derivative using operation sequence stored in f
    vector<double> H(6 * 2 * mb.nrDof() * 2 * mb.nrDof());
    vector<double> x(2 * mb.nrDof()); // domain space vector
    vector<double> dx(2 * mb.nrDof()); // domain space vector
    vector<double> ddx(2 * mb.nrDof()); // domain space vector
    for (int i = 0; i < mb.nrDof(); ++i) {
        x[i] = info.q(i, t);
        x[i + mb.nrDof()] = info.dqs[0](i, t);
        dx[i] = info.dqs[0](i, t);
        dx[i + mb.nrDof()] = info.dqs[1](i, t);
        ddx[i] = info.dqs[1](i, t);
        ddx[i + mb.nrDof()] = info.dqs[2](i, t);
    }
    H = f2.Jacobian(x);

    auto jacTNoAD = CastOutAD(jacT);
    auto HT = TransposeAD(H, 6 * 2 * mb.nrDof());
    auto Hdq = ProductAD(HT, dx, 6 * 2 * mb.nrDof());
    auto Hdqdq = ProductAD(Hdq, dx, 6);
    auto Jddq = ProductAD(jacTNoAD, ddx, 6);

    // DX = J * dQ
    auto DX = ProductAD(jacTNoAD, dx, 6);
    // D2X = J * d2Q + dJ * DQ
    auto D2X = SumAD(Jddq, Hdqdq);

    // print the results
    auto print = [](const char* title, const auto& v) {
        std::cout << title << " = [";
        for (decltype(v.size()) i = 0; i < v.size(); ++i)
            std::cout << " " << v[i];

        std::cout << " ]";
    };

    if (VERBOSITY) {
        Eigen::Vector6d v = tree.links[mb.nrBodies() - 1].motion()[0].vector();
        auto ADv = CastOutAD(CastOutAD(Y));
        print("\nAuto-diff:    vel", Y);
        print("\nRBDyn result: vel", mbc.bodyVelB[mb.nrBodies() - 1].vector());
        print("\nCMTM result:  vel", v);
        std::cout << "\nNorm error:   vel = " << (v - Eigen::Map<Eigen::Vector6d>(ADv.data())).norm();
        v = tree.links[mb.nrBodies() - 1].motion()[1].vector();
        print("\nAuto-diff:    acc", DX);
        print("\nRBDyn result: acc", mbc.bodyAccB[mb.nrBodies() - 1].vector());
        print("\nCMTM result:  acc", v);
        std::cout << "\nNorm error:   acc = " << (v - Eigen::Map<Eigen::Vector6d>(DX.data())).norm();
        v = tree.links[mb.nrBodies() - 1].motion()[2].vector();
        print("\nAuto-diff:   jerk", D2X);
        print("\nCMTM result: jerk", v);
        std::cout << "\nNorm error:  jerk = " << (v - Eigen::Map<Eigen::Vector6d>(D2X.data())).norm();
        std::cout << "\n[nrDof, Jacobian size, Hessian size]: ["
                  << mb.nrDof() << ", " << jac.size() << ", " << H.size() << "]";
    } else {
        auto jerk = tree.links[mb.nrBodies() - 1].motion()[2].vector();
        std::cout << "\nJerk test: " << (jerk - Eigen::Map<Eigen::Vector6d>(D2X.data())).norm();
    }
}

template <typename MI, typename Tree>
void checkFKAccAD(MI& info, Tree& tree)
{
    const auto& mb = info.model.mb;
    auto& mbc = info.model.mbc;
    int t = info.nt / 2;
    mbc.gravity.setZero(); // DO NOT USE GRAVITY HERE

    // Compute CMTM FK
    tree.init(t, info);
    FK(info, tree);
    // RBDyn fk
    mbc.q = rbd::vectorToParam(mb, info.q.col(t));
    mbc.alpha = rbd::vectorToDof(mb, info.dqs[0].col(t));
    mbc.alphaD = rbd::vectorToDof(mb, info.dqs[1].col(t));
    rbd::forwardKinematics(mb, mbc);
    rbd::forwardVelocity(mb, mbc);
    rbd::forwardAcceleration(mb, mbc);

    using CppAD::AD; // use AD as abbreviation for CppAD::AD
    using std::vector; // use vector as abbreviation for std::vector

    // domain space vector
    vector<AD<AD<double>>> X(3 * mb.nrDof()); // vector of domain space variables
    vector<AD<AD<double>>> Y(6);
    std::vector<JointAD<AD<AD<double>>>> jAD;
    for (int i = 0; i < mb.nrDof(); ++i) {
        X[i] = info.q(i, t);
        X[i + mb.nrDof()] = info.dqs[0](i, t);
        X[i + 2 * mb.nrDof()] = info.dqs[1](i, t);
    }

    // declare independent variables and start recording operation sequence
    CppAD::Independent(X);

    // range space vector
    Y = FKAccAD(X, info.model, jAD); // value during recording of operations

    // store operation sequence in f: X -> Y and stop recording
    CppAD::ADFun<AD<double>> f(X, Y);

    // compute derivative using operation sequence stored in f
    vector<AD<double>> jac(6 * 3 * mb.nrDof()); // Jacobian of f (m by n matrix)
    vector<AD<double>> X2(3 * mb.nrDof()); // domain space vector
    for (int i = 0; i < mb.nrDof(); ++i) {
        X2[i] = info.q(i, t);
        X2[i + mb.nrDof()] = info.dqs[0](i, t);
        X2[i + 2 * mb.nrDof()] = info.dqs[1](i, t);
    }

    // declare independent variables and start recording operation sequence
    CppAD::Independent(X2);

    jac = f.Jacobian(X2); // Jacobian for operation sequence
    auto jacT = TransposeAD(jac, 6);

    // store operation sequence in f: X -> Y and stop recording
    CppAD::ADFun<double> f2(X2, jacT);

    // compute 2nd-order derivative using operation sequence stored in f
    vector<double> H(6 * 3 * mb.nrDof() * 3 * mb.nrDof());
    vector<double> x(3 * mb.nrDof()); // domain space vector
    vector<double> dx(3 * mb.nrDof()); // domain space vector
    vector<double> ddx(3 * mb.nrDof()); // domain space vector
    for (int i = 0; i < mb.nrDof(); ++i) {
        x[i] = info.q(i, t);
        x[i + mb.nrDof()] = info.dqs[0](i, t);
        x[i + 2 * mb.nrDof()] = info.dqs[1](i, t);
        dx[i] = info.dqs[0](i, t);
        dx[i + mb.nrDof()] = info.dqs[1](i, t);
        dx[i + 2 * mb.nrDof()] = info.dqs[2](i, t);
        ddx[i] = info.dqs[1](i, t);
        ddx[i + mb.nrDof()] = info.dqs[2](i, t);
        ddx[i + 2 * mb.nrDof()] = info.dqs[3](i, t);
    }
    H = f2.Jacobian(x);

    auto jacTNoAD = CastOutAD(jacT);
    auto HT = TransposeAD(H, 6 * 3 * mb.nrDof());
    auto Hdq = ProductAD(HT, dx, 6 * 3 * mb.nrDof());
    auto Hdqdq = ProductAD(Hdq, dx, 6);
    auto Jddq = ProductAD(jacTNoAD, ddx, 6);

    // DX = J * dQ
    auto DX = ProductAD(jacTNoAD, dx, 6);
    // D2X = J * d2Q + dJ * DQ
    auto D2X = SumAD(Jddq, Hdqdq);

    // print the results
    auto print = [](const char* title, const auto& v) {
        std::cout << title << " = [";
        for (decltype(v.size()) i = 0; i < v.size(); ++i)
            std::cout << " " << v[i];

        std::cout << " ]";
    };

    if (VERBOSITY) {
        Eigen::Vector6d v = tree.links[mb.nrBodies() - 1].motion()[1].vector();
        auto ADv = CastOutAD(CastOutAD(Y));
        print("\nAuto-diff:    acc", Y);
        print("\nRBDyn result: acc", mbc.bodyAccB[mb.nrBodies() - 1].vector());
        print("\nCMTM result:  acc", v);
        std::cout << "\nNorm error:   acc = " << (v - Eigen::Map<Eigen::Vector6d>(ADv.data())).norm();
        v = tree.links[mb.nrBodies() - 1].motion()[2].vector();
        print("\nAuto-diff:   jerk", DX);
        // print("\nRBDyn result: acc", mbc.bodyAccB[mb.nrBodies() - 1].vector());
        print("\nCMTM result: jerk", v);
        std::cout << "\nNorm error:  jerk = " << (v - Eigen::Map<Eigen::Vector6d>(DX.data())).norm();
        v = tree.links[mb.nrBodies() - 1].motion()[3].vector();
        print("\nAuto-diff:   snap", D2X);
        print("\nCMTM result: snap", v);
        std::cout << "\nNorm error:  snap = " << (v - Eigen::Map<Eigen::Vector6d>(D2X.data())).norm();
        std::cout << "\n[nrDof, Jacobian size, Hessian size]: ["
                  << mb.nrDof() << ", " << jac.size() << ", " << H.size() << "]";
    } else {
        auto jerk = tree.links[mb.nrBodies() - 1].motion()[2].vector();
        std::cout << "\nJerk test: " << (jerk - Eigen::Map<Eigen::Vector6d>(D2X.data())).norm();
    }
}

// template <typename MI, typename Tree>
// void checkEigenFKAccAD(MI& info, Tree& tree)
// {
//     const auto& mb = info.model.mb;
//     auto& mbc = info.model.mbc;
//     int t = info.nt / 2;
//     mbc.gravity.setZero(); // DO NOT USE GRAVITY HERE

//     // Compute CMTM FK
//     tree.init(t, info);
//     FK(info, tree);
//     // RBDyn fk
//     mbc.q = rbd::vectorToParam(mb, info.q.col(t));
//     mbc.alpha = rbd::vectorToDof(mb, info.dqs[0].col(t));
//     mbc.alphaD = rbd::vectorToDof(mb, info.dqs[1].col(t));
//     rbd::forwardKinematics(mb, mbc);
//     rbd::forwardVelocity(mb, mbc);
//     rbd::forwardAcceleration(mb, mbc);

//     using CppAD::AD; // use AD as abbreviation for CppAD::AD

//     // domain space vector
//     Eigen::Matrix<AD<AD<double>>, Eigen::Dynamic, 1> X(2 * mb.nrDof()); // vector of domain space variables
//     // Eigen::Matrix<AD<AD<double>>, Eigen::Dynamic, 1> X(3 * mb.nrDof()); // vector of domain space variables
//     Eigen::Matrix<AD<AD<double>>, Eigen::Dynamic, 1> Y(6);
//     std::vector<EigenJointAD<AD<AD<double>>>> jAD;
//     for (int i = 0; i < mb.nrDof(); ++i) {
//         X(i) = info.q(i, t);
//         X(i + mb.nrDof()) = info.dqs[0](i, t);
//         // X(i + 2 * mb.nrDof()) = info.dqs[1](i, t);
//     }

//     // declare independent variables and start recording operation sequence
//     CppAD::Independent(X);

//     // range space vector
//     Y = EigenFKAD(X, info.model, jAD); // value during recording of operations

//     // store operation sequence in f: X -> Y and stop recording
//     CppAD::ADFun<AD<double>> f(X, Y);

//     // compute derivative using operation sequence stored in f
//     Eigen::Matrix<AD<double>, Eigen::Dynamic, 1> jac(6 * 2 * mb.nrDof()); // Jacobian of f (m by n matrix)
//     Eigen::Matrix<AD<double>, Eigen::Dynamic, 1> X2(2 * mb.nrDof()); // domain space vector
//     // Eigen::Matrix<AD<double>, Eigen::Dynamic, 1> jac(6 * 3 * mb.nrDof()); // Jacobian of f (m by n matrix)
//     // Eigen::Matrix<AD<double>, Eigen::Dynamic, 1> X2(3 * mb.nrDof()); // domain space vector
//     for (int i = 0; i < mb.nrDof(); ++i) {
//         X2(i) = info.q(i, t);
//         X2(i + mb.nrDof()) = info.dqs[0](i, t);
//         // X2(i + 2 * mb.nrDof()) = info.dqs[1](i, t);
//     }

//     // declare independent variables and start recording operation sequence
//     CppAD::Independent(X2);

//     jac = f.Jacobian(X2); // Jacobian for operation sequence
//     // auto jacT = TransposeAD(jac, 6);

//     // store operation sequence in f: X -> Y and stop recording
//     CppAD::ADFun<double> f2(X2, jac);

//     // compute 2nd-order derivative using operation sequence stored in f
//     Eigen::Matrix<double, Eigen::Dynamic, 1> H(6 * 2 * mb.nrDof() * 2 * mb.nrDof());
//     Eigen::Matrix<double, Eigen::Dynamic, 1> x(2 * mb.nrDof()); // domain space vector
//     Eigen::Matrix<double, Eigen::Dynamic, 1> dx(2 * mb.nrDof()); // domain space vector
//     Eigen::Matrix<double, Eigen::Dynamic, 1> ddx(2 * mb.nrDof()); // domain space vector
//     // Eigen::Matrix<double, Eigen::Dynamic, 1> H(6 * 3 * mb.nrDof() * 3 * mb.nrDof());
//     // Eigen::Matrix<double, Eigen::Dynamic, 1> x(3 * mb.nrDof()); // domain space vector
//     // Eigen::Matrix<double, Eigen::Dynamic, 1> dx(3 * mb.nrDof()); // domain space vector
//     // Eigen::Matrix<double, Eigen::Dynamic, 1> ddx(3 * mb.nrDof()); // domain space vector
//     for (int i = 0; i < mb.nrDof(); ++i) {
//         x(i) = info.q(i, t);
//         x(i + mb.nrDof()) = info.dqs[0](i, t);
//         // x(i + 2 * mb.nrDof()) = info.dqs[1](i, t);
//         dx(i) = info.dqs[0](i, t);
//         dx(i + mb.nrDof()) = info.dqs[1](i, t);
//         // dx(i + 2 * mb.nrDof()) = info.dqs[2](i, t);
//         ddx(i) = info.dqs[1](i, t);
//         ddx(i + mb.nrDof()) = info.dqs[2](i, t);
//         // ddx(i + 2 * mb.nrDof()) = info.dqs[3](i, t);
//     }
//     H = f2.Jacobian(x);

//     Eigen::Matrix<double, 6, Eigen::Dynamic, Eigen::RowMajor> jacTNoAD = Eigen::Map<Eigen::Matrix<double, 6, Eigen::Dynamic, Eigen::RowMajor>>(H.data(), 6, 2 * mb.nrDof() * 2 * mb.nrDof());
//     // Eigen::Matrix<double, 6, Eigen::Dynamic, Eigen::RowMajor> jacTNoAD = Eigen::Map<Eigen::Matrix<double, 6, Eigen::Dynamic, Eigen::RowMajor>>(H.data(), 6, 3 * mb.nrDof() * 3 * mb.nrDof());
//     Eigen::Matrix<double, 6, Eigen::Dynamic, Eigen::RowMajor> Hdq = H * dx;
//     Eigen::Matrix<double, 6, Eigen::Dynamic, Eigen::RowMajor> Hdqdq = Hdq * dx;
//     Eigen::Matrix<double, 6, Eigen::Dynamic, Eigen::RowMajor> Jddq = jacTNoAD * ddx;

//     // DX = J * dQ
//     auto DX = jacTNoAD * dx;
//     // D2X = J * d2Q + dJ * DQ
//     auto D2X = Jddq + Hdqdq;

//     // print the results
//     auto print = [](const char* title, const auto& v) {
//         std::cout << title << " = [";
//         for (decltype(v.size()) i = 0; i < v.size(); ++i)
//             std::cout << " " << v[i];

//         std::cout << " ]";
//     };

//     if (VERBOSITY) {
//         Eigen::Vector6d v = tree.links[mb.nrBodies() - 1].motion()[0].vector();
//         // Eigen::Vector6d v = tree.links[mb.nrBodies() - 1].motion()[1].vector();
//         auto ADv = CastOutAD(CastOutAD(Y));
//         print("\nAuto-diff:    acc", Y);
//         print("\nRBDyn result: acc", mbc.bodyAccB[mb.nrBodies() - 1].vector());
//         print("\nCMTM result:  acc", v);
//         std::cout << "\nNorm error:   acc = " << (v - Eigen::Map<Eigen::Vector6d>(ADv.data())).norm();
//         v = tree.links[mb.nrBodies() - 1].motion()[2].vector();
//         // v = tree.links[mb.nrBodies() - 1].motion()[2].vector();
//         print("\nAuto-diff:   jerk", DX);
//         // print("\nRBDyn result: acc", mbc.bodyAccB[mb.nrBodies() - 1].vector());
//         print("\nCMTM result: jerk", v);
//         std::cout << "\nNorm error:  jerk = " << (v - Eigen::Map<Eigen::Vector6d>(DX.data())).norm();
//         v = tree.links[mb.nrBodies() - 1].motion()[2].vector();
//         print("\nAuto-diff:   snap", D2X);
//         print("\nCMTM result: snap", v);
//         std::cout << "\nNorm error:  snap = " << (v - Eigen::Map<Eigen::Vector6d>(D2X.data())).norm();
//         std::cout << "\n[nrDof, Jacobian size, Hessian size]: ["
//                   << mb.nrDof() << ", " << jac.size() << ", " << H.size() << "]";
//     } else {
//         auto jerk = tree.links[mb.nrBodies() - 1].motion()[2].vector();
//         std::cout << "\nJerk test: " << (jerk - Eigen::Map<Eigen::Vector6d>(D2X.data())).norm();
//     }
// }

template <size_t Order, typename Fun>
double testSpeed(Fun fun, const std::string& funName, int nLoop, const std::string& model = "human", int nJoints = 5)
{
    std::cout << "\nBenching " << funName << " (order " << Order << ")";
    ModelInfo<Order> info(101, 1e-8);
    info.init(model, nJoints);
    Tree6d<Order> tree;

    const auto& mb = info.model.mb;
    auto& mbc = info.model.mbc;
    int t = info.nt / 2;
    mbc.gravity.setZero();
    tree.init(t, info);

    double mean = 0.;
    for (int i = 0; i < nLoop; ++i) {
        auto start = std::chrono::high_resolution_clock::now();
        fun(info, tree);
        auto end = std::chrono::high_resolution_clock::now();
        mean += std::chrono::duration<double>(end - start).count();
    }

    return mean / static_cast<double>(nLoop);
}

template <size_t Order>
double testADSpeed(int nLoop, int nrJoints)
{
    ModelInfo<Order> info(101, 1e-8);
    info.init("manipulator", nrJoints);
    const auto& mb = info.model.mb;
    auto& mbc = info.model.mbc;
    int t = info.nt / 2;
    mbc.gravity.setZero();

    using CppAD::AD; // use AD as abbreviation for CppAD::AD
    using std::vector; // use vector as abbreviation for std::vector

    // domain space vector
    vector<AD<AD<double>>> X(2 * mb.nrDof()); // vector of domain space variables
    vector<AD<AD<double>>> Y(6);
    // compute derivative using operation sequence stored in f
    vector<AD<double>> jac(6 * 2 * mb.nrDof()); // Jacobian of f (m by n matrix)
    vector<AD<double>> X2(2 * mb.nrDof()); // domain space vector
    // compute 2nd-order derivative using operation sequence stored in f
    vector<double> H(6 * 2 * mb.nrDof() * 2 * mb.nrDof());
    vector<double> x(2 * mb.nrDof()); // domain space vector
    vector<double> dx(2 * mb.nrDof()); // domain space vector
    vector<double> ddx(2 * mb.nrDof()); // domain space vector

    std::vector<JointAD<AD<AD<double>>>> jAD;
    for (int i = 0; i < mb.nrDof(); ++i) {
        X[i] = info.q(i, t);
        X[i + mb.nrDof()] = info.dqs[0](i, t);
        X2[i] = info.q(i, t);
        X2[i + mb.nrDof()] = info.dqs[0](i, t);
        x[i] = info.q(i, t);
        x[i + mb.nrDof()] = info.dqs[0](i, t);
        dx[i] = info.dqs[0](i, t);
        dx[i + mb.nrDof()] = info.dqs[1](i, t);
        ddx[i] = info.dqs[1](i, t);
        ddx[i + mb.nrDof()] = info.dqs[2](i, t);
    }

    double mean = 0.;
    std::cout << "\nTest AD nrJoints=" << nrJoints << ", nLoop=" << nLoop;
    int per = 0;
    for (int i = 0; i < nLoop; ++i) {
        if (nLoop >= 5 && i % (nLoop / 5) == 0) {
            std ::cout << " - " << per << "% ";
            per += 20;
        }

        auto start = std::chrono::high_resolution_clock::now();
        // declare independent variables and start recording operation sequence
        CppAD::Independent(X);

        // range space vector
        Y = FKAD(X, info.model, jAD); // value during recording of operations

        // store operation sequence in f: X -> Y and stop recording
        CppAD::ADFun<AD<double>> f(X, Y);

        // declare independent variables and start recording operation sequence
        CppAD::Independent(X2);

        // Get Jacobian
        jac = f.Jacobian(X2); // Jacobian for operation sequence
        auto jacT = TransposeAD(jac, 6);

        // store operation sequence in f: X -> Y and stop recording
        CppAD::ADFun<double> f2(X2, jacT);
        // Get Hessian
        H = f2.Jacobian(x);

        auto jacTNoAD = CastOutAD(jacT);
        auto HT = TransposeAD(H, 6 * 2 * mb.nrDof());
        auto Hdq = ProductAD(HT, dx, 6 * 2 * mb.nrDof());
        auto Hdqdq = ProductAD(Hdq, dx, 6);
        auto Jddq = ProductAD(jacTNoAD, ddx, 6);

        // DX = J * dQ
        auto DX = ProductAD(jacTNoAD, dx, 6);
        // D2X = J * d2Q + dJ * DQ
        auto D2X = SumAD(Jddq, Hdqdq);
        auto end = std::chrono::high_resolution_clock::now();
        mean += std::chrono::duration<double>(end - start).count();
    }

    return mean / static_cast<double>(nLoop);
}

template <size_t Order>
static void BM_rbd_human_FK(benchmark::State& state) {
    ModelInfo<Order> info(101, 1e-8);
    info.init("human");
    Tree6d<Order> tree;

    const auto& mb = info.model.mb;
    auto& mbc = info.model.mbc;
    int t = info.nt / 2;
    mbc.gravity.setZero();
    tree.init(t, info);

    for (auto _ : state) {
        rbd::forwardKinematics(mb, mbc);
        if constexpr (Order >= 1) {
            rbd::forwardVelocity(mb, mbc);
        }
        if constexpr (Order >= 2) {
            rbd::forwardAcceleration(mb, mbc);
        }
    }
}

template <size_t Order>
static void BM_cmtm_human_FK(benchmark::State& state) {
    ModelInfo<Order> info(101, 1e-8);
    info.init("human");
    Tree6d<Order> tree;

    const auto& mb = info.model.mb;
    auto& mbc = info.model.mbc;
    int t = info.nt / 2;
    mbc.gravity.setZero();
    tree.init(t, info);

    for (auto _ : state) {
        FK(info, tree);
    }
}

template <size_t Order>
static void BM_rbd_manipulator_FK(benchmark::State& state) {
    ModelInfo<Order> info(101, 1e-8);
    info.init("manipulator", state.range(0));
    Tree6d<Order> tree;

    const auto& mb = info.model.mb;
    auto& mbc = info.model.mbc;
    int t = info.nt / 2;
    mbc.gravity.setZero();
    tree.init(t, info);

    for (auto _ : state) {
        rbd::forwardKinematics(mb, mbc);
        if constexpr (Order >= 1) {
            rbd::forwardVelocity(mb, mbc);
        }
        if constexpr (Order >= 2) {
            rbd::forwardAcceleration(mb, mbc);
        }
    }
}

template <size_t Order>
static void BM_cmtm_manipulator_FK(benchmark::State& state) {
    ModelInfo<Order> info(101, 1e-8);
    info.init("manipulator", state.range(0));
    Tree6d<Order> tree;

    const auto& mb = info.model.mb;
    auto& mbc = info.model.mbc;
    int t = info.nt / 2;
    mbc.gravity.setZero();
    tree.init(t, info);

    for (auto _ : state) {
        FK(info, tree);
    }
}

static void BM_rbd_human_J(benchmark::State& state) {
    ModelInfo<2> info(101, 1e-8);
    info.init("human");
    Tree6d<2> tree;
    std::string bodyName = "RARM0";

    const auto& mb = info.model.mb;
    auto& mbc = info.model.mbc;
    int t = info.nt / 2;
    mbc.gravity.setZero();
    tree.init(t, info);
    rbd::forwardKinematics(mb, mbc);
    rbd::forwardVelocity(mb, mbc);
    rbd::forwardAcceleration(mb, mbc);

    for (auto _ : state) {
        Eigen::MatrixXd bj, bjd;
        rbd::Jacobian J(info.model.mb, bodyName);
        benchmark::DoNotOptimize(bj = J.bodyJacobian(info.model.mb, info.model.mbc));
        benchmark::DoNotOptimize(bjd = J.bodyJacobianDot(info.model.mb, info.model.mbc));
        benchmark::ClobberMemory();
    }
}

static void BM_cmtm_human_J(benchmark::State& state) {
    ModelInfo<2> info(101, 1e-8);
    info.init("human");
    Tree6d<2> tree;
    std::string bodyName = "RARM0";

    const auto& mb = info.model.mb;
    auto& mbc = info.model.mbc;
    int t = info.nt / 2;
    mbc.gravity.setZero();
    tree.init(t, info);
    rbd::forwardKinematics(mb, mbc);
    rbd::forwardVelocity(mb, mbc);
    rbd::forwardAcceleration(mb, mbc);

    for (auto _ : state) {
        Eigen::MatrixXd bj, bjd;
        benchmark::DoNotOptimize(bj = MotionJacobianOfOrder<0>(bodyName, info, tree));
        benchmark::DoNotOptimize(bjd = MotionJacobianOfOrder<1>(bodyName, info, tree));
        benchmark::ClobberMemory();
    }
}

static void BM_rbd_manipulator_J(benchmark::State& state) {
    ModelInfo<2> info(101, 1e-8);
    auto nDof = state.range(0);
    info.init("manipulator", nDof);
    Tree6d<2> tree;
    std::string bodyName = "b" + std::to_string(nDof);

    const auto& mb = info.model.mb;
    auto& mbc = info.model.mbc;
    int t = info.nt / 2;
    mbc.gravity.setZero();
    tree.init(t, info);
    rbd::forwardKinematics(mb, mbc);
    rbd::forwardVelocity(mb, mbc);
    rbd::forwardAcceleration(mb, mbc);

    for (auto _ : state) {
        rbd::Jacobian J;
        Eigen::MatrixXd bj, bjd;
        benchmark::DoNotOptimize(J = rbd::Jacobian(info.model.mb, bodyName));
        benchmark::DoNotOptimize(bj = J.bodyJacobian(info.model.mb, info.model.mbc));
        benchmark::DoNotOptimize(bjd = J.bodyJacobianDot(info.model.mb, info.model.mbc));
        benchmark::ClobberMemory();
    }
}

static void BM_cmtm_manipulator_J(benchmark::State& state) {
    ModelInfo<2> info(101, 1e-8);
    auto nDof = state.range(0);
    info.init("manipulator", nDof);
    Tree6d<2> tree;
    std::string bodyName = "b" + std::to_string(nDof);

    const auto& mb = info.model.mb;
    auto& mbc = info.model.mbc;
    int t = info.nt / 2;
    mbc.gravity.setZero();
    tree.init(t, info);
    rbd::forwardKinematics(mb, mbc);
    rbd::forwardVelocity(mb, mbc);
    rbd::forwardAcceleration(mb, mbc);

    for (auto _ : state) {
        Eigen::MatrixXd bj, bjd;
        benchmark::DoNotOptimize(bj = MotionJacobianOfOrder<0>(bodyName, info, tree));
        benchmark::DoNotOptimize(bjd = MotionJacobianOfOrder<1>(bodyName, info, tree));
        benchmark::ClobberMemory();
    }
}

static void BM_AD_manipulator_FK(benchmark::State& state)
{
    ModelInfo<3> info(101, 1e-8);
    auto nDof = state.range(0);
    info.init("manipulator", nDof);
    const auto& mb = info.model.mb;
    auto& mbc = info.model.mbc;
    int t = info.nt / 2;
    mbc.gravity.setZero();

    using CppAD::AD; // use AD as abbreviation for CppAD::AD

    // domain space vector
    Eigen::Matrix<AD<AD<double>>, Eigen::Dynamic, 1> X(2 * mb.nrDof()); // vector of domain space variables
    Eigen::Matrix<AD<AD<double>>, Eigen::Dynamic, 1> Y(6);
    // compute derivative using operation sequence stored in f
    Eigen::Matrix<AD<double>, 6, Eigen::Dynamic> jac(6, 2 * mb.nrDof()); // Jacobian of f (m by n matrix)
    Eigen::Matrix<AD<double>, Eigen::Dynamic, 1> jacLin(6 * 2 * mb.nrDof()); // Jacobian of f (m by n matrix)
    Eigen::Matrix<AD<double>, Eigen::Dynamic, 1> X2(2 * mb.nrDof()); // domain space vector
    // compute 2nd-order derivative using operation sequence stored in f
    Eigen::Matrix<double, Eigen::Dynamic, 1> HLin(6 * 2 * mb.nrDof() * 2 * mb.nrDof());
    Eigen::Matrix<double, 6, Eigen::Dynamic> dJ(6, 2 * mb.nrDof());
    Eigen::Matrix<double, Eigen::Dynamic, 1> dJLin(6 * 2 * mb.nrDof(), 1);
    Eigen::Matrix<double, Eigen::Dynamic, 1> x(2 * mb.nrDof()); // domain space vector
    Eigen::Matrix<double, Eigen::Dynamic, 1> dx(2 * mb.nrDof()); // domain space vector

    std::vector<EigenJointAD<AD<AD<double>>>> jAD;
    for (int i = 0; i < mb.nrDof(); ++i) {
        X[i] = info.q(i, t);
        X[i + mb.nrDof()] = info.dqs[0](i, t);
        X2[i] = info.q(i, t);
        X2[i + mb.nrDof()] = info.dqs[0](i, t);
        dx[i] = info.dqs[0](i, t);
        dx[i + mb.nrDof()] = info.dqs[1](i, t);
    }

    // Test as-if it were col-major
    for (auto _ : state) {
        // declare independent variables and start recording operation sequence
        CppAD::Independent(X);

        // range space vector
        Y = EigenFKAD(X, info.model, jAD); // value during recording of operations

        // store operation sequence in f: X -> Y and stop recording
        CppAD::ADFun<AD<double>> f(X, Y);

        // declare independent variables and start recording operation sequence
        CppAD::Independent(X2);

        // Get Jacobian
        jac = f.Jacobian(X2); // Jacobian for operation sequence
        jacLin = Eigen::Map<Eigen::Matrix<AD<double>, Eigen::Dynamic, 1>>(jac.data(), 6 * 2 * mb.nrDof(), 1);

        // store operation sequence in f: X -> Y and stop recording
        CppAD::ADFun<double> f2(X2, jacLin);

        // Get Hessian
        HLin = f2.Jacobian(x);
        dJLin = Eigen::Map<Eigen::MatrixXd>(HLin.data(), 6 * 2 * mb.nrDof(), 2 * mb.nrDof());
        dJ = Eigen::Map<Eigen::VectorXd>(dJLin.data(), 6, 2 * mb.nrDof());

        benchmark::ClobberMemory();
    }
}

struct CheckNum {
    enum : unsigned {
        FK = 1, // Forward Kinematics
        ID = 1 << 1, // Inverse Dynamics
        FD = 1 << 2, // Forward Dynamics
        J = 1 << 3, // Jacobians
        NUM = 1 << 4, // 2nd order Numerical differentiation
        AD = 1 << 5, // Automatic Differentiation
        BENCH = 1 << 6, // Automatic Differentiation
    };
};

int main(int argc, char** argv)
{
    char newLine = '\0';
    unsigned checkNum = 0;
    std::string firstNewLines = "";
    int benchArgc = 1;
    char** benchArgv = new char*[argc];
    auto copyArg = [benchArgv, argv](size_t to_index, size_t from_index) {
        auto length = strlen(argv[from_index]) + 1;
        benchArgv[to_index] = new char[length];
        memcpy(benchArgv[to_index], argv[from_index], length);
    };
    copyArg(0, 0);
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            std::cout << "Call with ./algo [options] where [options] can be"
                      << "\n- -v or --verbose for more verbosity"
                      << "\n- One of the followings: FK, ID, FD, J, NUM" << std::endl;
            delete[] benchArgv;
            return 0;
        }
        if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
            VERBOSITY = true;
            newLine = '\n';
        } else if (strcmp(argv[i], "FK") == 0) {
            checkNum |= CheckNum::FK;
        } else if (strcmp(argv[i], "ID") == 0) {
            checkNum |= CheckNum::ID;
        } else if (strcmp(argv[i], "J") == 0) {
            checkNum |= CheckNum::J;
        } else if (strcmp(argv[i], "FD") == 0) {
            checkNum |= CheckNum::FD;
        } else if (strcmp(argv[i], "NUM") == 0) {
            checkNum |= CheckNum::NUM;
        } else if (strcmp(argv[i], "AD") == 0) {
            checkNum |= CheckNum::AD;
        } else if (strcmp(argv[i], "BENCH") == 0) {
            checkNum |= CheckNum::BENCH;
        } else {
            copyArg(benchArgc++, i);
        }
    }
    if (!checkNum) {
        checkNum = -1;
        checkNum ^= CheckNum::BENCH;
    }
    if (!(checkNum ^ CheckNum::BENCH)) {
        checkNum |= CheckNum::AD;
        checkNum |= CheckNum::FK;
        checkNum |= CheckNum::J;
    }

    ModelInfo<5> info(101, 1e-8);
    info.init("human");
    if (VERBOSITY) {
        if (checkNum ^ CheckNum::AD) {
            std::cout << info;
        }
        firstNewLines = "\n\n";
    }
    std::string bodyName = "RARM0";

    Tree4d<5> tree4d;
    Tree6d<5> tree6d;

    // Check FK
    if (checkNum & CheckNum::FK && !(checkNum & CheckNum::BENCH)) {
        std::cout << firstNewLines << "Checking FK with space 4 CMTM" << newLine;
        checkFK(info, tree4d);
        std::cout << "\n\nChecking FK with space 6 CMTM" << newLine;
        checkFK(info, tree6d);
        firstNewLines = "\n\n";
    }

    // Check ID
    if (checkNum & CheckNum::ID && !(checkNum & CheckNum::BENCH)) {
        std::cout << firstNewLines << "Checking ID without gravity" << newLine;
        checkID(info, tree6d, false);
        std::cout << "\n\nChecking ID with gravity" << newLine;
        checkID(info, tree6d, true);
        firstNewLines = "\n\n";
    }

    // Check Jacobians
    if (checkNum & CheckNum::J && !(checkNum & CheckNum::BENCH)) {
        // Check motion jacobian
        std::cout << firstNewLines << "Checking Motion Jacobian for body " << bodyName << newLine;
        checkMotionJacobian(bodyName, info, tree6d);
        // Check link momentum jacobian
        std::cout << "\n\nChecking Link Momentum Jacobian for body " << bodyName << newLine;
        checkLinkMomentumJacobian(bodyName, info, tree6d);
        // Check link force jacobian
        std::cout << "\n\nChecking Link Force Jacobian for body " << bodyName << newLine;
        checkLinkForceJacobian(bodyName, info, tree6d);
        // Check joint momentum jacobian
        std::cout << "\n\nChecking Joint Momentum Jacobian for body " << bodyName << newLine;
        checkJointMomentumJacobian(bodyName, info, tree6d);
        // Check joint force jacobian
        std::cout << "\n\nChecking Joint Force Jacobian for body " << bodyName << newLine;
        checkJointForceJacobian(bodyName, info, tree6d);
        firstNewLines = "\n\n";
    }

    if (checkNum & CheckNum::FD && !(checkNum & CheckNum::BENCH)) {
        // Check FD
        std::cout << firstNewLines << "Checking FD" << newLine;
        checkFD(info, tree6d);
        firstNewLines = "\n\n";
    }

    if (checkNum & CheckNum::NUM && !(checkNum & CheckNum::BENCH)) {
        // Check double numerical differentiation
        std::cout << firstNewLines << "Checking 2nd order numerical derivative" << newLine;
        checkNum2Diff(info, tree6d);
        firstNewLines = "\n\n";
    }

    if (checkNum & CheckNum::AD && !(checkNum & CheckNum::BENCH)) {
        // Check double automatic differentiation
        info.init("manipulator"); // Simple model
        if (VERBOSITY) {
            std::cout << firstNewLines << info;
        }
        std::cout << firstNewLines << "Checking 2nd order automatic differentiation from velocity" << newLine;
        checkFKVecAD(info, tree4d);
        std::cout << firstNewLines << "Checking 2nd order automatic differentiation from acceleration" << newLine;
        checkFKAccAD(info, tree4d);
        firstNewLines = "\n\n";
    }

    // Test FK speed
    // auto printVec = [](const auto& v, const auto& prefix, const auto& suffix) {
    //     for (size_t i = 0; i < v.size(); ++i) {
    //         std::cout << "\n"
    //                   << prefix << i << ": " << v[i] << suffix;
    //     }
    // };
    // auto printWithComma = [](const auto& v, const auto& title) {
    //     std::cout << "\n"
    //               << title << ": ";
    //     for (size_t i = 0; i < v.size(); ++i) {
    //         std::cout << v[i] << ",";
    //     }
    // };
    if (checkNum & CheckNum::BENCH) {
        if (checkNum & CheckNum::FK) {
            benchmark::RegisterBenchmark("BM_rbd_human_FK_0", BM_rbd_human_FK<0>);
            benchmark::RegisterBenchmark("BM_rbd_human_FK_1", BM_rbd_human_FK<1>);
            benchmark::RegisterBenchmark("BM_rbd_human_FK_2", BM_rbd_human_FK<2>);
            benchmark::RegisterBenchmark("BM_cmtm_human_FK_0", BM_cmtm_human_FK<0>);
            benchmark::RegisterBenchmark("BM_cmtm_human_FK_1", BM_cmtm_human_FK<1>);
            benchmark::RegisterBenchmark("BM_cmtm_human_FK_2", BM_cmtm_human_FK<2>);
            benchmark::RegisterBenchmark("BM_cmtm_human_FK_3", BM_cmtm_human_FK<3>);
            benchmark::RegisterBenchmark("BM_cmtm_human_FK_4", BM_cmtm_human_FK<4>);
            benchmark::RegisterBenchmark("BM_cmtm_human_FK_5", BM_cmtm_human_FK<5>);
            benchmark::RegisterBenchmark("BM_cmtm_human_FK_6", BM_cmtm_human_FK<6>);
            benchmark::RegisterBenchmark("BM_cmtm_human_FK_7", BM_cmtm_human_FK<7>);
            benchmark::RegisterBenchmark("BM_cmtm_human_FK_8", BM_cmtm_human_FK<8>);
            benchmark::RegisterBenchmark("BM_cmtm_human_FK_9", BM_cmtm_human_FK<9>);
            benchmark::RegisterBenchmark("BM_cmtm_human_FK_10", BM_cmtm_human_FK<10>);
            benchmark::RegisterBenchmark("BM_cmtm_human_FK_11", BM_cmtm_human_FK<11>);
            benchmark::RegisterBenchmark("BM_cmtm_human_FK_12", BM_cmtm_human_FK<12>);
            benchmark::RegisterBenchmark("BM_rbd_manipulator_FK_0", BM_rbd_manipulator_FK<0>)->DenseRange(5, 200, 5);
            benchmark::RegisterBenchmark("BM_rbd_manipulator_FK_1", BM_rbd_manipulator_FK<1>)->DenseRange(5, 200, 5);
            benchmark::RegisterBenchmark("BM_rbd_manipulator_FK_2", BM_rbd_manipulator_FK<2>)->DenseRange(5, 200, 5);
            benchmark::RegisterBenchmark("BM_cmtm_manipulator_FK_0", BM_cmtm_manipulator_FK<0>)->DenseRange(5, 200, 5);
            benchmark::RegisterBenchmark("BM_cmtm_manipulator_FK_1", BM_cmtm_manipulator_FK<1>)->DenseRange(5, 200, 5);
            benchmark::RegisterBenchmark("BM_cmtm_manipulator_FK_2", BM_cmtm_manipulator_FK<2>)->DenseRange(5, 200, 5);
            benchmark::RegisterBenchmark("BM_cmtm_manipulator_FK_3", BM_cmtm_manipulator_FK<3>)->DenseRange(5, 200, 5);
            benchmark::RegisterBenchmark("BM_cmtm_manipulator_FK_4", BM_cmtm_manipulator_FK<4>)->DenseRange(5, 200, 5);
            benchmark::RegisterBenchmark("BM_cmtm_manipulator_FK_5", BM_cmtm_manipulator_FK<5>)->DenseRange(5, 200, 5);
            benchmark::RegisterBenchmark("BM_cmtm_manipulator_FK_6", BM_cmtm_manipulator_FK<6>)->DenseRange(5, 200, 5);
            benchmark::RegisterBenchmark("BM_cmtm_manipulator_FK_7", BM_cmtm_manipulator_FK<7>)->DenseRange(5, 200, 5);
            benchmark::RegisterBenchmark("BM_cmtm_manipulator_FK_8", BM_cmtm_manipulator_FK<8>)->DenseRange(5, 200, 5);
            benchmark::RegisterBenchmark("BM_cmtm_manipulator_FK_9", BM_cmtm_manipulator_FK<9>)->DenseRange(5, 200, 5);
            benchmark::RegisterBenchmark("BM_cmtm_manipulator_FK_10", BM_cmtm_manipulator_FK<10>)->DenseRange(5, 200, 5);
            benchmark::RegisterBenchmark("BM_cmtm_manipulator_FK_11", BM_cmtm_manipulator_FK<11>)->DenseRange(5, 200, 5);
            benchmark::RegisterBenchmark("BM_cmtm_manipulator_FK_12", BM_cmtm_manipulator_FK<12>)->DenseRange(5, 200, 5);
            // std::cout << "\nBenching FK";
            // std::vector<double> cmtm_fk_ms(13);
            // std::vector<double> rbd_fk_ms(3);
            // auto CMTM_FK = [](const auto& info, auto& tree) { FK(info, tree); };
            // auto RBD_FK0 = [](auto& info, auto& tree) { rbd::forwardKinematics(info.model.mb, info.model.mbc); };
            // auto RBD_FK1 = [](auto& info, auto& tree) { rbd::forwardKinematics(info.model.mb, info.model.mbc); rbd::forwardVelocity(info.model.mb, info.model.mbc); };
            // auto RBD_FK2 = [](auto& info, auto& tree) { rbd::forwardKinematics(info.model.mb, info.model.mbc); rbd::forwardVelocity(info.model.mb, info.model.mbc); rbd::forwardAcceleration(info.model.mb, info.model.mbc); };
            // cmtm_fk_ms[0] = testSpeed<0>(CMTM_FK, "CMTM FK", 100, "manipulator", 100) * 1000.;
            // cmtm_fk_ms[1] = testSpeed<1>(CMTM_FK, "CMTM FK", 100, "manipulator", 100) * 1000.;
            // cmtm_fk_ms[2] = testSpeed<2>(CMTM_FK, "CMTM FK", 100, "manipulator", 100) * 1000.;
            // cmtm_fk_ms[3] = testSpeed<3>(CMTM_FK, "CMTM FK", 100, "manipulator", 100) * 1000.;
            // cmtm_fk_ms[4] = testSpeed<4>(CMTM_FK, "CMTM FK", 100, "manipulator", 100) * 1000.;
            // cmtm_fk_ms[5] = testSpeed<5>(CMTM_FK, "CMTM FK", 100, "manipulator", 100) * 1000.;
            // cmtm_fk_ms[6] = testSpeed<6>(CMTM_FK, "CMTM FK", 100, "manipulator", 100) * 1000.;
            // cmtm_fk_ms[7] = testSpeed<7>(CMTM_FK, "CMTM FK", 100, "manipulator", 100) * 1000.;
            // cmtm_fk_ms[8] = testSpeed<8>(CMTM_FK, "CMTM FK", 100, "manipulator", 100) * 1000.;
            // cmtm_fk_ms[9] = testSpeed<9>(CMTM_FK, "CMTM FK", 100, "manipulator", 100) * 1000.;
            // cmtm_fk_ms[10] = testSpeed<10>(CMTM_FK, "CMTM FK", 100, "manipulator", 100) * 1000.;
            // cmtm_fk_ms[11] = testSpeed<11>(CMTM_FK, "CMTM FK", 100, "manipulator", 100) * 1000.;
            // cmtm_fk_ms[12] = testSpeed<12>(CMTM_FK, "CMTM FK", 100, "manipulator", 100) * 1000.;
            // rbd_fk_ms[0] = testSpeed<0>(RBD_FK0, "RBDyn FK", 100, "manipulator", 100) * 1000.;
            // rbd_fk_ms[1] = testSpeed<1>(RBD_FK1, "RBDyn FK", 100, "manipulator", 100) * 1000.;
            // rbd_fk_ms[2] = testSpeed<2>(RBD_FK2, "RBDyn FK", 100, "manipulator", 100) * 1000.;
            // printVec(cmtm_fk_ms, "FK (CMTM) order ", "ms");
            // printVec(rbd_fk_ms, "FK (RBDyn) order ", "ms");
            // printWithComma(cmtm_fk_ms, "CMTM");
            // printWithComma(rbd_fk_ms, "RBDyn");
        }
        if (checkNum & CheckNum::J) {
            benchmark::RegisterBenchmark("BM_rbd_human_J", BM_rbd_human_J);
            benchmark::RegisterBenchmark("BM_cmtm_human_J", BM_rbd_human_J);
            benchmark::RegisterBenchmark("BM_rbd_manipulator_J", BM_rbd_human_J)->DenseRange(5, 200, 5);
            benchmark::RegisterBenchmark("BM_cmtm_manipulator_J", BM_rbd_human_J)->DenseRange(5, 200, 5);

            // std::cout << "\nBenching J";
            // auto CMTM_J = [bodyName](const auto& info, auto& tree) { MotionJacobianOfOrder<0>(bodyName, info, tree); MotionJacobianOfOrder<1>(bodyName, info, tree); };
            // auto RBD_J = [bodyName](const auto& info, auto& tree) { rbd::Jacobian J(info.model.mb, bodyName); J.bodyJacobian(info.model.mb, info.model.mbc); J.bodyJacobianDot(info.model.mb, info.model.mbc); };
            // double cmtm_J = testSpeed<2>(CMTM_J, "CMTM Jacobian", 100) * 1000.;
            // double rbd_J = testSpeed<2>(RBD_J, "RBDyn Jacobian", 100) * 1000.;
            // std::cout << "\nJacobian computation (CMTM): " << cmtm_J << "ms";
            // std::cout << "\nJacobian computation (RBDyn): " << rbd_J << "ms";
        }
        if (checkNum & CheckNum::AD) {
            benchmark::RegisterBenchmark("BM_AD_manipulator_FK", BM_AD_manipulator_FK)->DenseRange(5, 200, 5);
            // std::cout << "\nBenching AD";
            // int checkSize = 40;
            // std::vector<double> cmtm_fk_ms(checkSize);
            // std::vector<double> ad_fk_ms(checkSize);
            // auto CMTM_FK = [](const auto& info, auto& tree) { FK(info, tree); };
            // for (int i = 0; i < checkSize; ++i) {
            //     ad_fk_ms[i] = testADSpeed<3>(10, 5 + i * 5) * 1000.;
            // }
            // for (int i = 0; i < checkSize; ++i) {
            //     std::cout << "\nTest CMTM nrJoints=" << 5 + i * 5;
            //     cmtm_fk_ms[i] = testSpeed<3>(CMTM_FK, "CMTM FK", 10, "manipulator", 5 + i * 5) * 1000.;
            // }
            // printVec(cmtm_fk_ms, "FK (CMTM) nrJoints 10 + 5 * ", "ms");
            // printVec(ad_fk_ms, "FK (AD) nrJoints 10 + 5 * ", "ms");
            // printWithComma(cmtm_fk_ms, "CMTM");
            // printWithComma(ad_fk_ms, "AD");
        }

        benchmark::Initialize(&benchArgc, benchArgv);
        benchmark::RunSpecifiedBenchmarks();
    }
    std::cout << std::endl;

    for (size_t i = 0; i < benchArgc; ++i) {
        delete[] benchArgv[i];
    }
    delete[] benchArgv;
    return 0;
}