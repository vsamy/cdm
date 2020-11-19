#pragma once

#include <cod/Model.hpp>
#include <cod/ModelConfig.hpp>
#include <rbdyn/EulerIntegration.h>
#include <rbdyn/MultiBody.h>
#include <rbdyn/MultiBodyConfig.h>

struct TrajectoryData {
    void setCurData(int t)
    {
        if (t < 0) {
            curData = 0;
        } else if (t > static_cast<int>(time.size())) {
            curData = static_cast<int>(time.size());
        } else {
            curData = t;
        }
    }

    int order;
    int curData;
    int nt;
    double dt;
    Eigen::VectorXd time;
    std::vector<Eigen::VectorXd> q;
    std::vector<Eigen::MatrixXd> dqs;
    Eigen::Vector3d gravity;
};

template <int Order>
TrajectoryData GenerateData(const rbd::MultiBody& mb, rbd::MultiBodyConfig& mbc, int nt, double dt = 1e-8)
{
    constexpr int order = Order > 2 ? Order : 2;
    TrajectoryData data;
    data.order = order;
    data.nt = nt;
    data.dt = dt;
    data.time = Eigen::VectorXd::LinSpaced(nt, 0., (nt - 1) * dt);
    data.gravity = Eigen::Vector3d(0, 0, 9.81);
    data.gravity.setZero();
    data.q.resize(nt, Eigen::VectorXd(mb.nrParams()));
    data.dqs.resize(nt);
    std::fill(data.dqs.begin(), data.dqs.end(), Eigen::MatrixXd(mb.nrDof(), order));

    data.q[0] = rbd::paramToVector(mb, mbc.q);
    // Compute velocity as a * sin(wt + b)
    Eigen::ArrayXd a = Eigen::ArrayXd::Random(mb.nrDof());
    Eigen::ArrayXd b = Eigen::ArrayXd::Random(mb.nrDof());
    Eigen::ArrayXd w = Eigen::ArrayXd::Random(mb.nrDof());
    for (int i = 0; i < nt; ++i) {
        size_t n = order / 2;
        for (size_t k = 0; k < n; ++k) {
            data.dqs[i].col(2 * k) = std::pow(-1., k) * a * w.pow(2 * k) * (w * data.time(i) + b).sin();
            data.dqs[i].col(2 * k + 1) = std::pow(-1., k) * a * w.pow(2 * k + 1) * (w * data.time(i) + b).cos();
        }
        if constexpr (order % 2 == 1) {
            data.dqs[i].col(order - 1) = std::pow(-1., n) * a * w.pow(2 * n) * (w * data.time(i) + b).sin();
        }
    }

    for (int i = 0; i < nt - 1; ++i) {
        auto alpha = rbd::vectorToDof(mb, data.dqs[i].col(0));
        auto alphaD = rbd::vectorToDof(mb, data.dqs[i].col(1));
        for (int j = 0; j < mb.nrJoints(); ++j)
            rbd::eulerJointIntegration(mb.joint(j).type(), alpha[j], alphaD[j], dt, mbc.q[j]);

        data.q[i + 1] = rbd::paramToVector(mb, mbc.q);
    }

    mbc.q = rbd::vectorToParam(mb, data.q[0]);
    mbc.alpha = rbd::vectorToDof(mb, data.dqs[0].col(0));
    mbc.alphaD = rbd::vectorToDof(mb, data.dqs[0].col(1));
    mbc.gravity = data.gravity;
    return data;
}

void Init(const TrajectoryData& data, const rbd::MultiBody& mb, rbd::MultiBodyConfig& mbc)
{
    for (int j = 0; j < mb.nrJoints(); ++j) {
        mbc.q = rbd::vectorToParam(mb, data.q[data.curData]);
        mbc.alpha = rbd::vectorToDof(mb, data.dqs[data.curData].col(0));
        mbc.alphaD = rbd::vectorToDof(mb, data.dqs[data.curData].col(1));
    }
}

void Init(const TrajectoryData& data, const cod::Model& m, cod::ModelConfig& mc)
{
    mc.order = data.order;
    mc.bodyMotions.resize(m.nLinks(), cod::CMTM(data.order));
    mc.jointMotions.resize(m.nLinks(), cod::CMTM(data.order));
    mc.jointMomentums.resize(m.nLinks(), cod::ForceVectorX(data.order));
    mc.jointForces.resize(m.nLinks(), cod::ForceVectorX(data.order));
    mc.bodyMomentums.resize(m.nLinks(), cod::ForceVectorX(data.order));
    mc.bodyForces.resize(m.nLinks(), cod::ForceVectorX(data.order));
    mc.jointTorques.resize(m.nLinks());
    mc.q = data.q[data.curData];
    mc.dqs = data.dqs[data.curData];

    for (int i = 0; i < m.nLinks(); ++i) {
        int p = m.jointPosInParam(i);
        cod::Transform dA;
        switch (m.joint(i).type()) {
        case cod::Joint::Type::Free: {
            dA.translation() = mc.q.template segment<3>(p + 4);
            [[fallthrough]];
        }
        case cod::Joint::Type::Spherical: {
            dA.rotation() = Eigen::Quaterniond{ mc.q(p), mc.q(p + 1), mc.q(p + 2), mc.q(p + 3) }.toRotationMatrix();
            break;
        }
        default:
            dA = expSE3(m.joint(i).S() * mc.q.segment(p, m.joint(i).dof()));
        }

        auto& jm = mc.jointMotions[i];
        jm.transform() = m.A0(i) * dA;
        for (int n = 0; n < data.order; ++n) {
            jm.motion()[n] = m.joint(i).S() * mc.dqs.col(n).segment(m.jointPosInDof(i), m.joint(i).dof());
        }

        jm.construct();
    }

    auto v = cod::MotionVectorX::Zero(data.order);
    v[1] = cod::MotionVector(Eigen::Vector3d::Zero(), data.gravity);
    mc.world.set(cod::Transform::Identity(), v);
}
