#pragma once

#include <cdm/ModelConfig.hpp>
#include <cdm/typedefs.hpp>
#include <RBDyn/EulerIntegration.h>
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>

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
    constexpr cdm::Index order = Order > 2 ? Order : 2;
    TrajectoryData data;
    data.order = order;
    data.nt = nt;
    data.dt = dt;
    data.time = Eigen::VectorXd::LinSpaced(nt, 0., (nt - 1) * dt);
    data.gravity = Eigen::Vector3d(0, 0, 9.81);
    data.q.resize(static_cast<size_t>(nt), Eigen::VectorXd(mb.nrParams()));
    data.dqs.resize(static_cast<size_t>(nt));
    std::fill(data.dqs.begin(), data.dqs.end(), Eigen::MatrixXd(mb.nrDof(), order));

    data.q[0] = rbd::paramToVector(mb, mbc.q);
    // Compute velocity as a * sin(wt + b)
    Eigen::ArrayXd a = Eigen::ArrayXd::Random(mb.nrDof());
    Eigen::ArrayXd b = Eigen::ArrayXd::Random(mb.nrDof());
    Eigen::ArrayXd w = Eigen::ArrayXd::Random(mb.nrDof());
    for (int i = 0; i < nt; ++i) {
        size_t ui = static_cast<size_t>(i);
        int n = order / 2;
        for (int k = 0; k < n; ++k) {
            data.dqs[ui].col(2 * k) = std::pow(-1., k) * a * w.pow(2 * k) * (w * data.time(i) + b).sin();
            data.dqs[ui].col(2 * k + 1) = std::pow(-1., k) * a * w.pow(2 * k + 1) * (w * data.time(i) + b).cos();
        }
        if constexpr (order % 2 == 1) {
            data.dqs[ui].col(order - 1) = std::pow(-1., n) * a * w.pow(2 * n) * (w * data.time(i) + b).sin();
        }
    }

    for (int i = 0; i < nt - 1; ++i) {
        size_t ui = static_cast<size_t>(i);
        auto alpha = rbd::vectorToDof(mb, data.dqs[ui].col(0));
        auto alphaD = rbd::vectorToDof(mb, data.dqs[ui].col(1));
        for (int j = 0; j < mb.nrJoints(); ++j) {
            size_t uj = static_cast<size_t>(j);
            rbd::eulerJointIntegration(mb.joint(j).type(), alpha[uj], alphaD[uj], dt, mbc.q[uj]);
        }

        data.q[ui + 1] = rbd::paramToVector(mb, mbc.q);
    }

    mbc.q = rbd::vectorToParam(mb, data.q[0]);
    mbc.alpha = rbd::vectorToDof(mb, data.dqs[0].col(0));
    mbc.alphaD = rbd::vectorToDof(mb, data.dqs[0].col(1));
    return data;
}

void Init(const TrajectoryData& data, const rbd::MultiBody& mb, rbd::MultiBodyConfig& mbc)
{
    mbc.gravity = data.gravity;
    for (int j = 0; j < mb.nrJoints(); ++j) {
        size_t dcd = static_cast<size_t>(data.curData);
        mbc.q = rbd::vectorToParam(mb, data.q[dcd]);
        mbc.alpha = rbd::vectorToDof(mb, data.dqs[dcd].col(0));
        mbc.alphaD = rbd::vectorToDof(mb, data.dqs[dcd].col(1));
    }
}

template <int Order>
void Init(const TrajectoryData& data, const cdm::Model& m, cdm::ModelConfig<Order>& mc)
{
    size_t un = static_cast<size_t>(m.nLinks());
    size_t dcd = static_cast<size_t>(data.curData);
    mc.bodyMotions.resize(un, cdm::CMTM<Order>(data.order));
    mc.jointMotions.resize(un, cdm::CMTM<Order>(data.order));
    mc.jointMomentums.resize(un, cdm::ForceVectorX<Order>(data.order));
    mc.jointForces.resize(un, cdm::ForceVectorX<Order>(data.order));
    mc.bodyMomentums.resize(un, cdm::ForceVectorX<Order>(data.order));
    mc.bodyForces.resize(un, cdm::ForceVectorX<Order>(data.order));
    mc.jointTorques.resize(un, Eigen::VectorXd(m.nDof()));
    mc.q = data.q[dcd];
    mc.dqs = data.dqs[dcd];

    for (cdm::Index i = 0; i < m.nLinks(); ++i) {
        cdm::Index p = m.jointPosInParam(i);
        cdm::Transform dA;
        switch (m.joint(i).type()) {
        case cdm::Joint::Type::Free: {
            dA.translation() = mc.q.template segment<3>(p + 4);
            [[fallthrough]];
        }
        case cdm::Joint::Type::Spherical: {
            dA.rotation() = Eigen::Quaterniond{ mc.q(p), mc.q(p + 1), mc.q(p + 2), mc.q(p + 3) }.toRotationMatrix();
            break;
        }
        default:
            dA = expSE3(m.joint(i).S() * mc.q.segment(p, m.joint(i).dof()));
        }

        auto& jm = mc.jointMotions[static_cast<size_t>(i)];
        jm.transform() = m.T0(i) * dA;
        for (cdm::Index n = 0; n < data.order; ++n) {
            jm.motion()[n] = m.joint(i).S() * mc.dqs.col(n).segment(m.jointPosInDof(i), m.joint(i).dof());
        }

        jm.construct();
    }

    auto v = cdm::MotionVectorX<Order>::Zero(data.order);
    v[1] = cdm::MotionVector(Eigen::Vector3d::Zero(), data.gravity);
    mc.world.set(cdm::Transform::Identity(), v);
}
