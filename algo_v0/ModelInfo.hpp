#pragma once

#include "SimpleHumanModel.hpp"
#include "SimpleManipulator.hpp"
#include "coma/Core"
#include <Eigen/Core>
#include <ostream>
#include <RBDyn/EulerIntegration.h>

struct Model {
    std::string name;
    rbd::MultiBody mb;
    rbd::MultiBodyConfig mbc;
    rbd::MultiBodyGraph mbg;
};

template <size_t Order>
struct ModelInfo {
    static constexpr int ord = Order >= 2 ? Order : 2;

    ModelInfo(int nrStep, double step)
        : model()
        , nt(nrStep)
        , dt(step)
        , t(Eigen::VectorXd::LinSpaced(nt, 0., (nt - 1) * dt))
    {
    }

    void init(const std::string& modelName, size_t nrJoints = 5)
    {
        model.name = modelName;
        if (modelName == "human") {
            auto [a, b, c] = makeHumanBody();
            model.mb = a;
            model.mbc = b;
            model.mbg = c;
        } else if (modelName == "manipulator") {
            auto [a, b, c] = makeManipulator(nrJoints);
            model.mb = a;
            model.mbc = b;
            model.mbg = c;
        }
        const auto& mb = model.mb;
        auto& mbc = model.mbc;

        q.resize(model.mb.nrParams(), nt);
        dqs.resize(ord);
        std::fill(dqs.begin(), dqs.end(), Eigen::MatrixXd(model.mb.nrDof(), nt));

        q.col(0) = rbd::paramToVector(mb, mbc.q);
        // Compute velocity as a * sin(wt + b)
        Eigen::ArrayXd a = Eigen::ArrayXd::Random(mb.nrDof());
        Eigen::ArrayXd b = Eigen::ArrayXd::Random(mb.nrDof());
        Eigen::ArrayXd w = Eigen::ArrayXd::Random(mb.nrDof());
        for (int i = 0; i < nt; ++i) {
            size_t n = ord / 2;
            for (size_t k = 0; k < n; ++k) {
                dqs[2 * k].col(i) = std::pow(-1., k) * a * w.pow(2 * k) * (w * t(i) + b).sin();
                dqs[2 * k + 1].col(i) = std::pow(-1., k) * a * w.pow(2 * k + 1) * (w * t(i) + b).cos();
            }
            if (ord % 2 == 1) {
                dqs[ord - 1].col(i) = std::pow(-1., n) * a * w.pow(2 * n) * (w * t(i) + b).sin();
            }
        }
        mbc.alpha = rbd::vectorToDof(mb, dqs[0].col(0));
        mbc.alphaD = rbd::vectorToDof(mb, dqs[1].col(0));

        for (int i = 0; i < nt - 1; ++i) {
            auto alpha = rbd::vectorToDof(mb, dqs[0].col(i));
            auto alphaD = rbd::vectorToDof(mb, dqs[1].col(i));
            for (int j = 0; j < mb.nrJoints(); ++j)
                rbd::eulerJointIntegration(mb.joint(j).type(), alpha[j], alphaD[j], dt, mbc.q[j]);

            q.col(i + 1) = rbd::paramToVector(mb, mbc.q);
        }
        mbc.q = rbd::vectorToParam(mb, q.col(0));
    }

    Eigen::VectorXd getAleph(int t) const
    {
        const auto& jointPosInDof = model.mb.jointsPosInDof();
        const auto& factors = coma::factorial_factors<double, ord>;
        Eigen::VectorXd v(ord * model.mb.nrDof());
        int curOrderPos = 0;
        for (int i = 0; i < model.mb.nrJoints(); ++i) {
            int dof = model.mb.joint(i).dof();
            for (size_t k = 0; k < ord; ++k) {
                v.segment(curOrderPos + k * dof, dof) = dqs[k].col(t).segment(jointPosInDof[i], dof) / factors[k];
            }
            curOrderPos += ord * dof;
        }

        return v;
    }

    friend std::ostream& operator<<(std::ostream& os, const ModelInfo& info)
    {
        const auto& mb = info.model.mb;
        os << "\nModel description"
           << "\n\tModel name:           " << info.model.name
           << "\n\tNr joints and bodies: " << mb.nrJoints()
           << "\n\tNr joint parameters:  " << mb.nrParams()
           << "\n\tNr DoF:               " << mb.nrDof();
        return os;
    }

    Model model;
    int nt;
    double dt;
    Eigen::VectorXd t;
    Eigen::MatrixXd q;
    std::vector<Eigen::MatrixXd> dqs;
};
