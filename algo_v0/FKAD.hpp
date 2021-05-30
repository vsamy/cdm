// #include "EigenAD.hpp"
#include "utilityAD.hpp"
#include "EigenAD.hpp"
#include <cppad/cppad.hpp> // the CppAD package http://www.coin-or.org/CppAD/
#include <vector>

template <typename Type>
using MatrixAD = Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>;
template <typename Type>
using VectorAD = Eigen::Matrix<Type, Eigen::Dynamic, 1>;

template <typename Type>
struct JointAD {
    std::vector<Type> S;
    std::vector<Type> R0;
    std::vector<Type> p0;
    std::vector<Type> R;
    std::vector<Type> p;
    std::vector<Type> w;
    std::vector<Type> v;
    std::vector<Type> dw;
    std::vector<Type> dv;
};

template <typename Type>
struct EigenJointAD {
    Eigen::Matrix<Type, 6, Eigen::Dynamic> S;
    Eigen::Matrix<Type, 3, 3> R0;
    Eigen::Matrix<Type, 3, 1> p0;
    Eigen::Matrix<Type, 3, 3> R;
    Eigen::Matrix<Type, 3, 1> p;
    Eigen::Matrix<Type, 3, 1> w;
    Eigen::Matrix<Type, 3, 1> v;
    Eigen::Matrix<Type, 3, 1> dw;
    Eigen::Matrix<Type, 3, 1> dv;
};

template <typename Model, typename Type>
std::vector<Type> FKVecAD(const std::vector<Type>& Q, const Model& model, std::vector<JointAD<Type>>& jointAD)
{
    const auto& mb = model.mb;
    const auto& mbc = model.mbc;
    const auto& Xt = mb.transforms();
    const std::vector<rbd::Joint>& joints = mb.joints();
    const std::vector<int>& pred = mb.predecessors();
    const std::vector<int>& succ = mb.successors();
    const std::vector<int>& jointPosInDof = mb.jointsPosInDof();
    jointAD.resize(mb.nrJoints());

    for (std::size_t i = 0; i < joints.size(); ++i) {
        jointAD[i].S = Eigen2AD<Type>(mbc.motionSubspace[i]);
        jointAD[i].R0 = Eigen2AD<Type>(Xt[i].rotation().inverse());
        jointAD[i].p0 = Eigen2AD<Type>(Xt[i].translation());

        std::vector<Type> q(Q.begin() + jointPosInDof[i], Q.begin() + jointPosInDof[i] + joints[i].dof());
        auto xj = ProductAD(jointAD[i].S, q, 6);
        std::vector<Type> dq(Q.begin() + mb.nrDof() + jointPosInDof[i], Q.begin() + mb.nrDof() + jointPosInDof[i] + joints[i].dof());
        auto dxj = ProductAD(jointAD[i].S, dq, 6);
        std::vector<Type> ax(xj.begin(), xj.begin() + 3);
        auto dR = exp3x3AD(ax);
        ax.assign(xj.begin() + 3, xj.end());

        auto Rr = ProductAD(jointAD[i].R0, dR, 3, 3);
        auto pr = SumAD(ProductAD(jointAD[i].R0, ax, 3), jointAD[i].p0);
        std::vector<Type> wr(dxj.begin(), dxj.begin() + 3);
        std::vector<Type> vr(dxj.begin() + 3, dxj.end());

        int parent = pred[i];
        std::vector<Type> Rp, pp;
        std::vector<Type> wp, vp;
        if (parent == -1) {
            Rp = IdMatAD<Type>(3);
            pp.assign(3, Type(0));
            wp.assign(3, Type(0));
            vp.assign(3, Type(0));
        } else {
            Rp = jointAD[parent].R;
            pp = jointAD[parent].p;
            wp = jointAD[parent].w;
            vp = jointAD[parent].v;
        }

        jointAD[i].R = ProductAD(Rp, Rr, 3, 3);
        jointAD[i].p = SumAD(ProductAD(Rp, pr, 3), pp);

        auto Rt = TransposeAD(Rr, 3);
        auto vTmp = SumAD(CrossAD(wp, pr), vp);
        jointAD[i].v = SumAD(ProductAD(Rt, vTmp, 3), vr);
        jointAD[i].w = SumAD(ProductAD(Rt, wp, 3), wr);
    }

    std::vector<Type> out = jointAD[mb.nrJoints() - 1].w;
    out.insert(out.end(), jointAD[mb.nrJoints() - 1].v.begin(), jointAD[mb.nrJoints() - 1].v.end());
    return out;
}

template <typename Model, typename Type>
std::vector<Type> FKAccAD(const std::vector<Type>& Q, const Model& model, std::vector<JointAD<Type>>& jointAD)
{
    const auto& mb = model.mb;
    const auto& mbc = model.mbc;
    const auto& Xt = mb.transforms();
    const std::vector<rbd::Joint>& joints = mb.joints();
    const std::vector<int>& pred = mb.predecessors();
    const std::vector<int>& succ = mb.successors();
    const std::vector<int>& jointPosInDof = mb.jointsPosInDof();
    jointAD.resize(mb.nrJoints());

    for (std::size_t i = 0; i < joints.size(); ++i) {
        jointAD[i].S = Eigen2AD<Type>(mbc.motionSubspace[i]);
        jointAD[i].R0 = Eigen2AD<Type>(Xt[i].rotation().inverse());
        jointAD[i].p0 = Eigen2AD<Type>(Xt[i].translation());

        std::vector<Type> q(Q.begin() + jointPosInDof[i], Q.begin() + jointPosInDof[i] + joints[i].dof());
        auto xj = ProductAD(jointAD[i].S, q, 6);
        std::vector<Type> dq(Q.begin() + mb.nrDof() + jointPosInDof[i], Q.begin() + mb.nrDof() + jointPosInDof[i] + joints[i].dof());
        auto dxj = ProductAD(jointAD[i].S, dq, 6);
        std::vector<Type> ddq(Q.begin() + 2 * mb.nrDof() + jointPosInDof[i], Q.begin() + 2 * mb.nrDof() + jointPosInDof[i] + joints[i].dof());
        auto ddxj = ProductAD(jointAD[i].S, ddq, 6);
        std::vector<Type> ax(xj.begin(), xj.begin() + 3);
        auto dR = exp3x3AD(ax);
        ax.assign(xj.begin() + 3, xj.end());

        auto Rr = ProductAD(jointAD[i].R0, dR, 3, 3);
        auto pr = SumAD(ProductAD(jointAD[i].R0, ax, 3), jointAD[i].p0);
        std::vector<Type> wr(dxj.begin(), dxj.begin() + 3);
        std::vector<Type> vr(dxj.begin() + 3, dxj.end());
        std::vector<Type> dwr(ddxj.begin(), ddxj.begin() + 3);
        std::vector<Type> dvr(ddxj.begin() + 3, ddxj.end());

        int parent = pred[i];
        std::vector<Type> Rp, pp;
        std::vector<Type> wp, vp;
        std::vector<Type> dwp, dvp;
        if (parent == -1) {
            Rp = IdMatAD<Type>(3);
            pp.assign(3, Type(0));
            wp.assign(3, Type(0));
            vp.assign(3, Type(0));
            dwp.assign(3, Type(0));
            dvp.assign(3, Type(0));
        } else {
            Rp = jointAD[parent].R;
            pp = jointAD[parent].p;
            wp = jointAD[parent].w;
            vp = jointAD[parent].v;
            dwp = jointAD[parent].dw;
            dvp = jointAD[parent].dv;
        }

        jointAD[i].R = ProductAD(Rp, Rr, 3, 3);
        jointAD[i].p = SumAD(ProductAD(Rp, pr, 3), pp);

        auto Rt = TransposeAD(Rr, 3);
        auto wTmp = ProductAD(Rt, wp, 3);
        auto vTmp = ProductAD(Rt, SumAD(CrossAD(wp, pr), vp), 3);
        jointAD[i].v = SumAD(vTmp, vr);
        jointAD[i].w = SumAD(wTmp, wr);

        auto tmp1 = CrossAD(wTmp, wr);
        auto tmp2 = ProductAD(Rt, dwp, 3);
        jointAD[i].dw = SumAD(SumAD(tmp1, tmp2), dwr);
        tmp1 = SumAD(CrossAD(vTmp, wr), CrossAD(wTmp, vr));
        tmp2 = ProductAD(Rt, SumAD(dvp, CrossAD(dwp, pr)), 3);
        jointAD[i].dv = SumAD(SumAD(tmp1, tmp2), dvr);
    }

    std::vector<Type> out = jointAD[mb.nrJoints() - 1].dw;
    out.insert(out.end(), jointAD[mb.nrJoints() - 1].dv.begin(), jointAD[mb.nrJoints() - 1].dv.end());
    return out;
}

template <typename Model, typename Type>
Eigen::Matrix<Type, 6, 1> EigenFKAD(const Eigen::Matrix<Type, Eigen::Dynamic, 1>& Q, const Model& model, std::vector<EigenJointAD<Type>>& jointAD)
{
    const auto& mb = model.mb;
    const auto& mbc = model.mbc;
    const auto& Xt = mb.transforms();
    const std::vector<rbd::Joint>& joints = mb.joints();
    const std::vector<int>& pred = mb.predecessors();
    const std::vector<int>& succ = mb.successors();
    const std::vector<int>& jointPosInDof = mb.jointsPosInDof();
    jointAD.resize(mb.nrJoints());

    for (std::size_t i = 0; i < joints.size(); ++i) {
        jointAD[i].S = mbc.motionSubspace[i].template cast<Type>();
        jointAD[i].R0 = Xt[i].rotation().inverse().template cast<Type>();
        jointAD[i].p0 = Xt[i].translation().template cast<Type>();

        Eigen::Matrix<Type, Eigen::Dynamic, 1> q = Q.segment(jointPosInDof[i], joints[i].dof());
        Eigen::Matrix<Type, 6, 1> xj = jointAD[i].S * q;
        Eigen::Matrix<Type, Eigen::Dynamic, 1> dq = Q.segment(mb.nrDof() + jointPosInDof[i], joints[i].dof());
        Eigen::Matrix<Type, 6, 1> dxj = jointAD[i].S * dq;
        // Eigen::Matrix<Type, Eigen::Dynamic, 1> ddq = Q.segment(2 * mb.nrDof() + jointPosInDof[i], joints[i].dof());
        // Eigen::Matrix<Type, 6, 1> ddxj = jointAD[i].S * ddq;
        Eigen::Matrix<Type, 3, 3> dR = exp3x3AD(Eigen::Matrix<Type, 3, 1>{ xj.template head<3>() });

        Eigen::Matrix<Type, 3, 3> Rr = jointAD[i].R0 * dR;
        Eigen::Matrix<Type, 3, 1> pr = jointAD[i].R0 * xj.template tail<3>() + jointAD[i].p0;
        Eigen::Matrix<Type, 3, 1> wr(dxj.data());
        Eigen::Matrix<Type, 3, 1> vr(dxj.data() + 3);
        // Eigen::Matrix<Type, 3, 1> dwr(ddxj.data());
        // Eigen::Matrix<Type, 3, 1> dvr(ddxj.data() + 3);

        int parent = pred[i];
        Eigen::Matrix<Type, 3, 3> Rp;
        Eigen::Matrix<Type, 3, 1> wp, vp, pp, dwp, dvp;
        if (parent != -1) {
            Rp = jointAD[parent].R;
            pp = jointAD[parent].p;
            wp = jointAD[parent].w;
            vp = jointAD[parent].v;
            // dwp = jointAD[parent].dw;
            // dvp = jointAD[parent].dv;
        } else {
            Rp.setIdentity();
            pp.setZero();
            wp.setZero();
            vp.setZero();
            // dwp.setZero();
            // dvp.setZero();
        }

        jointAD[i].R = Rp * Rr;
        jointAD[i].p = Rp * pr + pp;

        Eigen::Matrix<Type, 3, 1> wTmp = Rr.transpose() * wp;
        Eigen::Matrix<Type, 3, 1> vTmp = Rr.transpose() * (wp.cross(pr) + vp);
        jointAD[i].v = vTmp + vr;
        jointAD[i].w = wTmp + wr;

        // auto tmp1 = wTmp.cross(wr);
        // auto tmp2 = Rr.transpose() * dwp;
        // jointAD[i].dw = tmp1 + tmp2 + dwr;
        // tmp1 = vTmp.cross(wr) + wTmp.cross(vr);
        // tmp2 = Rr.transpose() * (dvp + dwp.cross(pr));
        // jointAD[i].dv = tmp1 + tmp2 + dvr;
    }

    Eigen::Matrix<Type, 6, 1> out;
    out << jointAD[mb.nrJoints() - 1].w, jointAD[mb.nrJoints() - 1].v;
    // out << jointAD[mb.nrJoints() - 1].dw, jointAD[mb.nrJoints() - 1].dv;
    return out;
}
