// #include "EigenAD.hpp"
#include "utilityAD.hpp"
#include "EigenAD.hpp"
#include <cppad/cppad.hpp> // the CppAD package http://www.coin-or.org/CppAD/
#include <vector>

template <typename Type>
using MatrixAD = Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>;
template <typename Type>
using VectorAD = Eigen::Matrix<Type, Eigen::Dynamic, 1>;

template <typename Model, typename Type>
std::vector<Type> IDAD(const std::vector<Type>& motion, const Model& model)
{
    const auto& mb = model.mb;
    const auto& mbc = model.mbc;
    const std::vector<rbd::Joint>& joints = mb.joints();
    const std::vector<rbd::Body>& bodies = mb.bodies();
    const std::vector<int>& pred = mb.predecessors();
    int nJoints = mb.nrJoints();
    std::vector<std::vector<Type>> f(nJoints);
    std::vector<Type> tau;

    sva::MotionVecd a_0(Eigen::Vector3d::Zero(), mbc.gravity);
    std::vector<sva::ForceVecd> f_(mb.nrJoints());
    std::vector<sva::MotionVecd> bodyAccB(mb.nrJoints());
    std::vector<std::vector<double>> jointTorque(mb.nrJoints());

    // First pass
    auto it = motion.begin();
    std::vector<std::vector<Type>> vb(nJoints);
    std::vector<std::vector<Type>> ab(nJoints);
    for (int i = 0; i < nJoints; ++i) {
        const sva::PTransformd & X_p_i = mbc.parentToSon[i];

        const sva::MotionVecd & vj_i = mbc.jointVelocity[i];
        sva::MotionVecd ai_tan = joints[i].tanAccel(mbc.alphaD[i]);

        const sva::MotionVecd & vb_i = mbc.bodyVelB[i];

        std::vector<Type> psi(it, it + joints[i].dof());
        std::vector<Type> dpsi(it + joints[i].dof(), it + 2 * joints[i].dof());
        it += 2 * joints[i].dof();

        auto S = Eigen2AD<Type>(mbc.motionSubspace[i]);
        auto X = Eigen2AD<Type>(X_p_i.matrix());
        auto vi = ProductAD(S, psi, 6, 1);
        auto ai = ProductAD(S, dpsi, 6, 1);
        if(pred[i] != -1) {
            vb[i] = SumAD(ProductAD(X, vb[pred[i]], 6), vi);
            auto tmp1 = ProductAD(X, ab[pred[i]], 6);
            auto tmp2 = Cross6AD(vb[i], vi);
            auto rbdTmp1 = X_p_i * bodyAccB[pred[i]];
            auto rbdTmp2 = vb_i.cross(vj_i);
            ab[i] = SumAD(SumAD(tmp1, tmp2), ai);
            bodyAccB[i] = rbdTmp1 + rbdTmp2 + ai_tan;
        } else {
            vb[i] = vi;
            auto tmp1 = ProductAD(X, Eigen2AD<Type>(a_0.vector()), 6);
            auto tmp2 = Cross6AD(vb[i], vi);
            auto rbdTmp1 = X_p_i * a_0;
            auto rbdTmp2 = vb_i.cross(vj_i);
            ab[i] = SumAD(SumAD(tmp1, tmp2), ai);
            bodyAccB[i] = rbdTmp1 + rbdTmp2 + ai_tan;
        }

        std::vector<Type> I = Eigen2AD<Type>(bodies[i].inertia().matrix());
        auto tmp1 = ProductAD(I, ab[i], 6);
        auto tmp2 = Cross6DAD(vb[i], ProductAD(I, vb[i], 6));
        // auto tmp3 = external force
        auto rbdTmp1 = bodies[i].inertia() * bodyAccB[i];
        auto rbdTmp2 = vb_i.crossDual(bodies[i].inertia() * vb_i);
        f[i] = SumAD(tmp1, tmp2);
        f_[i] = rbdTmp1 + rbdTmp2 - mbc.bodyPosW[i].dualMul(mbc.force[i]);
    }

    for (int i = nJoints - 1; i >= 0; --i) {
        auto ST = Eigen2AD<Type>(mbc.motionSubspace[i].transpose());
        auto tmp = ProductAD(ST, f[i], mbc.motionSubspace[i].cols(), 1);
        tau.insert(tau.begin(), tmp.begin(), tmp.end());

        jointTorque[i].resize(mbc.motionSubspace[i].cols());
        for (int j = 0; j < joints[i].dof(); ++j) {
            jointTorque[i][j] = mbc.motionSubspace[i].col(j).transpose() * f_[i].vector();
        }

        if(pred[i] != -1)
        {
            const sva::PTransformd & X_p_i = mbc.parentToSon[i];
            auto X = Eigen2AD<Type>(X_p_i.matrix().transpose());
            f[pred[i]] = SumAD(f[pred[i]], ProductAD(X, f[i], 6));
            f_[pred[i]] = f_[pred[i]] + X_p_i.transMul(f_[i]);
        }
    }

    return tau;
}
