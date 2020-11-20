#pragma once

#include <coma/Core>

namespace cdm {

using MotionVector = coma::MotionVector<double>;
using ForceVector = coma::ForceVector<double>;
using Transform = coma::Transform<double>;
using Cross = coma::Cross<double>;
using Inertia = coma::SpatialInertia<double>;
using MotionSubspace = coma::MotionSubspace<double>;

template <int NVec>
using CrossN = coma::CrossN<double, NVec>;
template <int Order>
using CMTM = coma::CMTM<double, 6, Order>;
template <int NVec>
using MotionVectorX = coma::MotionVectorX<double, NVec>;
template <int NVec>
using ForceVectorX = coma::ForceVectorX<double, NVec>;
template <int NVec>
using DiInertia = coma::DiInertia<double, NVec>;
template <int NVec>
using DiMotionSubspace = coma::DiMotionSubspace<double, NVec>;

} // namespace cdm
