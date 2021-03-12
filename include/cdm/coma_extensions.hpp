#pragma once

#include <coma/Core>

namespace coma {

namespace internal {

// typechecks

template <typename _Scalar, int _NVec>
struct traits<DBX<_Scalar, _NVec>> {
    static constexpr int n_vec = _NVec;
    using Scalar = _Scalar;
    using underlying_t = Eigen::MatrixXd;
};

} // namespace internal

template <typename Scalar, int NVec>
class DBX : public DiBlockT<DBX<Scalar, NVec>> {
public:
    DBX() = default;
};

} // namespace coma
