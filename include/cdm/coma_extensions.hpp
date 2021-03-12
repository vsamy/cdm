#pragma once

#include <coma/Core>

namespace coma {

template <typename, int>
class DBX;

template <typename, int>
class LBTM66;

namespace internal {

// typechecks

template <typename _Scalar, int _NVec>
struct traits<DBX<_Scalar, _NVec>> {
    static constexpr int n_vec = _NVec;
    using Scalar = _Scalar;
    using underlying_t = Eigen::MatrixXd;
};

template <typename _Scalar, int _Order>
struct traits<LBTM66<_Scalar, _Order>> {
    static constexpr int order = _Order;
    static constexpr int n_vec = order == Dynamic ? Dynamic : order + 1;
    using Scalar = _Scalar;
    using diagonal_t = Eigen::MatrixXd;
    using off_diagonal_t = Eigen::MatrixXd;
    using mat_t = Eigen::MatrixXd;
    using storage_t = Storage<off_diagonal_t, order>;
};

} // namespace internal

// Block triangular type

template <typename Derived>
class LowerBlockTriangularMatrixTT {
    friend Derived;
    using traits = typename internal::traits<Derived>;
    using Scalar = typename traits::Scalar;
    using diagonal_t = typename traits::diagonal_t;
    using off_diagonal_t = typename traits::off_diagonal_t;
    using mat_t = typename traits::mat_t;
    using storage_t = typename traits::storage_t;

public:
    LowerBlockTriangularMatrixTT() = default;

    const diagonal_t& diagPart() const noexcept { return m_diagonalPart; }
    diagonal_t& diagPart() noexcept { return m_diagonalPart; }
    const storage_t& offDiagonalPart() const noexcept { return m_offDiagonalPart; }
    storage_t& offDiagonalPart() noexcept { return m_offDiagonalPart; }
    const off_diagonal_t& operator[](Index i) const noexcept { return m_offDiagonalPart[i]; }
    off_diagonal_t& operator[](Index i) noexcept { return m_offDiagonalPart[i]; }
    const off_diagonal_t& at(Index i) const { return m_offDiagonalPart.at(i); }
    off_diagonal_t& at(Index i) { return m_offDiagonalPart.at(i); }

    Derived& setOrder(Index order)
    {
        m_offDiagonalPart.resize(order);
        return this->derived();
    }
    Derived& setZero()
    {
        COMA_STATIC_ASSERT((internal::has_setZero<diagonal_t, diagonal_t&()>::value), "diagonal_t has no setZero() function");
        COMA_STATIC_ASSERT((internal::has_setZero<off_diagonal_t, off_diagonal_t&()>::value), "diagonal_t has no setZero() function");
        m_diagonalPart.setZero();
        for (Index i = 0; i < m_offDiagonalPart.size(); ++i)
            m_offDiagonalPart[i].setZero();

        return this->derived();
    }
    Derived& setZero(Index rows, Index cols)
    {
        COMA_STATIC_ASSERT((internal::has_setZero<diagonal_t, diagonal_t&(Index, Index)>::value), "diagonal_t has no setZero() function");
        COMA_STATIC_ASSERT((internal::has_setZero<off_diagonal_t, off_diagonal_t&(Index, Index)>::value), "diagonal_t has no setZero() function");
        m_diagonalPart.setZero(rows, cols);
        for (Index i = 0; i < m_offDiagonalPart.size(); ++i)
            m_offDiagonalPart[i].setZero(rows, cols);

        return this->derived();
    }
    Derived& setIdentity()
    {
        COMA_STATIC_ASSERT((internal::has_setIdentity<diagonal_t, diagonal_t&()>::value), "diagonal_t has no setZero() function");
        COMA_STATIC_ASSERT((internal::has_setZero<off_diagonal_t, off_diagonal_t&()>::value), "diagonal_t has no setZero() function");
        m_diagonalPart.setIdentity();
        for (Index i = 0; i < m_offDiagonalPart.size(); ++i)
            m_offDiagonalPart[i].setZero();

        return this->derived();
    }
    Derived& setIdentity(Index rows, Index cols)
    {
        COMA_STATIC_ASSERT((internal::has_setIdentity<diagonal_t, diagonal_t&(Index, Index)>::value), "diagonal_t has no setZero() function");
        COMA_STATIC_ASSERT((internal::has_setZero<off_diagonal_t, off_diagonal_t&(Index, Index)>::value), "diagonal_t has no setZero() function");
        m_diagonalPart.setIdentity(rows, cols);
        for (Index i = 0; i < m_offDiagonalPart.size(); ++i)
            m_offDiagonalPart[i].setZero(rows, cols);

        return this->derived();
    }

    Derived inverse() const
    {
        // TODO: add assert in all operators
        Derived out;
        out.m_offDiagonalPart.resize(m_offDiagonalPart.size());
        out.m_diagonalPart = m_diagonalPart.inverse();
        for (Index i = 0; i < m_offDiagonalPart.size(); ++i) {
            out.m_offDiagonalPart[i].resizeLike(m_offDiagonalPart[i]);
            out.m_offDiagonalPart[i] = m_offDiagonalPart[i] * out.m_diagonalPart;
            for (Index j = 0; j < i; ++j)
                out.m_offDiagonalPart[i] += m_offDiagonalPart[j] * out.m_offDiagonalPart[i - 1 - j];

            out.m_offDiagonalPart[i] = -(out.m_diagonalPart * out.m_offDiagonalPart[i]);
        }

        return out;
    }

    // operators

    friend Derived operator*(const Derived& lhs, const Derived& rhs)
    {
        // TODO: add assert in all operators
        Derived out;
        out.m_diagonalPart = lhs.m_diagonalPart * rhs.m_diagonalPart;
        for (Index i = 0; i < lhs.m_offDiagonalPart.size(); ++i) {
            out[i] = lhs.m_diagonalPart * rhs.m_offDiagonalPart[i] + rhs.m_offDiagonalPart[i] * lhs.m_diagonalPart;
            for (Index j = 0; j < i; ++j)
                out.m_offDiagonalPart[i] += lhs.m_offDiagonalPart[j] * rhs.m_offDiagonalPart[i - 1 - j];
        }

        return out;
    }

    friend Derived operator-(Derived lhs, const Derived& rhs)
    {
        lhs -= rhs;
        return rhs;
    }

    Derived& operator-=(const Derived& rhs)
    {
        m_diagonalPart -= rhs.m_diagonalPart;
        for (Index i = 0; i < m_offDiagonalPart.size(); ++i) {
            m_offDiagonalPart[i] -= rhs.m_offDiagonalPart[i];
        }
        return this->derived();
    }

    friend Derived operator+(Derived lhs, const Derived& rhs)
    {
        lhs += rhs;
        return rhs;
    }
    Derived& operator+=(const Derived& rhs)
    {
        m_diagonalPart += rhs.m_diagonalPart;
        for (Index i = 0; i < m_offDiagonalPart.size(); ++i) {
            m_offDiagonalPart[i] += rhs.m_offDiagonalPart[i];
        }
        return this->derived();
    }

    template <int NVec>
    Derived& operator+=(const DiInertia<Scalar, NVec>& rhs)
    {
        COMA_ASSERT(m_diagonalPart.rows() == rhs.block().rows() && m_diagonalPart.cols() == rhs.block().cols(), "Diagonal matrix size mismatches");
        COMA_ASSERT(([&rhs, this]() {
            bool res = true;
            for (Index i = 0; i < m_offDiagonalPart.size(); ++i) {
                res = res && m_offDiagonalPart[i].rows() == rhs.block().rows() && m_offDiagonalPart[i].cols() == rhs.block().cols();
            }
            return res;
        }()), // Check sub-matrices size
            "Diagonal matrix size mismatches");
        m_diagonalPart += rhs.block().matrix();
        return this->derived();
    }
    template <int NVec>
    Derived& operator*=(const DiMotionSubspace<Scalar, NVec>& rhs)
    {
        COMA_ASSERT(m_diagonalPart.rows() == rhs.block().rows() && m_diagonalPart.cols() == rhs.block().cols(), "Diagonal matrix size mismatches");
        COMA_ASSERT(([&rhs, this]() {
            bool res = true;
            for (Index i = 0; i < m_offDiagonalPart.size(); ++i) {
                res = res && m_offDiagonalPart[i].rows() == rhs.block().rows() && m_offDiagonalPart[i].cols() == rhs.block().cols();
            }
            return res;
        }()), // Check sub-matrices size
            "Diagonal matrix size mismatches");
        m_diagonalPart *= rhs.block().matrix();
        return this->derived();
    }

    friend bool operator==(const Derived& lhs, const Derived& rhs) noexcept
    {
        return lhs.m_diagonalPart == rhs.m_diagonalPart && lhs.m_offDiagonalPart == rhs.m_offDiagonalPart;
    }
    friend bool operator!=(const Derived& lhs, const Derived& rhs) noexcept
    {
        return !(lhs == rhs);
    }
    bool isApprox(const Derived& rhs, Scalar prec = dummy_precision<Scalar>()) const noexcept
    {
        return m_diagonalPart.isApprox(rhs.m_diagonalPart, prec) && m_offDiagonalPart.isApprox(rhs.m_offDiagonalPart, prec);
    }

private:
    inline Derived& derived() { return static_cast<Derived&>(*this); }
    inline const Derived& derived() const { return static_cast<const Derived&>(*this); }

private:
    diagonal_t m_diagonalPart;
    storage_t m_offDiagonalPart;
};

// TODO: test operator on this
template <typename Scalar, int NVec>
class DBX : public DiBlockT<DBX<Scalar, NVec>> {
public:
    DBX() = default;
};

template <typename Scalar, int Order>
class LBTM66 : public LowerBlockTriangularMatrixTT<LBTM66<Scalar, Order>> {
public:
    BTMXX() = default;
};

// operators

template <typename Scalar, int Order, int NVec>
LBTM66<Scalar, Order> operator+(LBTM66<Scalar, Order> lhs, const DiInertia<Scalar, NVec>& rhs)
{
    lhs += rhs;
    return lhs;
}
template <typename Scalar, int Order, int NVec>
LBTM66<Scalar, Order> operator+(LBTM66<Scalar, Order> lhs, const DiMotionSubspace<Scalar, NVec>& rhs)
{
    lhs += rhs;
    return lhs;
}
template <typename Scalar, int Order, int NVec>
LBTM66<Scalar, Order> operator+(const DiMotionSubspace<Scalar, NVec>& lhs, LBTM66<Scalar, Order> rhs)
{
    rhs += lhs;
    return rhs;
}

template <typename Scalar, int Order, int NVec>
LBTM66<Scalar, Order> operator*(LBTM66<Scalar, Order> lhs, const DiMotionSubspace<Scalar, NVec>& rhs)
{
    lhs *= rhs;
    return lhs;
}

template <typename Scalar, int Order, int NVec>
BTMXX<Scalar, Order> operator*(const DiMotionSubspace<Scalar, NVec>& lhs, BTMXX<Scalar, Order> rhs)
{
    // TODO: add assert in all operators
    rhs.diagPart() = lhs.block().matrix() * rhs.diagPart();
    return rhs;
}

template <typename Scalar, int Order, typename Derived>
Eigen::VectorXd operator*(const BTMXX<Scalar, Order>& lhs, const Eigen::MatrixBase<Derived>& rhs)
{
    // TODO: add assert in all operators
    Eigen::VectorXd out{ rhs.size() };
    for (Index i = 0; i < Order; ++i) {
        out.segment<6>(6 * i) = lhs.diagPart() * rhs.segment<6>(6 * i);
        for (Index j = 0; j < i; ++j)
            out.segment<6>(6 * i) += lhs[i - j - 1] * rhs.segment<6>(6 * j);
    }

    return out;
}

template <typename Scalar, int Order, typename Derived>
LBTM66<Scalar, Order> DualMul(const CMTM<Scalar, 6, Order>& lhs, const LBTM66<Scalar, Order>& rhs)
{
    LBTM66<Scalar, Order> out;
    out.diagPart() = lhs.transformation().dualMatrix() * rhs.diagPart();
    for (Index i = 0; i < Order; ++i) {
        out[i] = lhs.transformation().dualMatrix() * rhs[i] + lhs[i] * rhs.diagPart();
        for (Index j = 0; j < i; ++j) {
            out[i] += lhs[i] * rhs[i - j];
        }
    }

    return out;
}

template <typename Scalar, int Order>
ForceVectorX<Scalar, Order> DualMul(const BTMXX<Scalar, Order>& lhs, const ForceVectorX<Scalar, Order>& rhs)
{
    return lhs * rhs;
}

} // namespace coma
