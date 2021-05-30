#pragma once

#include <Eigen/Core>
#include <cppad/cppad.hpp> // the CppAD package http://www.coin-or.org/CppAD/
#include <vector>

template <typename Type, typename Derived>
std::vector<Type> Eigen2AD(const Eigen::MatrixBase<Derived>& M)
{
    std::vector<Type> out;
    out.resize(M.size());

    Eigen::Index p = 0;
    for (Eigen::Index j = 0; j < M.cols(); ++j)
        for (Eigen::Index i = 0; i < M.rows(); ++i)
            out[p++] = M(i, j);

    return out;
}

template <typename Type>
std::vector<Type> SumAD(const std::vector<Type>& M1, const std::vector<Type>& M2)
{
    std::vector<Type> out(M1.size());

    for (size_t i = 0; i < M1.size(); ++i)
        out[i] = M1[i] + M2[i];

    return out;
}

template <typename Type>
std::vector<Type> SubAD(const std::vector<Type>& M1, const std::vector<Type>& M2)
{
    std::vector<Type> out(M1.size());

    for (size_t i = 0; i < M1.size(); ++i)
        out[i] = M1[i] - M2[i];

    return out;
}

template <typename Type>
std::vector<Type> ProductAD(const std::vector<Type>& M1, const std::vector<Type>& M2, size_t rowsOut, size_t colsOut = 1)
{
    std::vector<Type> out(rowsOut * colsOut, Type(0));

    size_t m = rowsOut;
    size_t n = M1.size() / rowsOut;
    size_t p = colsOut;
    for (size_t j = 0; j < p; ++j)
        for (size_t i = 0; i < m; ++i)
            for (size_t k = 0; k < n; ++k)
                out[i + j * p] += M1[i + k * m] * M2[k + j * p];

    return out;
}

template <typename Type>
Type DotProductAD(const std::vector<Type>& v1, const std::vector<Type>& v2)
{
    Type sum;

    sum = 0.;
    for (size_t i = 0; i < size_t(v1.size()); ++i)
        sum += v1[i] * v2[i];

    return sum;
}

template <typename Type>
std::vector<Type> CrossAD(const std::vector<Type>& v1, const std::vector<Type>& v2)
{
    std::vector<Type> out(3);
    out[0] = v1[1] * v2[2] - v1[2] * v2[1];
    out[1] = v1[2] * v2[0] - v1[0] * v2[2];
    out[2] = v1[0] * v2[1] - v1[1] * v2[0];
    return out;
}

template <typename Type>
std::vector<Type> Cross6AD(const std::vector<Type>& v1, const std::vector<Type>& v2)
{
    std::vector<Type> out(6);
    out[0] = v1[1] * v2[2] - v1[2] * v2[1];
    out[1] = v1[2] * v2[0] - v1[0] * v2[2];
    out[2] = v1[0] * v2[1] - v1[1] * v2[0];
    out[3] = v1[4] * v2[2] - v1[5] * v2[1] + v1[1] * v2[5] - v1[2] * v2[4];
    out[4] = v1[5] * v2[0] - v1[3] * v2[2] + v1[2] * v2[3] - v1[0] * v2[5];
    out[5] = v1[3] * v2[1] - v1[4] * v2[0] + v1[0] * v2[4] - v1[1] * v2[3];

    return out;
}

template <typename Type>
std::vector<Type> Cross6DAD(const std::vector<Type>& v1, const std::vector<Type>& v2)
{
    std::vector<Type> out(6);
    out[0] = v1[1] * v2[2] - v1[2] * v2[1] + v1[4] * v2[5] - v1[5] * v2[4];
    out[1] = v1[2] * v2[0] - v1[0] * v2[2] + v1[5] * v2[3] - v1[3] * v2[5];
    out[2] = v1[0] * v2[1] - v1[1] * v2[0] + v1[3] * v2[4] - v1[4] * v2[3];
    out[3] = v1[1] * v2[5] - v1[2] * v2[4];
    out[4] = v1[2] * v2[3] - v1[0] * v2[5];
    out[5] = v1[0] * v2[4] - v1[1] * v2[3];

    return out;
}

template <typename Type>
std::vector<Type> CrossAD(const std::vector<Type>& v)
{
    std::vector<Type> out(9);
    out[0] = Type(0);
    out[1] = v[2];
    out[2] = -v[1];
    out[3] = -v[2];
    out[4] = Type(0);
    out[5] = v[0];
    out[6] = v[1];
    out[7] = -v[0];
    out[8] = Type(0);
    return out;
}

template <typename Type>
std::vector<Type> TransposeAD(const std::vector<Type>& M, int rows)
{
    std::vector<Type> out(M.size());
    size_t cols = M.size() / rows;
    size_t p = 0;
    for (size_t j = 0; j < cols; ++j)
        for (size_t i = 0; i < rows; ++i)
            out[p++] = M[j + cols * i];

    return out;
}

template <typename Type>
std::vector<Type> IdMatAD(size_t size)
{
    std::vector<Type> out(size * size, Type(0));
    for (size_t i = 0; i < size; ++i)
        out[i + i * size] = Type(1);

    return out;
}

template <typename Type>
std::vector<Type> exp3x3AD(const std::vector<Type>& ax)
{
    std::vector<Type> out(9, Type(0));
    Type angle = CppAD::sqrt(DotProductAD(ax, ax));
    if (CppAD::abs(angle) < std::numeric_limits<double>::epsilon()) {
        out[0] = 1;
        out[4] = 1;
        out[8] = 1;
    } else {
        Type x = ax[0] / angle;
        Type y = ax[1] / angle;
        Type z = ax[2] / angle;
        Type sa = CppAD::sin(angle);
        Type ca = CppAD::cos(angle);

        out[0] = ca + (1 - ca) * x * x;
        out[1] = (1 - ca) * y * x + sa * z;
        out[2] = (1 - ca) * z * x - sa * y;
        out[3] = (1 - ca) * x * y - sa * z;
        out[4] = ca + (1 - ca) * y * y;
        out[5] = (1 - ca) * z * y + sa * x;
        out[6] = (1 - ca) * x * z + sa * y;
        out[7] = (1 - ca) * y * z - sa * x;
        out[8] = ca + (1 - ca) * z * z;
    }

    return out;
}

template <template <class> typename ADType, typename Type>
std::vector<Type> CastOutAD(const std::vector<ADType<Type>>& in)
{
    size_t s = in.size();
    std::vector<Type> out(s);
    for (size_t i = 0; i < s; ++i)
        out[i] = CppAD::Value(in[i]);

    return out;
}