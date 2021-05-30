#pragma once
#include <Eigen/Core>
#include <cppad/cppad.hpp>
#include <limits>

namespace Eigen {

template <class Base>
struct NumTraits<CppAD::AD<Base>> {
    // type that corresponds to the real part of an AD<Base> value
    typedef CppAD::AD<Base> Real;
    // type for AD<Base> operations that result in non-integer values
    typedef CppAD::AD<Base> NonInteger;
    //  type to use for numeric literals such as "2" or "0.5".
    typedef CppAD::AD<Base> Literal;
    // type for nested value inside an AD<Base> expression tree
    typedef CppAD::AD<Base> Nested;

    enum {
        // does not support complex Base types
        IsComplex             = 0 ,
        // does not support integer Base types
        IsInteger             = 0 ,
        // only support signed Base types
        IsSigned              = 1 ,
        // must initialize an AD<Base> object
        RequireInitialization = 1 ,
        // computational cost of the corresponding operations
        ReadCost              = 1 ,
        AddCost               = 2 ,
        MulCost               = 2
    };

    // machine epsilon with type of real part of x
    // (use assumption that Base is not complex)
    static CppAD::AD<Base> epsilon(void)
    {
        return CppAD::numeric_limits<CppAD::AD<Base>>::epsilon();
    }

    // relaxed version of machine epsilon for comparison of different
    // operations that should result in the same value
    static CppAD::AD<Base> dummy_precision(void)
    {
        return 100. * CppAD::numeric_limits<CppAD::AD<Base>>::epsilon();
    }

    // minimum normalized positive value
    static CppAD::AD<Base> lowest(void)
    {
        return CppAD::numeric_limits< CppAD::AD<Base> >::min();
    }

    // maximum finite value
    static CppAD::AD<Base> highest(void)
    {
        return CppAD::numeric_limits< CppAD::AD<Base> >::max();
    }

    // number of decimal digits that can be represented without change.
    static int digits10(void)
    {
        return CppAD::numeric_limits< CppAD::AD<Base> >::digits10;
    }
};

}

namespace CppAD {

// functions that return references
template <class Base>
const AD<Base>& conj(const AD<Base>& x) { return x; }
template <class Base>
const AD<Base>& real(const AD<Base>& x) { return x; }

// functions that return values (note abs is defined by cppad.hpp)
template <class Base> AD<Base>
imag(const AD<Base>& x) { return CppAD::AD<Base>(0.); }
template <class Base> AD<Base>
abs2(const AD<Base>& x) { return x * x; }

}

template <typename Type>
Eigen::Matrix<Type, 3, 3> exp3x3AD(const Eigen::Matrix<Type, 3, 1>& ax)
{
    Eigen::Matrix<Type, 3, 3> out;
    out.setZero();
    Type angle = CppAD::sqrt(ax.dot(ax));
    if (CppAD::abs(angle) < std::numeric_limits<double>::epsilon()) {
        out(0, 0) = Type(1);
        out(1, 1) = Type(1);
        out(2, 2) = Type(1);
    } else {
        Type x = ax(0) / angle;
        Type y = ax(1) / angle;
        Type z = ax(2) / angle;
        Type sa = CppAD::sin(angle);
        Type ca = CppAD::cos(angle);
        Type t = 1 - ca;

        out(0, 0) = ca + t * x * x;
        out(1, 0) = t * y * x + sa * z;
        out(2, 0) = t * z * x - sa * y;
        out(0, 1) = t * x * y - sa * z;
        out(1, 1) = ca + t * y * y;
        out(2, 1) = t * z * y + sa * x;
        out(0, 2) = t * x * z + sa * y;
        out(1, 2) = t * y * z - sa * x;
        out(2, 2) = ca + t * z * z;
    }

    return out;
}
