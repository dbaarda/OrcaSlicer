#ifndef BICUBIC_HPP
#define BICUBIC_HPP

#include <algorithm>
#include <vector>
#include <cmath>

#include <Eigen/Dense>

/***
 * We define some extensions to Eigen here to make it easier to define
 * Eigen::Vector and Eigen::RowVector types, and so we can put any type of
 * Eigen::Matrix inside another Eigen::Matrix or Eigen::Array. This means we
 * can interpolate not just scalar fields, but also vector fields.
 */
namespace Eigen {

// Define Eigen::Vector and Eigen::RowVector templates to make it easier to define these.
template<typename _Scalar, int _Rows, int _Options = AutoAlign | ColMajor, int _MaxRows = _Rows>
using Vector = Matrix<_Scalar, _Rows, 1, _Options, _MaxRows, 1>;
template<typename _Scalar, int _Cols, int _Options = AutoAlign | RowMajor, int _MaxCols = _Cols>
using RowVector = Matrix<_Scalar, 1, _Cols, _Options, 1, _MaxCols>;

// Define Eigen::NumTraits<Matrix<T,Y,X> so we can put Matrix and Vector
// types inside other Matrix or Array types. This copies the existing
// Eigen::NumTraits<Array<...>> definition.
template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
struct NumTraits<Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>>
{
    typedef Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>          MatrixType;
    typedef typename NumTraits<_Scalar>::Real                                    RealScalar;
    typedef Matrix<RealScalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>       Real;
    typedef typename NumTraits<_Scalar>::NonInteger                              NonIntegerScalar;
    typedef Matrix<NonIntegerScalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> NonInteger;
    typedef MatrixType&                                                          Nested;
    typedef typename NumTraits<_Scalar>::Literal                                 Literal;

    enum {
        IsComplex             = NumTraits<_Scalar>::IsComplex,
        IsInteger             = NumTraits<_Scalar>::IsInteger,
        IsSigned              = NumTraits<_Scalar>::IsSigned,
        RequireInitialization = 1,
        ReadCost = MatrixType::SizeAtCompileTime == Dynamic ? HugeCost : MatrixType::SizeAtCompileTime * NumTraits<_Scalar>::ReadCost,
        AddCost  = MatrixType::SizeAtCompileTime == Dynamic ? HugeCost : MatrixType::SizeAtCompileTime * NumTraits<_Scalar>::AddCost,
        MulCost  = MatrixType::SizeAtCompileTime == Dynamic ? HugeCost : MatrixType::SizeAtCompileTime * NumTraits<_Scalar>::MulCost
    };

    // This is the Eigen v5.0.0 version.
    // EIGEN_DEVICE_FUNC constexpr static RealScalar epsilon() { return NumTraits<RealScalar>::epsilon(); }
    // EIGEN_DEVICE_FUNC constexpr static RealScalar dummy_precision() { return NumTraits<RealScalar>::dummy_precision(); }
    // constexpr static int digits10() { return NumTraits<Scalar>::digits10(); }
    // constexpr static int max_digits10() { return NumTraits<Scalar>::max_digits10(); }

    EIGEN_DEVICE_FUNC static inline RealScalar epsilon() { return NumTraits<RealScalar>::epsilon(); }
    EIGEN_DEVICE_FUNC static inline RealScalar dummy_precision() { return NumTraits<RealScalar>::dummy_precision(); }
    static inline int                          digits10() { return NumTraits<_Scalar>::digits10(); }
};

} // namespace Eigen

namespace Slic3r::Geometry {

namespace BicubicInternal {
/****
 * These kernels have a 4x4 matrix that can be used to interpolate using four
 * equally spaced reference points f0,f1,f2,f3 and x in the range [0.0,1.0]
 * indicating the interpolation point relative distance between f1 and f2.
 * The interpolated point is calculated with;
 *
 * p = u * a * f
 *
 * Where;
 *
 *  u = RowVector{1, x, x^2, x^3}
 *  a = <4x4 matrix kernel>
 *  f = Vector{f0,f1,f2,f3}
 *
 * Note this means you can also calculate a vector of coefficients for a
 * given x that can be used to interpolate for the same x and different f
 * values. This is more efficient for bicubic and tricubic interpolation,
 * which requires multiple interpolations for the same z (16x) and y (4x)
 * offsets before interpolating x;
 *
 * cint = u * a;
 * p = cint * f;
 */

// Linear kernel, to be able to test cubic methods with hat kernels.
template<typename _Scalar> struct LinearKernel
{
    typedef _Scalar                     Scalar;
    typedef Eigen::Matrix<Scalar, 4, 4> Matrix4S;

    static inline const Matrix4S a = (Matrix4S() << 0, 1, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0).finished();
};

// Interpolation kernel aka Catmul-Rom aka Keyes kernel.
template<typename _Scalar> struct CubicCatmulRomKernel
{
    typedef _Scalar                     Scalar;
    typedef Eigen::Matrix<Scalar, 4, 4> Matrix4S;

    static inline const Matrix4S a = (Matrix4S() << 0, 2, 0, 0, -1, 0, 1, 0, 2, -5, 4, -1, -1, 3, -3, 1).finished() / Scalar(2);
};

// B-spline kernel
template<typename _Scalar> struct CubicBSplineKernel
{
    typedef _Scalar                     Scalar;
    typedef Eigen::Matrix<Scalar, 4, 4> Matrix4S;

    static inline const Matrix4S a = (Matrix4S() << 1, 4, 1, 0, -3, 0, 3, 0, 3, -6, 3, 0, -1, 3, -3, 1).finished() / Scalar(6);
};

} // namespace BicubicInternal

// CubicKernelWrapper implementing the interpolation for a given Kernel and ValueType.
//
// Note the Kernel is defined with a Scalar type used for the Kernel
// coefficients, multipliers, and offsets, which can be different to the
// ValueType of the points being interpolated between. By default ValueType
// is the same as the Kernel::Scalar, but can be explicitly set to something
// different like a vector. The only requirements are that ValueType can be
// multiplied by Scalar and added together.
template<typename _Kernel, typename _ValueType> struct CubicKernelWrapper
{
    typedef _Kernel                        Kernel;
    typedef _ValueType                     ValueType;
    typedef typename Kernel::Scalar        Scalar;
    typedef Eigen::RowVector<Scalar, 4>    RowVector4S;
    typedef Eigen::Vector<ValueType, 4>    Vector4V;
    typedef Eigen::Matrix<ValueType, 4, 4> Matrix4V;

    // The span number of points that this interpolates over.
    static constexpr size_t kernel_span = 4;

    // Get a u row vector for x of {1, x, x^2, x^3}.
    static const RowVector4S u(const Scalar x)
    {
        assert(0 <= x && x <= 1.0);
        RowVector4S v{1, x, x, x};
        v(2) *= x;
        v(3) *= v(2);
        return v;
    }

    // Get interpolation coefficients cint=u*a for a given x.
    static const RowVector4S cint(const Scalar x)
    {
        assert(0 <= x && x <= 1.0);
        return u(x) * Kernel::a;
    }

    // Get the piecewise interpolation kernel for x (used by Curves.cpp).
    static const Scalar kernel(Scalar x)
    {
        x = fabs(x);
        if (x >= Scalar(2))
            return 0.0f;
        if (x <= Scalar(1))
            return u(x) * Kernel::a.col(1);
        assert(x > Scalar(1) && x < Scalar(2));
        x -= Scalar(1);
        return u(x) * Kernel::a.col(0);
    }

    // Interpolate 1D with a 4 point vector and a remainder.
    static ValueType interpolate4(const Eigen::Ref<const Vector4V>& f, const Scalar x) { return u(x) * Kernel::a * f; }

    // Interpolate 1D with 4 points and a remainder.
    static ValueType interpolate4(ValueType f0, ValueType f1, ValueType f2, ValueType f3, Scalar x)
    { return interpolate4(Vector4V(f0, f1, f2, f3), x); }

    // Get a 4x1 vector from a vector for 1D interpolating with.
    template<typename Derived> static const Vector4V fblock(const Eigen::MatrixBase<Derived>& f, const Eigen::Index ix)
    {
        const Eigen::Index w = f.size();
        // If inside the fully interpolated region, just use segment()
        if (ix > 0 && ix + 2 < w)
            return f.template segment<4>(ix - 1);
        // Overlaps outside the matrix, use a NullaryExpr to extend outside the matrix.
        return Vector4V::NullaryExpr([&f, ix, w](const Eigen::Index x) { return f(std::clamp<Eigen::Index>(x + ix - 1, 0, w - 1)); })
            .template segment<4>(0);
    }

    // Get a 4x4 block from a matrix for D2 interpolating with.
    template<typename Derived>
    static const Matrix4V fblock(const Eigen::MatrixBase<Derived>& f, const Eigen::Index iy, const Eigen::Index ix)
    {
        const Eigen::Index w = f.cols();
        const Eigen::Index h = f.rows();
        // If inside the fully interpolated region, just use block()
        if (ix > 0 && ix + 2 < w && iy > 0 && iy + 2 < h)
            return f.template block<4, 4>(ix - 1, iy - 1);
        // Overlaps outside the matrix, use a NullaryExpr to extend outside the matrix.
        return Matrix4V::NullaryExpr([&f, ix, iy, w, h](const Eigen::Index y, const Eigen::Index x) {
                   return f(std::clamp<Eigen::Index>(y + iy - 1, 0, h - 1), std::clamp<Eigen::Index>(x + ix - 1, 0, w - 1));
               })
            .template block<4, 4>(0, 0);
    }

    // Interpolate 1D with a vector and an x position.
    template<typename Derived> static ValueType interpolate(const Eigen::MatrixBase<Derived>& f, const Scalar x)
    {
        const Eigen::Index ix = std::floor(x);
        const Scalar       rx = x - Scalar(ix);
        // return u(rx) * Kernel::a * fblock(f, ix);
        return interpolate4(fblock(F, ix), rx);
    }

    // Interpolate a 2D scalar field in a matrix. Note x is the column and y is
    // the row, so the Matrix indexing order is m(y,x). The point p is any
    // vector of Scalars, but only .x() and .y() are used.
    template<typename Derived, typename Point>
    static ValueType interpolate(const Eigen::MatrixBase<Derived>& f, const Eigen::MatrixBase<Point>& p)
    {
        const Eigen::Index ix = std::floor(p.x());
        const Eigen::Index iy = std::floor(p.y());
        const Scalar       rx = p.x() - Scalar(ix);
        const Scalar       ry = p.y() - Scalar(iy);
        // return u(rx) * (u(ry) * Kernel::a * fblock(f, iy, ix).transpose());
        return interpolate4((cint(ry) * fblock(F, iy, ix)).transpose(), rx);
    }
};

// Linear splines
template<typename Scalar, typename ValueType = Scalar>
using LinearKernel = CubicKernelWrapper<BicubicInternal::LinearKernel<Scalar>, ValueType>;

// Catmul-Rom splines
template<typename Scalar, typename ValueType = Scalar>
using CubicCatmulRomKernel = CubicKernelWrapper<BicubicInternal::CubicCatmulRomKernel<Scalar>, ValueType>;

// Cubic B-splines
template<typename Scalar, typename ValueType = Scalar>
using CubicBSplineKernel = CubicKernelWrapper<BicubicInternal::CubicBSplineKernel<Scalar>, ValueType>;

} // namespace Slic3r::Geometry

#endif /* BICUBIC_HPP */
