#ifndef BICUBIC_HPP
#define BICUBIC_HPP

#include <algorithm>
#include <vector>
#include <cmath>

#include <Eigen/Dense>

namespace Slic3r::Geometry {

namespace BicubicInternal {
/****
 * These kernels have a 4x4 matrix that can be used to interpolate using four
 * equally spaced reference points {f0,f1,f2,f3} and x in the range [0.0,1.0]
 * indicating the interpolation point relative distance between f1 and f2.
 * The interpolated point is calculated with;
 *
 * p(x) = u(x) * A * f
 *
 * Where;
 *
 *  u(x) = RowVector4S{1, x, x^2, x^3}
 *  A = Matrix4S{<4x4 matrix kernel>}
 *  f = Vector4V{f0,f1,f2,f3}
 *
 * Note this means you can also calculate a vector of coefficients for a
 * given x that can be used to interpolate for the same x and different f
 * values. This requires 18+4*i multiplications, which is more efficient for
 * bicubic and tricubic interpolation which requires multiple interpolations
 * for the same z (i=16) and y (i=4) offsets before interpolating x;
 *
 * cint(x) = u(x) * A
 * p(f) = cint(x) * f
 *
 * It is also possible to precalculate weighted points to interpolate for
 * multiple different `x` values between the same `f` reference points, which
 * can be more efficient for iterating to find a crossing point. This
 * requires 16+6*i multiplications, compared to 22*i for interpolating from
 * scratch, so it's cheaper if doing more than 1 interpolation;
 *
 * fint(f) = A * f
 * p(x) = u(x) * fint(f)
 *
 * If we want to bicubic interpolate an (x,y) grid we can use a Matrix4V `f`
 * value with a 4x4 grid of f(y,x) reference points to interpolate. Again the
 * cint(x).translate() and cint(y) values can be pre-computed for
 * interpolating multiple reference points with the same offsets, requiring
 * 36 + 20*i multiplications;
 *
 * p(x,y) = cint(y) * f * cint(x).translate()
 *
 * If we want to bicubic interpolate multiple x,y values between the same
 * f1-f2 reference points, we can pre-calculate a fint(f) value. This
 * requires 128+24i multiplications, compared to 56*i multiplications for
 * interpolating each point from scratch, so it is only cheaper if
 * interpolating more than 4 points;
 *
 *  p(x,y) = cint(y) * f * cint(x).translate()
 *         = (u(y) * A) * f * (u(x) * A).translate()
 *         = u(y) * (A * f * A.translate()) * u(x).translate()
 *
 * fint(f) = A * f * A.translate()
 * p(x,y) = u(y) * fint(f) * u(x).translate()
 *
 * If the points we want to interpolate are actually P(x,y,z) vectors, there
 * are a few ways to do this;
 *
 * * Interpolate each vector dimension independently.
 * * Embed the vectors using ScalarWrapper into the Matrixes.
 * * Use eigen Tensors with multiple dimensions.
 * * Flatten/map the multiple dimensions into matrixes.
 *
 * Because we normally want to interpolate a whole vector, data locality
 * suggests we want the vector dimensions closely packed, which suggests not
 * storing and interpolating each vector dimension independently.
 *
 * For flattening/mapping it all into Matrixes, this involves mapping a
 * logical f(d,x,y) index to a float vector of "d" dimensions at a point
 * (x,y) in a 2D vector field into a matrix(r,c) index. Ideally we want the
 * in-memory layout to reflect the locality order priority of d,x,y. This
 * means the indexing-order and row-vs-col major are important. When
 * flattening down to a 2d Matrix any two adjacent indexes can be merged into
 * a single index. This gives a few options;
 *
 * * row-major and indexed as f(y,x,d). We can flatten this to m(y*xsize+x,d)
 *   or m(y,x*dsize+d). This means roughly col=x,row=y, but relies on
 *   row-major layout and the point vector is a row which are not the eigen
 *   defaults.
 *
 * * col-major with a layout of f(d,x,y) which can be flattened into
 *   m(d,y*xsize+x) or m(r=x*dsize+d,y). This uses the default column-major
 *   mode, and the point vectors are columns, but it does mean x=row,y=col,
 *   which is not normally how it's visualized.
 *
 * * row-major and indexed as f(y,d,x) which can be flattened into
 *   m(y*dmax+d,x) or m(y,d*xsize+x). This doesn't meet the preferred
 *   locality, but has the big advantage that we can collapse the vector
 *   dimensions into the row or column, which simplifys the matrix
 *   operations.
 *
 * The way we want to map it depends on the operation we want to do. To
 * interpolate a vector from an f(y,d,x) field with size f(4,3,4), putting
 * interpolation points in a matrix f(4*d,4) we use;
 *
 * P(x,y) = (f * cint(x).translate()).remap(d,4) * cint(y).translate()
 */

// Linear kernel, to be able to test cubic methods with hat kernels.
template<typename _Scalar>
struct LinearKernel
{
    typedef _Scalar                     Scalar;
    typedef Eigen::Matrix<Scalar, 4, 4> Matrix4S;

    static inline const Matrix4S A{{0, 1, 0, 0}, {0, -1, 1, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
};

// Interpolation kernel aka Catmul-Rom aka Keyes kernel.
template<typename _Scalar>
struct CubicCatmulRomKernel
{
    typedef _Scalar                     Scalar;
    typedef Eigen::Matrix<Scalar, 4, 4> Matrix4S;

    static inline const Matrix4S A = Matrix4S({{0, 2, 0, 0}, {-1, 0, 1, 0}, {2, -5, 4, -1}, {-1, 3, -3, 1}}) / Scalar(2);
};

// B-spline kernel.
template<typename _Scalar>
struct CubicBSplineKernel
{
    typedef _Scalar                     Scalar;
    typedef Eigen::Matrix<Scalar, 4, 4> Matrix4S;

    static inline const Matrix4S A = Matrix4S({{1, 4, 1, 0}, {-3, 0, 3, 0}, {3, -6, 3, 0}, {-1, 3, -3, 1}}) / Scalar(6);
};

} // namespace BicubicInternal

// CubicKernelWrapper implementing interpolation for a Kernel and ValueType.
//
// Note the Kernel is defined with a Scalar type used for the Kernel
// coefficients, multipliers, and offsets, which can be different to the
// ValueType of the points being interpolated between. By default ValueType
// is the same as the Kernel::Scalar, but can be explicitly set to something
// different like ScalarWrapper(Vector3f). The only requirements are that
// ValueType's can be added together and multiplied by Scalar in Matrix
// operations.
template<typename _Kernel, typename _ValueType>
struct CubicKernelWrapper
{
    typedef _Kernel                        Kernel;
    typedef _ValueType                     ValueType;
    typedef typename Kernel::Scalar        Scalar;
    typedef Eigen::Vector<Scalar, 4>       Vector4S;
    typedef Eigen::RowVector<Scalar, 4>    RowVector4S;
    typedef Eigen::Vector<ValueType, 4>    Vector4V;
    typedef Eigen::Matrix<ValueType, 4, 4> Matrix4V;

    // The span number of points that this interpolates over.
    static constexpr Eigen::Index kernel_span = 4;

    // Get a u row vector for x of {1, x, x^2, x^3}.
    static RowVector4S u(const Scalar x)
    {
        assert(0.0 <= x && x <= 1.0);
        RowVector4S v{1, x, x, x};
        v(2) *= x;
        v(3) *= v(2);
        return v;
    }

    // Get interpolation coefficients cint=u*A for a given x.
    static RowVector4S cint(const Scalar x) { return u(x) * Kernel::A; }

    // Get the piecewise interpolation kernel for x (used by Curves.cpp).
    static Scalar kernel(Scalar x)
    {
        x = fabs(x);
        if (x >= Scalar(2))
            return 0.0f;
        if (x <= Scalar(1))
            return u(x) * Kernel::A.col(1);
        assert(x > Scalar(1) && x < Scalar(2));
        return u(x - 1) * Kernel::A.col(0);
    }

    // Interpolate 1D with a 4 point vector and a remainder.
    template<typename Derived>
    static ValueType interpolate4(const Eigen::MatrixBase<Derived>& f, const Scalar x)
    { return cint(x) * f; }

    // Interpolate 1D with 4 points and a remainder.
    static ValueType interpolate4(const ValueType& f0, const ValueType& f1, const ValueType& f2, const ValueType& f3, const Scalar x)
    { return interpolate4(Vector4V(f0, f1, f2, f3), x); }

    // Get a 4x1 vector from a long vector for 1D interpolating at ix.
    template<typename Derived>
    static const Eigen::Ref<const Vector4V> fblock(const Eigen::MatrixBase<Derived>& f, const Eigen::Index ix)
    {
        const Eigen::Index w = f.size();
        // If inside the fully interpolated region, just use segment().
        if (ix > 0 && ix + 2 < w)
            return f.template segment<4>(ix - 1);
        // Overlaps outside the matrix, use clamped indexing to extend outside the matrix.
        auto cx = [&w](const Eigen::Index x) { return std::clamp<Eigen::Index>(x, 0, w - 1); };
        return f({cx(ix - 1), cx(ix), cx(ix + 1), cx(ix + 2)});
    }

    // Get a 4x4 block from a large matrix for D2 interpolating at ix,iy
    template<typename Derived>
    static const Eigen::Ref<const Matrix4V> fblock(const Eigen::MatrixBase<Derived>& f, const Eigen::Index iy, const Eigen::Index ix)
    {
        const Eigen::Index w = f.cols();
        const Eigen::Index h = f.rows();
        // If inside the fully interpolated region, just use block().
        if (ix > 0 && ix + 2 < w && iy > 0 && iy + 2 < h)
            return f.template block<4, 4>(ix - 1, iy - 1);
        // Overlaps outside the matrix, use clamped indexing to extend outside the matrix.
        auto cx = [&w](const Eigen::Index x) { return std::clamp<Eigen::Index>(x, 0, w - 1); };
        auto cy = [&h](const Eigen::Index y) { return std::clamp<Eigen::Index>(y, 0, h - 1); };
        return f({cy(iy - 1), cy(iy), cy(iy + 1), cy(iy + 2)}, {cx(ix - 1), cx(ix), cx(ix + 1), cx(ix + 2)});
    }

    // Interpolate 1D with a vector and an x position.
    template<typename Derived>
    static ValueType interpolate(const Eigen::MatrixBase<Derived>& f, const Scalar x)
    {
        const Eigen::Index ix = std::floor(x);
        const Scalar       rx = x - Scalar(ix);
        return cint(rx) * fblock(f, ix);
    }

    // Interpolate a Matrix 2D scalar field at (x,y). Note x is the column
    // and y is the row, so the Matrix indexing order is m(y,x).
    template<typename Derived>
    static ValueType interpolate(const Eigen::MatrixBase<Derived>& f, const Scalar y, const Scalar x)
    {
        const Eigen::Index ix = std::floor(x);
        const Eigen::Index iy = std::floor(y);
        const Scalar       rx = x - Scalar(ix);
        const Scalar       ry = y - Scalar(iy);
        return cint(ry) * fblock(f, iy, ix) * cint(rx).transpose();
    }

    // Interpolate a 2D scalar field in a matrix at a point p. The point p is
    // any vector of Scalars, but only .x() and .y() are used.
    template<typename Derived, typename Point>
    static ValueType interpolate(const Eigen::MatrixBase<Derived>& f, const Eigen::MatrixBase<Point>& p)
    { return interpolate(f, p.y(), p.x()); }
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
