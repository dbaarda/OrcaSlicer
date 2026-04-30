#ifndef libslic3r_Geometry_Fields_hpp_
#define libslic3r_Geometry_Fields_hpp_

#include <algorithm>
#include <vector>
#include <cmath>

#include <Eigen/Dense>
#include "Bicubic.hpp"
#include "ScalarWrapper.hpp"

namespace Slic3r::Geometry {

/************************************************
 * 2D and 3D Fields of type T.
 *
 * This implements 2D and 3D fields of type T as grids of values. The type
 * T can be any numeric type or vector.
 *
 * Note that spacial data is normaly processed first in the x-dimension in
 * the inner-most loop, then the y-dimension, then the z-dimension in the
 * outer-most loop. This means we normally want the data packed in that order
 * to optimize data locality. Eigen Matrixes default to column-major storage,
 * so we exlicitly use row-major storage to ensure better data locality when
 * iterating through x.
 *
 * Note since Eigen Matrix indexing uses m(row, col) and 2D points use x for
 * the col and y for the row, a 2D point p(x,y) is indexed in a 2D field
 * using f2(y,x). A 3D field is implemented as a Vector of 2D fields, and is
 * indexed with f3(z)(y,x). This is why we reverse the more typical x,y,z
 * index order everywhere.
 *
 * Other x,y,z packing and indexing is possible, just use the coordinates
 * consistently in the desired packing order whenever instantiating
 * templates, calling functions, or indexing the fields.
 *
 */

using namespace Eigen;

// Strip const and only const off a variable.
template<typename T>
inline T& unconst(const T& v)
{ return const_cast<T&>(v); }
template<typename T>
inline T* unconst(const T* v)
{ return const_cast<T*>(v); }
template<typename T>
inline T& unconst(T& v)
{ return v; }

// Add const and only const to a variable.
template<typename T>
inline const T& enconst(T& v)
{ return const_cast<const T&>(v); }
template<typename T>
inline const T* enconst(T* v)
{ return const_cast<const T*>(v); }
template<typename T>
inline const T& enconst(const T& v)
{ return v; }

// Helper class for a list of index values with compile-time values, some of
// which can be Dynamic and set on initialization.
template<Index... Vals>
struct Indexes
{
    static constexpr auto size() { return sizeof...(Vals); }
    using ValuesType                                 = std::array<Index, size()>;
    static constexpr ValuesType ValuesAtCompileTime  = {Vals...};
    static constexpr Index      ProductAtCompileTime = getProduct(ValuesAtCompileTime);
    const ValuesType            values;

    // Init from a list of values which can be shorter than size(). When the
    // list is shorter they set the last values, and the earlier values are
    // initialized from ValuesAtCompileTime. The final values must match any
    // non-Dynamic ValuesAtCompileTime, and no values can be Dynamic. This
    // means you can have static sizes for the minor indexes and only need to
    // set the dynamic major indexes.
    constexpr Indexes(std::initializer_list<Index> vals) : values{getValues(vals)} {}
    template<typename... Args>
    constexpr Indexes(Args... args) : Indexes({static_cast<Index>(args)...})
    {}

    auto           begin() { return values.begin(); }
    auto           end() { return values.end(); }
    constexpr auto operator[](size_t i) { return values[i]; }
    constexpr auto operator()(size_t i) { return values[i]; }
    // get the "total size" product of all the values.
    constexpr auto product() const { return getProduct(values); }

private:
    constexpr ValuesType getValues(std::initializer_list<Index> vals) const
    {
        assert(vals.size() <= size());
        ValuesType values;
        size_t     ndef = ValuesAtCompileTime.size() - vals.size();
        std::copy_n(ValuesAtCompileTime.begin(), ndef, values.begin());
        std::copy_n(vals.begin(), vals.size(), values.begin() + ndef);
        for (size_t i = 0; i < size(); i++) {
            assert(values[i] > 0 && (values[i] == ValuesAtCompileTime[i] || ValuesAtCompileTime[i] == Dynamic));
        }
        return values;
    }
    // Get the product of all values, returning Dynamic if any are Dynamic.
    constexpr Index getProduct(const ValuesType& vals) const
    {
        Index size = 1;
        for (auto v : vals) {
            if (v == Dynamic)
                return Dynamic;
            size *= v;
        }
        return size;
    }
};

// get an Indexes list of strides from an indes list of sizes.
template<typename Accumulator, Index Stride, Index... Sizes>
struct getstrides;

template<Index... Strides, Index Stride>
struct getstrides<Indexes<Strides...>, Stride>
{
    using type = Indexes<Strides...>;
};

template<Index... Strides, Index Stride, Index Next, Index... Sizes>
struct getstrides<Indexes<Strides...>, Stride, Next, Sizes...>
{
    using type = getstrides<Indexes<Strides..., Stride>, internal::size_at_compile_time(Stride, Next), Sizes...>::type;
};

template<Index... Sizes>
using getstrides_t = getstrides<Indexes<>, 1, Sizes...>::type;

// A helper for getting strides for N dimensions of given sizes.
// This uses ColMajor (little-endian) indexing, so the first index has the smallest stride.
// TODO: Hmmm... this doesn't support slices, windows etc. maybe there's a better way?
template<Index... Dims>
struct Strides
{
    using SizesType   = Indexes<Dims...>;
    using StridesType = getstrides_t<Dims...>;
    const SizesType         sizes;
    const StridesType       strides;
    static constexpr size_t Rank                 = SizesType::size();
    static constexpr auto&  SizesAtCompileTime   = SizesType::ValuesAtCompileTime;
    static constexpr auto&  StridesAtCompileTime = StridesType::ValuesAtCompileTime;

    constexpr Strides(std::initializer_list<Index> values) : sizes{values}, strides{_getStrides(sizes)} {}
    template<typename... Args>
    constexpr Strides(Args... args) : Strides({static_cast<Index>(args)...})
    {}

    // Size() with no args gives the total size of the linear array.
    constexpr Index Size() const { return sizes[Rank - 1] * strides[Rank - 1]; };
    // Stride() with no args gives the stride for overall linear array.
    constexpr Index Stride() const { return strides[0]; };
    // Size(dim) gives the size of the corresponding index.
    constexpr Index Size(size_t dim) const
    {
        assert(dim < Rank);
        return sizes[dim];
    };
    // Stride(dim) gives the stride for the corresponding index.
    constexpr Index Stride(size_t dim) const
    {
        assert(dim < Rank);
        return strides[dim];
    };

    // Get the linear offset index for a sequence of dimension indexes.
    //
    // This can take a shortened list of dimensions and missing indexes are
    // equal to zero. Because indexing is ColMajor (little-endian) the first
    // indexes are dropped first, so for Strides(X,Y,Z) using Index(y,z) is the
    // same as Index(0,y,z).
    constexpr Index Offset(std::initializer_list<Index> indexes) const
    {
        assert(indexes.size() <= Rank);
        Index  idx = 0;
        size_t i   = Rank - indexes.size();
        ;
        for (auto d : indexes) {
            idx += d * strides[i++];
        }
        return idx;
    }
    template<typename... Args>
    constexpr Index Offset(Args... args) const
    { return Offset({static_cast<Index>(args)...}); }

    // This gets the dimension index values for a given offset.
    constexpr auto Indexes(Index i)
    {
        assert(i < Size());
        std::array<Index, Rank> indexes;
        for (auto s : sizes) {
            indexes[i++] = i % s;
            i /= s;
        }
        return indexes;
    }

private:
    // Calculate the strides from the sizes. Note a size of Dynamic in the list
    // of sizes means all more major index strides are also Dynamic.
    static constexpr auto _getStrides(const SizesType& sizes)
    {
        std::array<Index, sizes.size()> strides;
        Index                           stride = 1;
        size_t                          i      = 0;
        for (auto s : sizes) {
            strides[i++] = stride;
            stride       = (stride == Dynamic || s == Dynamic) ? Dynamic : stride * s;
        }
        return strides;
    }
};

// An p(x,y,z) Vector of type T;
template<typename T>
using Vec3 = Eigen::Vector<T, 3>;

template<typename T>
using SVec3 = ScalarWrapper<Vec3<T>>;

// A 2D Field of type T as a matrix indexed with (x,y)
template<typename T, int X = Dynamic, int Y = Dynamic>
using Field2D = Matrix<T, X, Y>;

template<typename T, int X = Dynamic, int Y = Dynamic>
using SField2D = ScalarWrapper<Field2D<T, X, Y>>;

template<typename T, int X, int Y>
auto getX(Field2D<T, X, Y>& f, const int x)
{ return f.row(x); }

template<typename T, int X, int Y>
auto getY(Field2D<T, X, Y>& f, const int y)
{ return f.col(y); }

// A 3D Field of type T as a matrix indexed with (y*X + x, z).
template<typename T, int X = Dynamic, int Y = Dynamic, int Z = Dynamic>
using _Field3D = Matrix<T, internal::size_at_compile_time(X, Y), Z>;

// A 3D Field of type T.
//
// Note this is implemented as a subclass of _Field3D because we want to add
// some extra Field3D specific methods.
template<typename T, int X = Dynamic, int Y = Dynamic, int Z = Dynamic>
class Field3D
{
public:
    _Field3D<T, X, Y, Z>                    m_data;
    internal::variable_if_dynamic<Index, X> m_xsize;

    Field3D(const Index x = X, const Index y = Y, const Index z = Z) : m_data(x * y, z), m_xsize(x) {}

    template<typename OtherDerived>
    Field3D(const MatrixBase<OtherDerived>& other, Index xsize = X) : m_data(other), m_xsize(xsize)
    { assert(0 < m_xsize && m_xsize <= m_data.rows() && (m_data.rows() % m_xsize) == 0); }

    Field3D& operator=(const Field3D& other)
    {
        m_xsize.setValue(other.m_xsize);
        m_data = other.m_data;
        return *this;
    }

    template<typename OtherDerived>
    Field3D& operator=(const MatrixBase<OtherDerived>& other)
    {
        assert((m_xsize <= other.rows()) && ((other.rows() % m_xsize) == 0));
        m_data = other;
        return *this;
    }
    /*
        Field3D& operator*=(const T& v)
        {
            for (Index i = 0; i < this->size(); i++)
                (*this)(i) *= v;
            return *this;
        }
    */
    inline Index xsize() const { return m_xsize; }
    inline Index ysize() const { return m_data.rows() / m_xsize; }
    inline Index zsize() const { return m_data.cols(); }

    // auto sum() { return m_data.colwise().sum().reshape(m_xsize, AutoSize); }

    Field3D& setRandom()
    {
        m_data.setRandom();
        return (*this);
    }

    Field3D& setConstant(const T& val)
    {
        m_data.setConstant(val);
        return (*this);
    }

    Field3D& setZero() { return setConstant(T(0)); }

    template<typename NewType>
    const auto cast() const
    {
        Field3D<NewType, X, Y, Z> ans(xsize(), ysize(), zsize());
        ans.m_data = m_data.template cast<NewType>();
        return ans;
    }
    template<>
    const auto cast<T>() const
    { return *this; }

    inline const auto getX(Index x) const
    {
        using StrideType = Eigen::Stride<internal::size_at_compile_time(X, Y), X>;
        assert(0 <= x && x < xsize());
        return Map<Field2D<T, Y, Z>, 0, StrideType>(unconst(m_data.data()) + x, ysize(), zsize(), StrideType(m_data.rows(), xsize()));
    }

    inline const auto getY(const Index y) const
    {
        using StrideType = OuterStride<internal::size_at_compile_time(X, Y)>;
        assert(0 <= y && y < ysize());
        return Map<Field2D<T, X, Z>, 0, StrideType>(unconst(m_data.data()) + y * xsize(), xsize(), zsize(), StrideType(m_data.rows()));
    }

    inline const auto getZ(const Index z) const
    {
        assert(0 <= z && z < zsize());
        return Map<Field2D<T, X, Y>>(unconst(m_data.data()) + z * m_data.rows(), xsize(), ysize());
        // return m_data(z).reshape(xsize(), AutoSize);
    }

    inline auto getX(const Index x) { return unconst(enconst(*this).getX(x)); }
    inline auto getY(const Index y) { return unconst(enconst(*this).getY(y)); }
    inline auto getZ(const Index z) { return unconst(enconst(*this).getZ(z)); }

    inline auto        operator()(const Index z) { return getZ(z); }
    inline auto        operator()(const Index y, const Index z) { return getZ(z).col(y); }
    inline auto&       operator()(const Index x, const Index y, const Index z) { return getZ(z)(x, y); }
    inline const auto  operator()(const Index z) const { return getZ(z); }
    inline const auto  operator()(const Index y, const Index z) const { return getZ(z).col(y); }
    inline const auto& operator()(const Index x, const Index y, const Index z) const { return getZ(z)(x, y); }
};

// A 2D Vector Field of Vec3 of type T.
template<typename T, int X = Dynamic, int Y = Dynamic>
using Vec3Field2D = Field2D<SVec3<T>, X, Y>;

// A 3D Vector Field of Vec3 of type T.
template<typename T, int X = Dynamic, int Y = Dynamic, int Z = Dynamic>
using Vec3Field3D = Field3D<SVec3<T>, X, Y, Z>;

/*********************************************
 * Interpolate Fields.
 *
 * These do linear or cubic interpolation to get field values at arbitrary x,
 * p(x,y), or p(x,y,z) points within the field grid. The x,y,z values should
 * be floats with values in the range of [0, <size>) where <size> is f.zsize()
 * for z, f.rows() for y, or f.cols() for x. It will interpolate a value
 * between the enclosing grid points.
 */
// Old Eigen versions don't have this defined.
// using RowVector2d = Matrix<double, 1, 2>;

// A simple linear interpolation "kernel" with a similar API to Bicubic.hpp.
namespace LerpKernel {
const RowVector2d cint(const double x) { return RowVector2d(1 - x, x); }

// Get a 2x1 vector from a vector for lerp.
template<typename Derived>
const auto fblock(const EigenBase<Derived>& F, const Index x)
{
    assert(0 <= x && x < F.size());
    return F.template segment<2>(x);
}

// Get a 2x2 block from a matrix for lerp.
template<typename Derived>
const auto fblock(const EigenBase<Derived>& F, const Index x, const Index y)
{
    assert(0 <= x && x < F.rows());
    assert(0 <= y && y < F.cols());
    return F.template block<2, 2>(x, y);
}

} // namespace LerpKernel

// linearly interpolate within a Field2D to get the value at p(x,y).
template<typename P, typename T, int X, int Y>
const T lerp(const Field2D<T, X, Y>& f, const P& p)
{
    const Index  iy = std::floor(p.y());
    const Index  ix = std::floor(p.x());
    const double ry = p.y() - iy;
    const double rx = p.x() - ix;
    return LerpKernel::cint(rx) * LerpKernel::fblock(f, ix, iy) * LerpKernel::cint(ry).transpose();
}

// linearly interpolate within a Field3D to get the value at p(x,y,z).
template<typename P, typename T, int X, int Y, int Z>
const T lerp(const Field3D<T, X, Y, Z>& f, const P& p)
{
    const Index  iz = std::floor(p.z());
    const Index  iy = std::floor(p.y());
    const Index  ix = std::floor(p.x());
    const double rz = p.z() - iz;
    const double ry = p.y() - iy;
    const double rx = p.x() - ix;
    assert(0 <= iz && iz < f.zsize());
    assert(0 <= iy && iy < f.ysize());
    assert(0 <= ix && ix < f.xsize());
    Vector2d cy = LerpKernel::cint(ry);
    Vector2d cx = LerpKernel::cint(rx).transpose();
    double   z0 = cy * LerpKernel::fblock(f(iz), ix, iy) * cx;
    double   z1 = cy * LerpKernel::fblock(f(iz + 1), ix, iy) * cx;
    return LerpKernel::cint(rz) * Vector2d(z0, z1);
}

// cubic interpolate within a Field2D to get the value at p(x,y).
template<typename P, typename T, int X, int Y>
const T cubic(const Field2D<T, X, Y>& f, const P& p)
{
    using Scalar = typename P::Scalar;
    using Kernel = CubicCatmulRomKernel<Scalar, T>;
    return Kernel::interpolate(f, p);
}

// cubic interpolate within a Field3D to get the value at p(x,y,z).
template<typename P, typename T, int X, int Y, int Z>
const T cubic(const Field3D<T, X, Y, Z>& f, const P& p)
{
    using Scalar    = typename P::Scalar;
    using Kernel    = CubicCatmulRomKernel<Scalar, T>;
    const Index  iz = std::floor(p.z());
    const Index  iy = std::floor(p.y());
    const Index  ix = std::floor(p.x());
    const Scalar rz = p.z() - iz;
    const Scalar ry = p.y() - iy;
    const Scalar rx = p.x() - ix;
    assert(0 <= iz && iz < f.zsize());
    assert(0 <= iy && iy < f.ysize());
    assert(0 <= ix && ix < f.xsize());
    const typename Kernel::Vector4S    cy = Kernel::cint(ry);
    const typename Kernel::RowVector4S cx = Kernel::cint(rx).transpose();
    typename Kernel::Vector4V          zf;
    for (int z = 0; z < 4; z++)
        zf(z) = cx * Kernel::fblock(f(std::clamp<Index>(iz + z - 1, 0, f.zsize() - 1)), ix, iy) * cy;
    return Kernel::interpolate4(zf, rz);
}

// cubic interpolate a Field3D into a Field2D at a given z.
template<typename S, typename T, int X, int Y, int Z>
const Field2D<T, X, Y> cubic_z(const Field3D<T, X, Y, Z>& f, const S z)
{
    using Kernel = CubicCatmulRomKernel<S, ScalarWrapper<Eigen::Vector<T, internal::size_at_compile_time(X, Y)>>>;
    return Kernel::interpolate(AsScalarCols(f.m_data), z).derived().reshaped(f.m_xsize, AutoSize);
}

/*********************************************
 * GausianSmooth (almost) Fields.
 *
 * These do multiple iterations of width=3 (3x3 for 2D or 3x3x3 for 3D)
 * average smoothing, which approaches Gausian smoothing as the iterations
 * increase. After 5 to 6 iterations it is close enough for most
 * applications. See the following for details;
 *
 * https: *www.peterkovesi.com/papers/FastGaussianSmoothing.pdf
 *
 * For n iterations of width 3 averaging the GausianSmooting stddev is
 * roughly;
 *
 * stddev = sqrt((n*3^2 - n)/12) = sqrt(n) * 0.816; or stddev=2.0 for n=6.
 *
 * You need at least 5 or 6 iterations before the approximation gets accurate
 * enough if you want to use the derivative. For 6 iterations the original
 * values are spread 68% within radius 2 and 95% within radius 4 over the
 * adjacent cells.
 */

// Do a single iteration of summing 3x3 adjacent values of Field2D.
template<typename T, int X, int Y>
void sum33(Ref<Field2D<T, X, Y>>& f)
{
    const Index      x = f.rows();
    const Index      y = f.cols();
    Field2D<T, X, Y> buf(x, y);
    // TODO: this could be done using only 3 temporary rows as a circular
    // buffer and cycling through them. Maybe this would be cheaper?
    // First fill buf with the sums of 3 adjacent columns.
    // Set outer-edge rows with the edge value extended.
    buf.row(0)     = 2.0 * f.row(0) + f.row(1);
    buf.row(x - 1) = 2.0 * f.row(x - 1) + f.row(x - 2);
    // Set inner block to sum of 3 row values.
    buf.middleRows(1, x - 2) = f.topRows(x - 2) + f.middleRows(1, x - 2) + f.bottomRows(x - 2);
    // next fill f with the sums of 3 cols from buf.
    // Set outer-edge cols with the edge value extended
    f.col(0)     = 2.0 * buf.col(0) + buf.col(1);
    f.col(y - 1) = 2.0 * buf.col(y - 1) + buf.col(y - 2);
    // Set inner block to sum of 3 col values.
    f.middleCols(1, y - 2) = buf.leftCols(y - 2) + buf.middleCols(1, y - 2) + buf.rightCols(y - 2);
}

// GausianSmooth (almost) a Field2D with multiple iterations of 3x3 averaging.
template<typename D>
void smooth(MatrixBase<D>& f, const int n = 6)
// template<typename F, typename CleanF=std::remove_reference_t<F>, typename T=CleanF::Scalar, int X=CleanF::RowsAtCompileTime, int
// Y=CleanF::ColsAtCompileTime, typename=std::enable_if_t<std::is_same_v<CleanF, Field2D<T, X, Y> >, int> > void smooth(F&& f, const int n = 6)
{
    double m = 1.0;
    for (int i = 0; i < n; i++) {
        sum33<D::Scalar, D::RowsAtCompileTime, D::ColsAtCompileTime>(f);
        m *= 9.0;
    }
    f *= D::Scalar(1.0 / m);
}
/*
template<typename T, int X, int Y>
void smooth(Ref<Field2D<T,X,Y>>&& f, const int n = 6) {
  smooth<T,X,Y>(f, n);
}
*/
template<typename D>
void smooth(D&& f, const int n = 6)
{ smooth(f, n); }

// Do a single iteration of summing 3 adjacent values of a vector.
// This version uses block-sums.
template<typename T, int X = Dynamic>
void _sum3_1(Ref<Eigen::Vector<T, X>>& v)
{
    const Index x       = v.size();
    T           b0      = 2 * v(0) + v(1);
    T           bn      = v(x - 2) + 2 * v(x - 1);
    v.segment(1, x - 2) = (v.head(x - 2) + v.segment(1, x - 2) + v.tail(x - 2)).eval();
    v(0)                = b0;
    v(x - 1)            = bn;
}

// This version uses a rolling buffer-sum.
template<typename D>
void _sum3_2(MatrixBase<D> const& v)
{
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(D);
    const Index x = v.size();
    // Initialize buf with edge value 0 repeated.
    // Used to hold 3 values i-1, i, i+1 for summing, with value i at index i%3.
    Eigen::Vector<typename D::Scalar, 3> buf(3);
    buf(0) = v(0);
    buf(2) = v(0);
    // Go through values copying the i+1 value into buf and put the sum into v.
    for (Index i = 0; i < x - 1; i++) {
        buf((i + 1) % 3)                 = v(i + 1);
        const_cast<MatrixBase<D>&>(v)(i) = buf.sum();
    }
    // Set buf with edge value x-1 value extended and calc last value.
    buf(x % 3)                           = v(x - 1);
    const_cast<MatrixBase<D>&>(v)(x - 1) = buf.sum();
}

// This version uses pair-sums.
template<typename D, typename = std::enable_if_t<std::remove_reference_t<D>::IsVectorAtCompileTime, int>>
void _sum3_3(D&& v)
{
    using CleanD = std::remove_reference_t<D>;
    //EIGEN_STATIC_ASSERT_VECTOR_ONLY(D);
    const Index x = v.size();
    // If v contains ScalarWrappers of dynamic sized vectors, we can't just create an instance because we don't know the sizes of those
    // vectors. So instead we have to copy it.
    // Eigen::Vector<typename CleanD::Scalar, CleanD::SizeAtCompileTime> buf(x);
    using VectorType=typename internal::plain_matrix_type<CleanD>::type;
    VectorType buf = v;
    // This works by summing input pairs into the even indexed buffer values,
    // then calculating the odd indexed entries as the previous sum plus the
    // next value, and then incrementing the even indexed sums by the
    // previous value. This reuse of the shared-partial sum uses 1.5 adds per
    // output value, instead of 2.
    // buf(seq(fix<0>, indexing::last-fix<1>, fix<2>)) = v(seq(fix<0>, indexing::last-fix<1>, fix<2>)) + v(seq(fix<1>, indexing::last, fix<2>));
    buf(seq(fix<0>, indexing::last - fix<1>, fix<2>)) += buf(seq(fix<1>, indexing::last, fix<2>));
    buf(seq(fix<1>, indexing::last - fix<1>, fix<2>)) = buf(seq(fix<0>, indexing::last - fix<2>, fix<2>)) +
                                                        v(seq(fix<2>, indexing::last, fix<2>));
    buf(indexing::last)                               = (x % 2 ? buf(indexing::last - fix<1>) : v(indexing::last)) + v(indexing::last);
    buf(seq(fix<2>, indexing::last, fix<2>)) += v(seq(fix<1>, indexing::last - fix<1>, fix<2>));
    buf(0) += v(0);
    v = buf;
}

// GausianSmooth (almost) a vector with multiple iterations of 3 value averaging.
template<typename D, typename = std::enable_if_t<std::remove_reference_t<D>::IsVectorAtCompileTime, int>>
void _smoothF1(D&& v, const int n = 6)
{
    //using CleanD = std::remove_reference_t<D>;
    //EIGEN_STATIC_ASSERT_VECTOR_ONLY(D);
    double m = 1.0;
    for (int i = 0; i < n; i++) {
        _sum3_3(v);
        m *= 3.0;
    }
    v = v / m;
}

// GausianSmooth (almost) a Field2D with multiple iterations of 3x3 averaging.
template<typename T, int X = Dynamic, int Y = Dynamic>
void _smoothF2(Field2D<T, X, Y>& f, const int n = 6)
{
    _smoothF1(AsScalarCols(f), n);
    _smoothF1(AsScalarRows(f), n);
    /*
    for (auto r : f.rowwise()) {
        _smooth(r, n);
    }
    for (auto c : f.colwise()) {
        _smooth(c, n);
    }
    */
}

template<typename T, int X = Dynamic, int Y = Dynamic>
void _smoothF2(Field2D<T, X, Y>&& f, const int n = 6)
{ _smoothF2<T, X, Y>(f, n); }

// GausianSmooth (almost) a Field3D with multiple iterations of 3x3x3 averaging.
template<typename T, int X, int Y, int Z>
void _smoothF3(Field3D<T, X, Y, Z>& f, const int n = 6)
{
    const Index zlen = f.zsize();
    // First smooth each layer n times.
    for (Index z = 0; z < zlen; z++)
        _smoothF2<T, X, Y>(f(z), n);
    // Now smooth each z-vector n times.
    for (int x = 0; x < f.xsize(); x++) {
        for (int y = 0; y < f.ysize(); y++) {
            _smoothF1(AsScalarCols(f.m_data), n);
        }
    }
}

// GausianSmooth (almost) a Field3D with multiple iterations of 3x3x3 averaging.
template<typename T, int X, int Y, int Z>
void smooth(Field3D<T, X, Y, Z>& f, const int n = 6)
{
    double      m    = 1.0;
    const Index zlen = f.zsize();
    // First smooth each layer n times.
    // for (Index z = 0; z < zlen; z++)
    //    smooth(f(z), n);
    // Now sum 3 adjacent layers n times.
    _smoothF3(f, n);
    /*
    for (int i = 0; i < n; i++) {
        sum3<T, X, Y, Z>(f);
        m *= 3.0;
    }
    // Now divide the sums by 3^n to get the average.
    f = f * T(1 / m);
     */
}

/*********************************************
 * Integrate Fields.
 *
 * These do trapezoidal integration of 2D and 3D fields in the x, y, z or xyz
 * (all three into a 3D vector) directions.
 *
 * Optional dz,dy,dx arguments are the grid-delta distances in each direction
 * for the integration. They normally default to 1.0, but for zyx integration
 * the defaults are dz=1.0, dy=dz, dx=dy, so specifying only one value sets
 * them all to that value, and specifying only two the second value sets both
 * dy and dx.
 *
 * Optional z0,y0,x0 arguments are the z,y,x grid indexes for the integration
 * zero point in the corresponding dimensions. These default to -1 which
 * means the mid-point of the field in that dimention.
 *
 */

// Integrate a Field2D along the x axis.
template<typename T, int X, int Y>
Field2D<T, X, Y> integratex(Field2D<T, X, Y>& f, const float dx = 1.0, Index x0 = -1)
{
    const Index      ylen = f.rows();
    const Index      xlen = f.cols();
    Field2D<T, X, Y> buf(ylen, xlen);
    if (x0 < 0)
        x0 = xlen / 2;
    // set p(y,x)=0 for x=0.
    buf.col(x0).setZero();
    // trapezoidal integrate f from x0 to xlen-1.
    for (Index x = x0 + 1; x < xlen; x++)
        buf.col(x) = buf.col(x - 1) + (0.5 * dx) * (f.col(x - 1) + f.col(x));
    // reverse-trapezoidal integrate f from x0 to 0.
    for (Index x = x0 - 1; x >= 0; x--)
        buf.col(x) = buf.col(x + 1) - (0.5 * dx) * (f.col(x + 1) + f.col(x));
    return buf;
}

// Integrate a Field2D along the y axis.
template<typename T, int X, int Y>
Field2D<T, X, Y> integratey(Field2D<T, X, Y>& f, const float dy = 1.0, Index y0 = -1)
{
    const Index      ylen = f.rows();
    const Index      xlen = f.cols();
    Field2D<T, X, Y> buf(ylen, xlen);
    if (y0 < 0)
        y0 = ylen / 2;
    // set p(y,x)=0 for y=y0.
    buf.row(y0).setZero();
    // trapezoidal integrate f from y0 to ylen-1.
    for (Index y = y0 + 1; y < ylen; y++)
        buf.row(y) = buf.row(y - 1) + (0.5 * dy) * (f.row(y - 1) + f.row(y));
    // reverse-trapezoidal integrate f from y0 to 0.
    for (Index y = y0 - 1; y >= 0; y--)
        buf.row(y) = buf.row(y + 1) - (0.5 * dy) * (f.row(y + 1) + f.row(y));
    return buf;
}

// Integrate a Field3D along the x axis.
template<typename T, int X, int Y, int Z>
Field3D<T, X, Y, Z> integratex(Field3D<T, X, Y, Z>& f, const float dx = 1.0, Index x0 = -1)
{
    const Index         zlen = f.zsize();
    Field3D<T, X, Y, Z> buf(zlen);
    // xintegrate all layers.
    for (Index z = 0; z < zlen; z++)
        buf(z) = integratex(f(z), dx, x0);
    return buf;
}

// Integrate a Field3D along the y axis.
template<typename T, int X, int Y, int Z>
Field3D<T, X, Y, Z> integratey(Field3D<T, X, Y, Z>& f, const float dy = 1.0, Index y0 = -1)
{
    const Index         zlen = f.zsize();
    Field3D<T, X, Y, Z> buf(zlen);
    // yintegrate all layers.
    for (Index z = 0; z < zlen; z++)
        buf(z) = integratey(f(z), dy, y0);
    return buf;
}

// Integrate a Field3D along the z axis.
template<typename T, int X, int Y, int Z>
Field3D<T, X, Y, Z> integratez(Field3D<T, X, Y, Z>& f, const float dz = 1.0, Index z0 = -1)
{
    const Index         zlen = f.zsize();
    const Index         ylen = f.ysize();
    const Index         xlen = f.xsize();
    Field3D<T, X, Y, Z> buf(zlen);
    // set p(z,y,x)=0 for z=0.
    buf(0) = Field2D<T, X, Y>::Zero(ylen, xlen);
    if (z0 < 0)
        z0 = zlen / 2;
    // set p(z,y,x)=0 for z=z0.
    buf(z0).setZero();
    // trapezoidal integrate f from z0 to zlen-1.
    for (Index z = z0 + 1; z < zlen; z++)
        buf(z) = buf(z - 1) + (0.5 * dz) * (f(z - 1) + f(z));
    // reverse-trapezoidal integrate f from z0 to 0.
    for (Index z = z0 - 1; z >= 0; z--)
        buf(z) = buf(z + 1) - (0.5 * dz) * (f(z + 1) + f(z));
    return buf;
}

// Integrate a Field3D along the z, y, and x axies into a Vec3Field3D.
template<typename T, int X, int Y, int Z>
Vec3Field3D<T, X, Y, Z> integratexyz(
    Field3D<T, X, Y, Z>& f, const float dz = 1.0, float dy = -1, float dx = -1, Index z0 = -1, Index y0 = -1, Index x0 = -1)
{
    const Index zlen = f.zsize();
    const Index ylen = f.ysize();
    const Index xlen = f.xsize();
    if (dy < 0)
        dy = dz;
    if (dx < 0)
        dx = dy;
    Vec3Field3D<T, X, Y, Z> buf(zlen);
    // Integratez for whole field, but for x and y we can do it by layer.
    auto bufz              = integratez(f, dz, z0);
    using MapVec3Field2Dim = Map<decltype(f(0)), 0, InnerStride<3>>;
    for (Index z = 0; z < zlen; z++) {
        buf(z)                                          = decltype(buf(0))(ylen, xlen);
        MapVec3Field2Dim(buf(z).data().x(), ylen, xlen) = integratex(f(z), dx, x0);
        MapVec3Field2Dim(buf(z).data().y(), ylen, xlen) = integratey(f(z), dy, y0);
        MapVec3Field2Dim(buf(z).data().z(), ylen, xlen) = bufz(z);
    }
    return buf;
}

} // namespace Slic3r::Geometry

#endif /* libslic3r_Geometry_Fields_hpp_ */
