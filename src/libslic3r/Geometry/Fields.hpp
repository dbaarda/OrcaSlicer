#ifndef libslic3r_Geometry_Fields_hpp_
#define libslic3r_Geometry_Fields_hpp_

#include <algorithm>
#include <vector>
#include <cmath>

#include <Eigen/Dense>
#include "Bicubic.hpp"

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

// An p(x,y,z) Vector of type T;
template<typename T> using Vec3 = Eigen::Matrix<T, 3, 1>;

// A 2D Field of type T.
template<typename T, int Y = Eigen::Dynamic, int X = Eigen::Dynamic> using Field2 = Eigen::Matrix<T, Y, X, Eigen::RowMajor>;

// A 3D Field of type T.
//
// Note this is implemented as a subclass because nesting a matrix as a
// custom scalar type inside a vector didn't just work. So I created this
// subclass and just added redefinitions of any methods and operators that
// didn't compile or work properly. There is almost certainly a better way to
// do this, maybe using the unsupported Eigen::Tensor stuff, but that isn't
// included with the old v3.3.7 Eigen included in OrcaSlicer.
template<typename T, int Z = Eigen::Dynamic, int Y = Eigen::Dynamic, int X = Eigen::Dynamic>
class Field3 : public Eigen::Vector<Field2<T, Y, X>, Z>
{
public:
    Field3(void) : Eigen::Vector<Field2<T, Y, X>, Z>() {}

    Field3(const Eigen::Index z, const Eigen::Index y = 0, const Eigen::Index x = 0) : Eigen::Vector<Field2<T, Y, X>, Z>(z)
    {
        if (y > 0 && x > 0) {
            for (Eigen::Index i = 0; i < z; i++)
                (*this)(i) = Field2<T, Y, X>(y, x);
        }
    }

    template<typename OtherDerived> Field3(const Eigen::MatrixBase<OtherDerived>& other) : Eigen::Vector<Field2<T, Y, X>, Z>(other) {}

    template<typename OtherDerived> Field3& operator=(const Eigen::MatrixBase<OtherDerived>& other)
    {
        this->Eigen::Vector<Field2<T, Y, X>, Z>::operator=(other);
        return *this;
    }

    Field3& operator*=(const T v)
    {
        for (Eigen::Index i = 0; i < this->size(); i++)
            (*this)(i) *= v;
        return *this;
    }

    auto sum()
    {
        auto s = (*this)(0);
        for (Eigen::Index i = 1; i < this->size(); i++)
            s += (*this)(i);
        return s;
    }
};

// A 2D Vector Field of Vec3 of type T.
template<typename T, int Y = Eigen::Dynamic, int X = Eigen::Dynamic> using Vec3Field2 = Field2<Vec3<T>, Y, X>;

// A 3D Vector Field of Vec3 of type T.
template<typename T, int Z = Eigen::Dynamic, int Y = Eigen::Dynamic, int X = Eigen::Dynamic> using Vec3Field3 = Field3<Vec3<T>, Z, Y, X>;

/*********************************************
 * Interpolate Fields.
 *
 * These do linear or cubic interpolation to get field values at arbitrary x,
 * p(x,y), or p(x,y,z) points within the field grid. The x,y,z values should
 * be floats with values in the range of [0, <size>) where <size> is f.size()
 * for z, f.rows() for y, or f.cols() for x. It will interpolate a value
 * between the enclosing grid points.
 */
// Old Eigen versions don't have this defined.
// using RowVector2d = Eigen::Matrix<double, 1, 2>;

// A simple linear interpolation "kernel" with a similar API to Bicubic.hpp.
namespace LerpKernel {

const Eigen::RowVector2d cint(const double x) { return Eigen::RowVector2d(1 - x, x); }

// Get a 2x1 vector from a vector for lerp.
template<typename Derived> const auto fblock(const Eigen::EigenBase<Derived>& F, const Eigen::Index x)
{
    assert(0 <= x && x < F.size());
    return F.template segment<2>(x);
}

// Get a 2x2 block from a matrix for lerp.
template<typename Derived> const auto fblock(const Eigen::EigenBase<Derived>& F, const Eigen::Index y, const Eigen::Index x)
{
    assert(0 <= y && y < F.rows());
    assert(0 <= x && x < F.cols());
    return F.template block<2, 2>(x, y);
}

} // namespace LerpKernel

// linearly interpolate within a Field2 to get the value at p(x,y).
template<typename P, typename T, int Y = Eigen::Dynamic, int X = Eigen::Dynamic> const T lerp(const Field2<T, Y, X>& f, const P& p)
{
    const Eigen::Index iy = std::floor(p.y());
    const Eigen::Index ix = std::floor(p.x());
    const double       ry = p.y() - iy;
    const double       rx = p.x() - ix;
    return LerpKernel::cint(ry) * LerpKernel::fblock(f, iy, ix) * LerpKernel::cint(rx).transpose();
}

// linearly interpolate within a Field3 to get the value at p(x,y,z).
template<typename P, typename T, int Z = Eigen::Dynamic, int Y = Eigen::Dynamic, int X = Eigen::Dynamic>
const T lerp(const Field3<T, Z, Y, X>& f, const P& p)
{
    const Eigen::Index iz = std::floor(p.z());
    const Eigen::Index iy = std::floor(p.y());
    const Eigen::Index ix = std::floor(p.x());
    const double       rz = p.z() - iz;
    const double       ry = p.y() - iy;
    const double       rx = p.x() - ix;
    assert(0 <= iz && iz < f.size());
    assert(0 <= iy && iy < f(0).rows());
    assert(0 <= ix && ix < f(0).cols());
    Eigen::Vector2d cy = LerpKernel::cint(ry);
    Eigen::Vector2d cx = LerpKernel::cint(rx).transpose();
    double          z0 = cy * LerpKernel::fblock(f(iz), iy, ix) * cx;
    double          z1 = cy * LerpKernel::fblock(f(iz + 1), iy, ix) * cx;
    return LerpKernel::cint(rz) * Eigen::Vector2d(z0, z1);
}

// cubic interpolate within a Field2 to get the value at p(x,y).
template<typename P, typename T, int Y = Eigen::Dynamic, int X = Eigen::Dynamic> const T cubic(const Field2<T, Y, X>& f, const P& p)
{
    using Scalar          = typename P::Scalar;
    using Kernel          = CubicCatmulRomKernel<Scalar, T>;
    const Eigen::Index iy = std::floor(p.y());
    const Eigen::Index ix = std::floor(p.x());
    const double       ry = p.y() - iy;
    const double       rx = p.x() - ix;
    assert(0 <= iy && iy < f.rows());
    assert(0 <= ix && ix < f.cols());
    auto cy = Kernel::cint(ry);
    auto cx = Kernel::cint(rx).transpose().eval();
    return cy * Kernel::fblock(f, iy, ix) * cx;
}

// cubic interpolate within a Field3 to get the value at p(x,y,z).
template<typename P, typename T, int Z = Eigen::Dynamic, int Y = Eigen::Dynamic, int X = Eigen::Dynamic>
const T cubic(const Field3<T, Z, Y, X>& f, const P& p)
{
    using Scalar          = typename P::Scalar;
    using Kernel          = CubicCatmulRomKernel<Scalar, T>;
    const Eigen::Index iz = std::floor(p.z());
    const Eigen::Index iy = std::floor(p.y());
    const Eigen::Index ix = std::floor(p.x());
    const Scalar       rz = p.z() - iz;
    const Scalar       ry = p.y() - iy;
    const Scalar       rx = p.x() - ix;
    assert(0 <= iz && iz < f.size());
    assert(0 <= iy && iy < f(0).rows());
    assert(0 <= ix && ix < f(0).cols());
    const auto                cy = Kernel::cint(ry);
    const auto                cx = Kernel::cint(rx).transpose().eval();
    typename Kernel::Vector4V zf;
    for (int z = 0; z < 4; z++)
        zf(z) = cy * Kernel::fblock(f(std::clamp<Eigen::Index>(iz + z - 1, 0, f.size() - 1)), iy, ix) * cx;
    return Kernel::cint(rz) * zf;
}

// cubic interpolate a Field3 into a Field2 at a given z.
template<typename S, typename T, int Z = Eigen::Dynamic, int Y = Eigen::Dynamic, int X = Eigen::Dynamic>
const Field2<T, Y, X> cubic_z(const Field3<T, Z, Y, X>& f, const S z)
{
    using Kernel = CubicCatmulRomKernel<S, Field2<T, Y, X>>;
    return Kernel::interpolate(f, z);
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

// Do a single iteration of summing 3x3 adjacent values of Field2.
template<typename T, int Y, int X> void sum33(Field2<T, Y, X>& f)
{
    const Eigen::Index y = f.rows();
    const Eigen::Index x = f.cols();
    Field2<T, Y, X>    buf(y, x);
    // TODO: this could be done using only 3 temporary rows as a circular
    // buffer and cycling through them. Maybe this would be cheaper?
    // First fill buf with the sums of 3 adjacent columns.
    // Set outer-edge cols with the edge value extended.
    buf.col(0)     = 2.0 * f.col(0) + f.col(1);
    buf.col(x - 1) = 2.0 * f.col(x - 1) + f.col(x - 2);
    // Set inner block to sum of 3 col values.
    buf.middleCols(1, x - 2) = f.leftCols(x - 2) + f.middleCols(1, x - 2) + f.rightCols(x - 2);
    // next fill f with the sums of 3 rows from buf.
    // Set outer-edge rows with the edge value extended
    f.row(0)     = 2.0 * buf.row(0) + buf.row(1);
    f.row(y - 1) = 2.0 * buf.row(y - 1) + buf.row(y - 2);
    // Set inner block to sum of 3 row values.
    f.middleRows(1, y - 2) = buf.topRows(y - 2) + buf.middleRows(1, y - 2) + buf.bottomRows(y - 2);
}

// GausianSmooth (almost) a Field2 with multiple iterations of 3x3 averaging.
template<typename T, int Y = Eigen::Dynamic, int X = Eigen::Dynamic> void smooth(Field2<T, Y, X>& f, const int n = 6)
{
    double m = 1.0;
    for (int i = 0; i < n; i++) {
        sum33<T, Y, X>(f);
        m *= 9.0;
    }
    f *= T(1.0 / m);
}

// Do a single iteration of summing 3 adjacent layers of Field3.
template<typename T, int Z, int Y, int X> void sum3(Field3<T, Z, Y, X>& f)
{
    const Eigen::Index z = f.size();
    const Eigen::Index y = f(0).rows();
    const Eigen::Index x = f(0).cols();
    assert(z > 0);
    assert(y > 0);
    assert(x > 0);
    // Used to hold 3 layers i-1, i, i+1 for summing, with layer i at index i%3.
    Field3<T, 3, Y, X> buf(3, y, x);
    // Save the first 2 layers into buf for later sums.
    buf.head(2) = f.head(2);
    // Set outer-edge layer 0 with the edge value extended.
    f(0) = 2.0 * f(0) + f(1);
    // Go through layers copying the i+1 layer into buf and put the sum/3 into f.
    for (Eigen::Index i = 1; i < z - 1; i++) {
        buf((i + 1) % 3) = f(i + 1);
        f(i)             = buf.sum();
    }
    // Set outer-edge layer z-1 with the edge value extended.
    f(z - 1) = 2.0 * f(z - 1) + buf((z - 2) % 3);
}

// GausianSmooth (almost) a Field3 with multiple iterations of 3x3x3 averaging.
template<typename T, int Z = Eigen::Dynamic, int Y = Eigen::Dynamic, int X = Eigen::Dynamic>
void smooth(Field3<T, Z, Y, X>& f, const int n = 6)
{
    double             m    = 1.0;
    const Eigen::Index zlen = f.size();
    // First smooth each layer n times.
    for (Eigen::Index z = 0; z < zlen; z++)
        smooth<T, Y, X>(f(z), n);
    // Now sum 3 adjacent layers n times.
    for (int i = 0; i < n; i++) {
        sum3<T, Z, Y, X>(f);
        m *= 3.0;
    }
    // Now divide the sums by 3^n to get the average.
    f *= T(1 / m);
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

// Integrate a Field2 along the x axis.
template<typename T, int Y = Eigen::Dynamic, int X = Eigen::Dynamic>
Field2<T, Y, X> integratex(Field2<T, Y, X>& f, const float dx = 1.0, Eigen::Index x0 = -1)
{
    const Eigen::Index ylen = f.rows();
    const Eigen::Index xlen = f.cols();
    Field2<T, Y, X>    buf(ylen, xlen);
    if (x0 < 0)
        x0 = xlen / 2;
    // set p(y,x)=0 for x=0.
    buf.col(x0).setZero();
    // trapezoidal integrate f from x0 to xlen-1.
    for (Eigen::Index x = x0 + 1; x < xlen; x++)
        buf.col(x) = buf.col(x - 1) + (0.5 * dx) * (f.col(x - 1) + f.col(x));
    // reverse-trapezoidal integrate f from x0 to 0.
    for (Eigen::Index x = x0 - 1; x >= 0; x--)
        buf.col(x) = buf.col(x + 1) - (0.5 * dx) * (f.col(x + 1) + f.col(x));
    return buf;
}

// Integrate a Field2 along the y axis.
template<typename T, int Y = Eigen::Dynamic, int X = Eigen::Dynamic>
Field2<T, Y, X> integratey(Field2<T, Y, X>& f, const float dy = 1.0, Eigen::Index y0 = -1)
{
    const Eigen::Index ylen = f.rows();
    const Eigen::Index xlen = f.cols();
    Field2<T, Y, X>    buf(ylen, xlen);
    if (y0 < 0)
        y0 = ylen / 2;
    // set p(y,x)=0 for y=y0.
    buf.row(y0).setZero();
    // trapezoidal integrate f from y0 to ylen-1.
    for (Eigen::Index y = y0 + 1; y < ylen; y++)
        buf.row(y) = buf.row(y - 1) + (0.5 * dy) * (f.row(y - 1) + f.row(y));
    // reverse-trapezoidal integrate f from y0 to 0.
    for (Eigen::Index y = y0 - 1; y >= 0; y--)
        buf.row(y) = buf.row(y + 1) - (0.5 * dy) * (f.row(y + 1) + f.row(y));
    return buf;
}

// Integrate a Field3 along the x axis.
template<typename T, int Z = Eigen::Dynamic, int Y = Eigen::Dynamic, int X = Eigen::Dynamic>
Field3<T, Z, Y, X> integratex(Field3<T, Z, Y, X>& f, const float dx = 1.0, Eigen::Index x0 = -1)
{
    const Eigen::Index zlen = f.size();
    Field3<T, Z, Y, X> buf(zlen);
    // xintegrate all layers.
    for (Eigen::Index z = 0; z < zlen; z++)
        buf(z) = integratex(f(z), dx, x0);
    return buf;
}

// Integrate a Field3 along the y axis.
template<typename T, int Z = Eigen::Dynamic, int Y = Eigen::Dynamic, int X = Eigen::Dynamic>
Field3<T, Z, Y, X> integratey(Field3<T, Z, Y, X>& f, const float dy = 1.0, Eigen::Index y0 = -1)
{
    const Eigen::Index zlen = f.size();
    Field3<T, Z, Y, X> buf(zlen);
    // yintegrate all layers.
    for (Eigen::Index z = 0; z < zlen; z++)
        buf(z) = integratey(f(z), dy, y0);
    return buf;
}

// Integrate a Field3 along the z axis.
template<typename T, int Z = Eigen::Dynamic, int Y = Eigen::Dynamic, int X = Eigen::Dynamic>
Field3<T, Z, Y, X> integratez(Field3<T, Z, Y, X>& f, const float dz = 1.0, Eigen::Index z0 = -1)
{
    const Eigen::Index zlen = f.size();
    const Eigen::Index ylen = f(0).rows();
    const Eigen::Index xlen = f(0).cols();
    Field3<T, Z, Y, X> buf(zlen);
    // set p(z,y,x)=0 for z=0.
    buf(0) = Field2<T, Y, X>::Zero(ylen, xlen);
    if (z0 < 0)
        z0 = zlen / 2;
    // set p(z,y,x)=0 for z=z0.
    buf(z0).setZero();
    // trapezoidal integrate f from z0 to zlen-1.
    for (Eigen::Index z = z0 + 1; z < zlen; z++)
        buf(z) = buf(z - 1) + (0.5 * dz) * (f(z - 1) + f(z));
    // reverse-trapezoidal integrate f from z0 to 0.
    for (Eigen::Index z = z0 - 1; z >= 0; z--)
        buf(z) = buf(z + 1) - (0.5 * dz) * (f(z + 1) + f(z));
    return buf;
}

// Integrate a Field3 along the z, y, and x axies into a Vec3Field3.
template<typename T, int Z = Eigen::Dynamic, int Y = Eigen::Dynamic, int X = Eigen::Dynamic>
Vec3Field3<T, Z, Y, X> integratexyz(Field3<T, Z, Y, X>& f,
                                    const float         dz = 1.0,
                                    float               dy = -1,
                                    float               dx = -1,
                                    Eigen::Index        z0 = -1,
                                    Eigen::Index        y0 = -1,
                                    Eigen::Index        x0 = -1)
{
    const Eigen::Index zlen = f.size();
    const Eigen::Index ylen = f(0).rows();
    const Eigen::Index xlen = f(0).cols();
    if (dy < 0)
        dy = dz;
    if (dx < 0)
        dx = dy;
    Vec3Field3<T, Z, Y, X> buf(zlen);
    // Integratez for whole field, but for x and y we can do it by layer.
    auto bufz              = integratez(f, dz, z0);
    using MapVec3Field2Dim = Eigen::Map<decltype(f(0)), 0, Eigen::InnerStride<3>>;
    for (Eigen::Index z = 0; z < zlen; z++) {
        buf(z)                                          = decltype(buf(0))(ylen, xlen);
        MapVec3Field2Dim(buf(z).data().x(), ylen, xlen) = integratex(f(z), dx, x0);
        MapVec3Field2Dim(buf(z).data().y(), ylen, xlen) = integratey(f(z), dy, y0);
        MapVec3Field2Dim(buf(z).data().z(), ylen, xlen) = bufz(z);
    }
    return buf;
}

} // namespace Slic3r::Geometry

#endif /* libslic3r_Geometry_Fields_hpp_ */
