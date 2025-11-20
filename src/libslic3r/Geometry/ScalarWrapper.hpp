#ifndef libslic3r_Geometry_ScalarWrapper_hpp_
#define libslic3r_Geometry_ScalarWrapper_hpp_

/***
 * An Eigen::ScalarWrapper for using eigen types like scalars.
 *
 * This means Matrix's and Array's can be treated like scalars and put inside
 * other Matrix's or Array's. This means you can use an eigen
 * expression designed for Bicubic interpolation of a Real scalar field to do
 * Bicubic interpolation of a Vector Field.
 */

#include <algorithm>
#include <vector>
#include <cmath>

#include <Eigen/Dense>
#include <Eigen/src/Core/IO.h>

/***
 * We define some extensions to Eigen here so we can put any type of
 * Eigen object inside another Eigen::Matrix or Eigen::Array. This means we
 * can interpolate not just scalar fields, but also vector fields.
 */
namespace Eigen {

template<typename ExpressionType> class ScalarWrapper;

// Define NumTraits<ScalarWrapper<> so we can put ScalarWrapper
// types inside other Matrix or Array types. This copies the existing
// NumTraits<Array<...>> definition.
template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
struct NumTraits<ScalarWrapper<Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>>>
{
    typedef _Scalar                                                              MatrixScalar;
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

    EIGEN_DEVICE_FUNC constexpr static RealScalar epsilon() { return NumTraits<RealScalar>::epsilon(); }
    EIGEN_DEVICE_FUNC constexpr static RealScalar dummy_precision() { return NumTraits<RealScalar>::dummy_precision(); }
    constexpr static int                          digits10() { return NumTraits<_Scalar>::digits10(); }
    constexpr static int                          max_digits10() { return NumTraits<_Scalar>::max_digits10(); }
};

template<typename ExpressionType, typename BinaryOp>
struct ScalarBinaryOpTraits<ScalarWrapper<ExpressionType>, typename ExpressionType::Scalar, BinaryOp>
{
    typedef ScalarWrapper<ExpressionType> ReturnType;
};
template<typename ExpressionType, typename BinaryOp>
struct ScalarBinaryOpTraits<typename ExpressionType::Scalar, ScalarWrapper<ExpressionType>, BinaryOp>
{
    typedef ScalarWrapper<ExpressionType> ReturnType;
};

/** \class ScalarWrapper
 * \ingroup Core_Module
 *
 * \brief Wrapper around Matrix of Array expressions so they can be used like a scalar.
 *
 * This class is the return type of scalar(), and most of the time this is the only way it is used.
 *
 * \sa scalar(), class MatrixWrapper, class ArrayWrapper.
 */

// Note we intentionally don't inherit from MatrixBase because we don't want
// to confuse scalar vs matrix template matching.
// class ScalarWrapper : public MatrixBase<ScalarWrapper<ExpressionType> > {
template<typename ExpressionType> class ScalarWrapper
{
public:
    // typedef MatrixBase<ScalarWrapper<ExpressionType> > Base;
    typedef typename ExpressionType::Scalar                                                      Scalar;
    typedef std::conditional_t<internal::is_lvalue<ExpressionType>::value, Scalar, const Scalar> ScalarWithConstIfNotLvalue;

    // Make all template instantiations friends so they can access m_expression of each other.
    template<typename OtherExpressionType> friend class ScalarWrapper;

    // EIGEN_DENSE_PUBLIC_INTERFACE(ScalarWrapper)
    // EIGEN_INHERIT_ASSIGNMENT_OPERATORS(ScalarWrapper)

    EIGEN_DEVICE_FUNC inline ScalarWrapper(void) : m_expression(){};
    // EIGEN_DEVICE_FUNC inline ScalarWrapper(const ScalarWrapper&) = default;
    template<typename OtherExpressionType>
    EIGEN_DEVICE_FUNC inline ScalarWrapper(const ScalarWrapper<OtherExpressionType>& other) : m_expression(other.m_expression)
    {}
    template<typename OtherDerived>
    explicit EIGEN_DEVICE_FUNC inline ScalarWrapper(const DenseBase<OtherDerived>& other) : m_expression(other.derived())
    {}

    // We need to be able to initialize from a Scalar for Zero(), Ones(), etc.
    explicit ScalarWrapper(const Scalar& value) : m_expression{ExpressionType::Constant(value)} {}

    template<typename OtherExpressionType>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ScalarWrapper& operator=(const ScalarWrapper<OtherExpressionType>& other)
    {
        m_expression = other.m_expression;
        return *this;
    }
    template<typename OtherExpressionType> EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ScalarWrapper& operator=(const OtherExpressionType& other)
    {
        m_expression = other;
        return *this;
    }

    EIGEN_DEVICE_FUNC constexpr Index rows() const noexcept { return m_expression.rows(); }
    EIGEN_DEVICE_FUNC constexpr Index cols() const noexcept { return m_expression.cols(); }
    // EIGEN_DEVICE_FUNC constexpr Index outerStride() const noexcept { return m_expression.outerStride(); }
    // EIGEN_DEVICE_FUNC constexpr Index innerStride() const noexcept { return m_expression.innerStride(); }

    EIGEN_DEVICE_FUNC constexpr ScalarWithConstIfNotLvalue* data() { return m_expression.data(); }
    EIGEN_DEVICE_FUNC constexpr const Scalar*               data() const { return m_expression.data(); }

    EIGEN_DEVICE_FUNC inline const Scalar& coeffRef(Index rowId, Index colId) const { return m_expression.coeffRef(rowId, colId); }
    EIGEN_DEVICE_FUNC inline const Scalar& coeffRef(Index index) const { return m_expression.coeffRef(index); }
    EIGEN_DEVICE_FUNC inline const Scalar& coeff(Index rowId, Index colId) const { return m_expression.coeff(rowId, colId); }
    EIGEN_DEVICE_FUNC inline const Scalar& coeff(Index index) const { return m_expression.coeff(index); }

    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const ExpressionType& nestedExpression() const { return m_expression; }
    // EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const internal::remove_all_t<NestedExpressionType>& matrix() const { return m_expression; }
    // EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const auto array() const { return m_expression.array(); }

    // EIGEN_DEVICE_FUNC void resize(Index newSize) { m_expression.resize(newSize); }
    // EIGEN_DEVICE_FUNC void resize(Index rows, Index cols) { m_expression.resize(rows, cols); }

    auto operator-() const { return ScalarWrapper(-m_expression); }

    template<typename OtherExpressionType> auto operator+(const ScalarWrapper<OtherExpressionType>& other) const
    { return ScalarWrapper(m_expression + other.m_expression); }
    auto        operator+(const Scalar& s) const { return ScalarWrapper(m_expression + s); }
    friend auto operator+(const Scalar& s, const ScalarWrapper& w) { return ScalarWrapper(s + w.m_expression); }

    template<typename OtherExpressionType> auto operator-(const ScalarWrapper<OtherExpressionType>& other) const
    { return ScalarWrapper(m_expression - other.m_expression); }
    auto        operator-(const Scalar& s) const { return ScalarWrapper(m_expression - s); }
    friend auto operator-(const Scalar& s, const ScalarWrapper& w) { return ScalarWrapper(s - w.m_expression); }

    template<typename OtherExpressionType> auto operator*(const ScalarWrapper<OtherExpressionType>& other) const
    { return ScalarWrapper(m_expression * other.m_expression); }
    auto        operator*(const Scalar& s) const { return ScalarWrapper(m_expression * s); }
    friend auto operator*(const Scalar& s, const ScalarWrapper& w) { return ScalarWrapper(s * w.m_expression); }

    template<typename OtherExpressionType> auto operator/(const ScalarWrapper<OtherExpressionType>& other) const
    { return ScalarWrapper(m_expression / other.m_expression); }
    auto        operator/(const Scalar& s) const { return ScalarWrapper(m_expression / s); }
    friend auto operator/(const Scalar& s, const ScalarWrapper& w) { return ScalarWrapper(s / w.m_expression); }

    template<typename OtherExpressionType> auto& operator+=(const ScalarWrapper<OtherExpressionType>& other)
    {
        m_expression += other.m_expression;
        return *this;
    }
    auto& operator+=(const Scalar& s)
    {
        m_expression += s;
        return *this;
    }

    template<typename OtherExpressionType> auto& operator-=(const ScalarWrapper<OtherExpressionType>& other)
    {
        m_expression -= other.m_expression;
        return *this;
    }
    auto& operator-=(const Scalar& s)
    {
        m_expression -= s;
        return *this;
    }

    template<typename OtherExpressionType> auto& operator*=(const ScalarWrapper<OtherExpressionType>& other)
    {
        m_expression *= other.m_expression;
        return *this;
    }
    auto& operator*=(const Scalar& s)
    {
        m_expression *= s;
        return *this;
    }

    template<typename OtherExpressionType> auto& operator/=(const ScalarWrapper<OtherExpressionType>& other)
    {
        m_expression /= other.m_expression;
        return *this;
    }
    auto& operator/=(const Scalar& s)
    {
        m_expression /= s;
        return *this;
    }

    template<typename OtherExpressionType> auto operator==(const ScalarWrapper<OtherExpressionType>& other) const
    { return m_expression == other.m_expression; }
    auto operator==(const Scalar& s) const { return m_expression == s; }

    template<typename OtherExpressionType> auto operator!=(const ScalarWrapper<OtherExpressionType>& other) const
    { return m_expression != other.m_expression; }
    auto operator!=(const Scalar& s) const { return m_expression != s; }

    template<typename OtherExpressionType> auto operator<=(const ScalarWrapper<OtherExpressionType>& other) const
    { return m_expression <= other.m_expression; }
    auto operator<=(const Scalar& s) const { return m_expression <= s; }

    template<typename OtherExpressionType> auto operator<(const ScalarWrapper<OtherExpressionType>& other) const
    { return m_expression < other.m_expression; }
    auto operator<(const Scalar& s) const { return m_expression < s; }

    template<typename OtherExpressionType> auto operator>=(const ScalarWrapper<OtherExpressionType>& other) const
    { return m_expression >= other.m_expression; }
    auto operator>=(const Scalar& s) const { return m_expression >= s; }

    template<typename OtherExpressionType> auto operator>(const ScalarWrapper<OtherExpressionType>& other) const
    { return m_expression > other.m_expression; }
    auto operator>(const Scalar& s) const { return m_expression > s; }

    const WithFormat<ExpressionType> format(const IOFormat& fmt) const { return m_expression.format(fmt); }

protected:
    ExpressionType m_expression;
};

// These wrap an expression as a scalar with ScalarWrapper.
template<typename ExpressionType> EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ScalarWrapper<ExpressionType> scalar(ExpressionType& m)
{ return ScalarWrapper<ExpressionType>(m); }
template<typename ExpressionType>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const ScalarWrapper<const ExpressionType> scalar(const ExpressionType& m)
{ return ScalarWrapper<const ExpressionType>(m); }

template<typename ExpressionType> auto real(const ScalarWrapper<ExpressionType>& x) { return scalar(real(x.m_expression)); }
template<typename ExpressionType> auto imag(const ScalarWrapper<ExpressionType>& x) { return scalar(imag(x.m_expression)); }
template<typename ExpressionType> auto conj(const ScalarWrapper<ExpressionType>& x) { return scalar(conj(x.m_expression)); }
template<typename ExpressionType> auto sqrt(const ScalarWrapper<ExpressionType>& x) { return scalar(sqrt(x.m_expression)); }
template<typename ExpressionType> auto abs(const ScalarWrapper<ExpressionType>& x) { return scalar(abs(x.m_expression)); }
template<typename ExpressionType> auto cos(const ScalarWrapper<ExpressionType>& x) { return scalar(cos(x.m_expression)); }
template<typename ExpressionType> auto sin(const ScalarWrapper<ExpressionType>& x) { return scalar(sin(x.m_expression)); }
template<typename ExpressionType> auto acos(const ScalarWrapper<ExpressionType>& x) { return scalar(acos(x.m_expression)); }
template<typename ExpressionType> auto atan2(const ScalarWrapper<ExpressionType>& y, const ScalarWrapper<ExpressionType>& x)
{ return scalar(atan2(y.m_expression, x.m_expression)); }

template<typename ExpressionType> std::ostream& operator<<(std::ostream& stream, const ScalarWrapper<ExpressionType>& x)
{
    IOFormat fmt;
    if (x.rows() == 1) {
        fmt = IOFormat(StreamPrecision, DontAlignCols, ",", ",", "{", "}", "", "");
    } else if (x.cols() == 1) {
        fmt = IOFormat(StreamPrecision, DontAlignCols, ",", ",", "", "", "{", "}");
    } else {
        fmt = IOFormat(StreamPrecision, DontAlignCols, ",", ",", "{", "}", "{", "}");
    }
    stream << x.format(fmt);
    return stream;
}

} // namespace Eigen

#endif /* libslic3r_Geometry_ScalarWrapper_hpp_ */
