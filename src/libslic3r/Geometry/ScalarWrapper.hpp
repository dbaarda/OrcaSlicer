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

namespace internal {
// TODO: we need to make plain_matrix_type<D> work for ScalarWrappers around expressions.
}

template<typename ExpressionType>
class ScalarWrapper;

// Used for getting a type with fallback to a default.
template<typename T, typename Default>
using get_type_t = typename std::conditional_t<std::is_void_v<T>, Default, T>;

// Used for getting the Eigen Base classes for a derived class.
template<typename D>
using get_dense_base_t = DenseBase<std::decay_t<D>>;
template<typename D>
using get_matrix_base_t = MatrixBase<std::decay_t<D>>;
template<typename D>
using get_array_base_t = ArrayBase<std::decay_t<D>>;

/* These are used for constraining the types of forwarding reference template args.
 *
 *  template<T, B=enable_dense_t<T>>
 *  auto func(T&& v) {...};
 */
// enable if T is derived from Base, and return Base.
template<typename T, typename Base>
using enable_base_t = typename std::enable_if_t<std::is_base_of_v<Base, std::remove_reference_t<T>>, Base>;

// enable if T is the same as Type, and return Type.
template<typename T, typename Type>
using enable_type_t = std::enable_if_t<std::is_same_v<Type, std::remove_reference_t<T>>, Type>;

// Enable if T has the right eigen dense, matrix, or array base class.
template<typename T>
using enable_dense_t = enable_base_t<T, get_dense_base_t<T>>;
template<typename T>
using enable_matrix_t = enable_base_t<T, get_matrix_base_t<T>>;
template<typename T>
using enable_array_t = enable_base_t<T, get_array_base_t<T>>;

// Helpers for finding the inner-most Scalar type for multiply-nested ScalarWrappers.
template<typename Scalar, typename Inner = void>
struct get_inner_scalar
{
    using type = get_type_t<Inner, Scalar>;
};

template<typename Scalar>
struct get_inner_scalar<Scalar, typename NumTraits<Scalar>::XprInnerScalar>;

template<typename Scalar>
using get_inner_scalar_t = get_inner_scalar<Scalar>::type;

// Define NumTraits<ScalarWrapper<> so we can put ScalarWrapper types inside
// other Matrix or Array types. This is based on the existing
// NumTraits<Array<...>> definition.
//
// Note ScalarWrapper can wrap reference types, so we sometimes need to strip
// the reference for some traits lookups.
template<typename _XprType>
struct NumTraits<ScalarWrapper<_XprType>>
{
    typedef _XprType                                                                     XprType;
    typedef XprType::Scalar                                                              XprScalar;
    typedef typename NumTraits<XprScalar>::Real                                          XprScalarReal;
    typedef typename NumTraits<XprScalar>::NonInteger                                    XprScalarNonInteger;
    typedef get_inner_scalar_t<XprScalar>                                                XprInnerScalar;
    typedef ScalarWrapper<typename XprType::RealReturnType>                              Real;
    typedef ScalarWrapper<typename XprType::template CastXpr<XprScalarNonInteger>::Type> NonInteger;
    typedef ScalarWrapper<XprType>&                                                      Nested;
    typedef typename NumTraits<XprScalar>::Literal                                       Literal;

    enum {
        IsComplex             = NumTraits<XprScalar>::IsComplex,
        IsInteger             = NumTraits<XprScalar>::IsInteger,
        IsSigned              = NumTraits<XprScalar>::IsSigned,
        RequireInitialization = 1,
        ReadCost = XprType::SizeAtCompileTime == Dynamic ? HugeCost : XprType::SizeAtCompileTime * int(NumTraits<XprScalar>::ReadCost),
        AddCost  = XprType::SizeAtCompileTime == Dynamic ? HugeCost : XprType::SizeAtCompileTime * int(NumTraits<XprScalar>::AddCost),
        MulCost  = XprType::SizeAtCompileTime == Dynamic ? HugeCost : XprType::SizeAtCompileTime * int(NumTraits<XprScalar>::MulCost)
    };

    EIGEN_DEVICE_FUNC constexpr static XprScalarReal epsilon() { return NumTraits<XprScalarReal>::epsilon(); }
    EIGEN_DEVICE_FUNC constexpr static XprScalarReal dummy_precision() { return NumTraits<XprScalarReal>::dummy_precision(); }

    constexpr static int digits10() { return NumTraits<XprScalar>::digits10(); }
    constexpr static int max_digits10() { return NumTraits<XprScalar>::max_digits10(); }
};

// The wrapped reference NumTraits inherits the non-reference version.
template<typename _XprType>
struct NumTraits<ScalarWrapper<_XprType&>> : NumTraits<ScalarWrapper<_XprType>>
{};

// Add support for scalar binary operations between a non-reference ScalarWrapper and a Scalar.
template<typename ExpressionType, typename BinaryOp>
struct ScalarBinaryOpTraits<ScalarWrapper<ExpressionType>, typename ExpressionType::Scalar, BinaryOp>
{
    typedef ScalarWrapper<ExpressionType> ReturnType;
};

// Add support for scalar binary operations between a Scalar and a non-reference ScalarWrapper.
template<typename ExpressionType, typename BinaryOp>
struct ScalarBinaryOpTraits<typename ExpressionType::Scalar, ScalarWrapper<ExpressionType>, BinaryOp>
{
    typedef ScalarWrapper<ExpressionType> ReturnType;
};

// Add support for scalar binary operations between a reference ScalarWrapper and a Scalar.
template<typename XprType, typename BinaryOp>
struct ScalarBinaryOpTraits<ScalarWrapper<XprType&>, typename XprType::Scalar, BinaryOp>
    : ScalarBinaryOpTraits<ScalarWrapper<XprType>, typename XprType::Scalar, BinaryOp>
{};

// Add support for scalar binary operations between a Scalar and a reference ScalarWrapper.
template<typename XprType, typename BinaryOp>
struct ScalarBinaryOpTraits<typename XprType::Scalar, ScalarWrapper<XprType&>, BinaryOp>
    : ScalarBinaryOpTraits<typename XprType::Scalar, ScalarWrapper<XprType>, BinaryOp>
{};

/** \class ScalarWrapper
 * \ingroup Core_Module
 *
 * \brief Wrapper around Matrix or Array expressions so they can be used like a scalar.
 *
 * This class is the return type of scalar(), and most of the time this is the only way it is used.
 *
 * There are the following kinds of wrappers;
 *
 *   * a non-const reference wrapper: a multable reference to another lvalue expression.
 *   * a const reference wrapper: a non-mutable reference to another lvalue expression.
 *   * a non-cost instance wrapper: a mutable wrapped expression.
 *   * a const instance wrapper: a non-multable wrapped expression.
 *
 * Generally, reference wrappers are views of another expression or matrix,
 * and instance wrappers have their own data. However, eigen expressions can
 * themselves be views of other expressions, so an instance view of a block
 * expression will be a view of another expression.
 *
 * In general, wrappers created from lvalue expressions will be reference
 * wrappers, and wrappers created from rvalues will be instance wrappers.
 * This should ensure that wrappers created around temporary rvalue
 * expressions will contain (moved) copies of those expressions, but wrappers
 * around persistent lvalue expressions will be lightweight references.
 *
 * Wrappers created by scalar expressions of other wrappers will be instance
 * wrappers around the eigen expressions those operations generate. These
 * expressions are not yet evaluated, and will be evaluated when the wrapper
 * is assigned to another wrapper around an expression to contain the result.
 *
 * Once a wrapper is created its expression type can't be changed,
 * Assignments and mutations will be applyed to the instance or reference
 * wrapped expression.
 *
 * \sa scalar(), class MatrixWrapper, class ArrayWrapper, class NestByValue.
 */

// Note we intentionally don't inherit from MatrixBase because we don't want
// to confuse scalar vs matrix template matching.
template<typename ExpressionType>
class ScalarWrapper
{
public:
    // typedef MatrixBase<ScalarWrapper<ExpressionType> > Base;
    typedef std::remove_reference_t<ExpressionType>   NonRefExpressionType;
    typedef std::remove_const_t<NonRefExpressionType> CleanExpressionType;
    typedef typename NonRefExpressionType::Scalar     Scalar;

    // Make this template a friend so all template instances can access m_expression of each other.
    template<typename OtherExpressionType>
    friend class ScalarWrapper;

    // The default empty constructor.
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE constexpr ScalarWrapper() = default;
    // The default copy constructor.
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ScalarWrapper(const ScalarWrapper&) = default;
    // Initialize from an initializer list.
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ScalarWrapper(std::initializer_list<Scalar> list) : m_expression{list} {}

    // Initialize implicitly from the right ExpressionType.
    // EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ScalarWrapper(const ExpressionType& expr) : m_expression(expr) {}

    // Initialize from an expression. This is a forwarding reference
    // initializer for references to lvalues and copies of rvalues.
    template<typename ArgXprType, typename BaseXprType = enable_dense_t<ArgXprType>>
    EIGEN_DEVICE_FUNC explicit EIGEN_STRONG_INLINE ScalarWrapper(ArgXprType&& expr) : m_expression(expr)
    {}

    // Copy a different ScalarWrapper.
    template<typename OtherExpressionType>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ScalarWrapper(const ScalarWrapper<OtherExpressionType>& other) : m_expression(other.m_expression)
    {}

    // We need to be able to initialize from a Scalar for Zero(), Ones(), etc.
    explicit ScalarWrapper(const Scalar& value)
        : m_expression{ExpressionType::Constant(std::max<Eigen::Index>(1, ExpressionType::RowsAtCompileTime),
                                                std::max<Eigen::Index>(1, ExpressionType::ColsAtCompileTime),
                                                value)}
    {}
    explicit ScalarWrapper(Eigen::Index size, const Scalar& value) : m_expression{ExpressionType::Constant(size, value)} {}
    explicit ScalarWrapper(Eigen::Index rows, Eigen::Index cols, const Scalar& value)
        : m_expression{ExpressionType::Constant(rows, cols, value)}
    {}

    // Default assignment operators.
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ScalarWrapper& operator=(const ScalarWrapper& other)
    {
        m_expression = other.m_expression;
        return *this;
    }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ScalarWrapper& operator=(ScalarWrapper&& other)
    {
        m_expression = other.m_expression;
        return *this;
    }

    // Assign from a ScalarWrapper with a different ExpressionType.
    template<typename OtherExpressionType>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ScalarWrapper& operator=(const ScalarWrapper<OtherExpressionType>& other)
    {
        m_expression = other.m_expression;
        return *this;
    }
    // Assign From any kind of ExpressionType.
    template<typename OtherExpressionType>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ScalarWrapper& operator=(const OtherExpressionType& other)
    {
        m_expression = other;
        return *this;
    }

    EIGEN_DEVICE_FUNC constexpr Index rows() const noexcept { return m_expression.rows(); }
    EIGEN_DEVICE_FUNC constexpr Index cols() const noexcept { return m_expression.cols(); }
    inline const auto&                operator()(const Eigen::Index r) const { return m_expression(r); }
    inline const auto&                operator()(const Eigen::Index r, const Eigen::Index c) const { return m_expression(r, c); }
    inline auto&                      operator()(const Eigen::Index r) { return m_expression(r); }
    inline auto&                      operator()(const Eigen::Index r, const Eigen::Index c) { return m_expression(r, c); }

    // These are various different ways to unwrap the expression.
                                                      operator const ExpressionType&() const { return m_expression; }
    EIGEN_DEVICE_FUNC const CleanExpressionType&      nestedExpression() const { return m_expression; }
    EIGEN_DEVICE_FUNC CleanExpressionType&            nestedExpression() { return m_expression; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const auto& derived() const { return m_expression.derived(); }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE auto&       derived() { return m_expression.derived(); }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto) matrix() const { return m_expression.matrix(); }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto) matrix() { return m_expression.matrix(); }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto) array() const { return m_expression.array(); }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto) array() { return m_expression.array(); }

    auto operator-() const { return scalar(-m_expression); }

    template<typename OtherExpressionType>
    auto operator+(const ScalarWrapper<OtherExpressionType>& other) const
    { return scalar(m_expression + other.m_expression); }
    auto        operator+(const Scalar& s) const { return scalar(m_expression + s); }
    friend auto operator+(const Scalar& s, const ScalarWrapper& w) { return scalar(s + w.m_expression); }

    template<typename OtherExpressionType>
    auto operator-(const ScalarWrapper<OtherExpressionType>& other) const
    { return scalar(m_expression - other.m_expression); }
    auto        operator-(const Scalar& s) const { return scalar(m_expression - s); }
    friend auto operator-(const Scalar& s, const ScalarWrapper& w) { return scalar(s - w.m_expression); }

    template<typename OtherExpressionType>
    auto operator*(const ScalarWrapper<OtherExpressionType>& other) const
    { return scalar(m_expression * other.m_expression); }
    auto        operator*(const Scalar& s) const { return scalar(m_expression * s); }
    friend auto operator*(const Scalar& s, const ScalarWrapper& w) { return scalar(s * w.m_expression); }

    template<typename OtherExpressionType>
    auto operator/(const ScalarWrapper<OtherExpressionType>& other) const
    { return scalar(m_expression / other.m_expression); }
    auto        operator/(const Scalar& s) const { return scalar(m_expression / s); }
    friend auto operator/(const Scalar& s, const ScalarWrapper& w) { return scalar(s / w.m_expression); }

    template<typename OtherExpressionType>
    auto& operator+=(const ScalarWrapper<OtherExpressionType>& other)
    {
        m_expression += other.m_expression;
        return *this;
    }
    auto& operator+=(const Scalar& s)
    {
        m_expression += s;
        return *this;
    }

    template<typename OtherExpressionType>
    auto& operator-=(const ScalarWrapper<OtherExpressionType>& other)
    {
        m_expression -= other.m_expression;
        return *this;
    }
    auto& operator-=(const Scalar& s)
    {
        m_expression -= s;
        return *this;
    }

    template<typename OtherExpressionType>
    auto& operator*=(const ScalarWrapper<OtherExpressionType>& other)
    {
        m_expression *= other.m_expression;
        return *this;
    }
    template<typename OtherType>
    auto& operator*=(const OtherType& s)
    {
        m_expression *= s;
        return *this;
    }

    template<typename OtherExpressionType>
    auto& operator/=(const ScalarWrapper<OtherExpressionType>& other)
    {
        m_expression /= other.m_expression;
        return *this;
    }
    auto& operator/=(const Scalar& s)
    {
        m_expression /= s;
        return *this;
    }

    // comparison operators of Scalars must always return a bool, so we need to
    // use matrix for equality and array coefwise compares and .all().
    template<typename OtherExpressionType>
    bool operator==(const ScalarWrapper<OtherExpressionType>& other) const
    { return m_expression.matrix() == other.m_expression.matrix(); }
    template<typename OtherType>
    bool operator==(const OtherType& other) const
    { return m_expression.matrix() == other.matrix(); }
    // bool operator==(const Scalar& s) const { return m_expression == s; }

    template<typename OtherExpressionType>
    bool operator!=(const ScalarWrapper<OtherExpressionType>& other) const
    { return m_expression.matrix() != other.m_expression.matrix(); }
    bool operator!=(const Scalar& s) const { return m_expression.matrix() != s; }

    template<typename OtherExpressionType>
    bool operator<=(const ScalarWrapper<OtherExpressionType>& other) const
    { return (m_expression.array() <= other.m_expression.array()).all(); }
    bool operator<=(const Scalar& s) const { return (m_expression.array() <= s).all(); }

    template<typename OtherExpressionType>
    bool operator<(const ScalarWrapper<OtherExpressionType>& other) const
    { return (m_expression.array() < other.m_expression.array()).all(); }
    bool operator<(const Scalar& s) const { return (m_expression.array() < s).all(); }

    template<typename OtherExpressionType>
    bool operator>=(const ScalarWrapper<OtherExpressionType>& other) const
    { return (m_expression.array() >= other.m_expression.array()).all(); }
    bool operator>=(const Scalar& s) const { return (m_expression.array() >= s).all(); }

    template<typename OtherExpressionType>
    bool operator>(const ScalarWrapper<OtherExpressionType>& other) const
    { return (m_expression.array() > other.m_expression.array()).all(); }
    bool operator>(const Scalar& s) const { return (m_expression.array() > s).all(); }

    const auto format(const IOFormat& fmt) const { return m_expression.format(fmt); }

protected:
    ExpressionType m_expression;
};

/*
// Add the scalar() method to DenseBase for lvalue and rvalue.
template<typename Derived>
EIGEN_DEVICE_FUNC inline ScalarWrapper<Derived&> DenseBase<Derived>::scalar() &
{ return ScalarWrapper<Derived&>(derived()); }
template<typename Derived>
EIGEN_DEVICE_FUNC inline ScalarWrapper<Derived> DenseBase<Derived>::scalar() &&
{ return ScalarWrapper<Derived>(derived()); }
*/

// This will wrap an expression as a scalar with ScalarWrapper. This is an
// alternative to using a CTAD deduction guide for ScalarWrapper().
template<typename ExpressionType>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ScalarWrapper<ExpressionType> scalar(ExpressionType&& m)
{ return ScalarWrapper<ExpressionType>(m); }

// This is a C++17 "CTAD" deduction guide to make ScalarWrapper use a
// reference for an lvalue and a copy for an rvalue.
template<typename ArgXprType>
ScalarWrapper(ArgXprType&& expr) -> ScalarWrapper<ArgXprType>;

template<typename ExpressionType>
auto real(const ScalarWrapper<ExpressionType>& x)
{ return scalar(real(x.m_expression)); }
template<typename ExpressionType>
auto imag(const ScalarWrapper<ExpressionType>& x)
{ return scalar(imag(x.m_expression)); }
template<typename ExpressionType>
auto conj(const ScalarWrapper<ExpressionType>& x)
{ return scalar(conj(x.m_expression)); }
template<typename ExpressionType>
auto sqrt(const ScalarWrapper<ExpressionType>& x)
{ return scalar(sqrt(x.m_expression)); }
template<typename ExpressionType>
auto abs(const ScalarWrapper<ExpressionType>& x)
{ return scalar(abs(x.m_expression)); }
template<typename ExpressionType>
auto cos(const ScalarWrapper<ExpressionType>& x)
{ return scalar(cos(x.m_expression)); }
template<typename ExpressionType>
auto sin(const ScalarWrapper<ExpressionType>& x)
{ return scalar(sin(x.m_expression)); }
template<typename ExpressionType>
auto acos(const ScalarWrapper<ExpressionType>& x)
{ return scalar(acos(x.m_expression)); }
template<typename ExpressionType>
auto atan2(const ScalarWrapper<ExpressionType>& y, const ScalarWrapper<ExpressionType>& x)
{ return scalar(atan2(y.m_expression, x.m_expression)); }

template<typename ExpressionType>
std::ostream& operator<<(std::ostream& stream, const ScalarWrapper<ExpressionType>& x)
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

// Convert a matrix(r,c) into a Vector<SVector> indexed with (r)(c).
template<typename D>
auto AsScalarRows(D&& m)
{
    using CleanD = std::remove_reference_t<D>;
    using VectT  = Eigen::Vector<typename CleanD::Scalar, CleanD::ColsAtCompileTime>;
    using VecSVT = Eigen::Vector<ScalarWrapper<VectT>, CleanD::RowsAtCompileTime>;
    VecSVT ans(m.rows());
    // TODO: Sadly NullararyExpr are non-mutable so we need to create a
    // vector of wrapped rows. It would be good to figure out something
    // better than this.
    for (auto i = 0; i < m.rows(); i++) {
        ans(i) = scalar(m.row(i));
    }
    return ans;
}

// Convert a matrix(r,c) into a Vector<SVector> indexed with (c)(r).
template<typename D>
auto AsScalarCols(D&& m)
{
    using CleanD = std::remove_reference_t<D>;
    using VectT  = Vector<typename CleanD::Scalar, CleanD::RowsAtCompileTime>;
    using VecSVT = Vector<ScalarWrapper<VectT>, CleanD::ColsAtCompileTime>;
    VecSVT ans(m.rows());
    for (auto i = 0; i < m.cols(); i++) {
        ans(i) = scalar(m.col(i));
    }
    return ans;
}

} // namespace Eigen

#endif /* libslic3r_Geometry_ScalarWrapper_hpp_ */
