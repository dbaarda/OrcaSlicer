#ifndef libslic3r_Geometry_EigenNesting_hpp_
#define libslic3r_Geometry_EigenNesting_hpp_

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
template<typename Scalar_, int Rows_, int Options_ = AutoAlign | ColMajor, int MaxRows_ = Rows_>
using Vector = Matrix<Scalar_, Rows_, 1, Options_, MaxRows_, 1>;
template<typename Scalar_, int Cols_, int Options_ = AutoAlign | RowMajor, int MaxCols_ = Cols_>
using RowVector = Matrix<Scalar_, 1, Cols_, Options_, 1, MaxCols_>;

// Define Eigen::ScalarMatrix subclass of Matrix for embedding in other Arrays or Matrixs.
template<typename Matrix_> class ScalarMatrix : public Matrix_
{
public:
    typedef Matrix_ Base;
    // using Base::Scalar;

    ScalarMatrix(void) : Base() {}
    // We need to be able to initialize from a Scalar for Zero(), Ones(), etc.
    explicit ScalarMatrix(const typename Base::Scalar& value) : Base() { Base::setConstant(value); }

    template<typename OtherDerived> ScalarMatrix(const MatrixBase<OtherDerived>& other) : Base(other) {}

    template<typename OtherDerived> ScalarMatrix& operator=(const MatrixBase<OtherDerived>& other)
    {
        this->Base::operator=(other);
        return *this;
    }

    Matrix_& matrix() { return static_cast<Matrix_&>(*this); }
};

template<typename Matrix_> ScalarMatrix<Matrix_>& as_scalar(Matrix_& m) { return static_cast<ScalarMatrix<Matrix_>&>(m); }
template<typename Matrix_> Matrix_&               as_matrix(ScalarMatrix<Matrix_>& m) { return static_cast<Matrix_&>(m); }

// Define NumTraits<ScalarMatrix<> so we can put ScalarMatrix
// types inside other Matrix or Array types. This copies the existing
// NumTraits<Array<...>> definition.
template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
struct NumTraits<ScalarMatrix<Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>>>
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

    // This is the Eigen v5.0.0 version.
    // EIGEN_DEVICE_FUNC constexpr static RealScalar epsilon() { return NumTraits<RealScalar>::epsilon(); }
    // EIGEN_DEVICE_FUNC constexpr static RealScalar dummy_precision() { return NumTraits<RealScalar>::dummy_precision(); }
    // constexpr static int digits10() { return NumTraits<Scalar>::digits10(); }
    // constexpr static int max_digits10() { return NumTraits<Scalar>::max_digits10(); }

    EIGEN_DEVICE_FUNC static inline RealScalar epsilon() { return NumTraits<RealScalar>::epsilon(); }
    EIGEN_DEVICE_FUNC static inline RealScalar dummy_precision() { return NumTraits<RealScalar>::dummy_precision(); }
    static inline int                          digits10() { return NumTraits<_Scalar>::digits10(); }
};

template<typename BinaryOp, typename Matrix_> struct ScalarBinaryOpTraits<ScalarMatrix<Matrix_>, ScalarMatrix<Matrix_>, BinaryOp>
{
    typedef ScalarMatrix<Matrix_> ReturnType;
};

template<typename Matrix_>
struct ScalarBinaryOpTraits<ScalarMatrix<Matrix_>,
                            typename Matrix_::Scalar,
                            internal::scalar_product_op<ScalarMatrix<Matrix_>, typename Matrix_::Scalar>>
    : internal::scalar_product_traits<ScalarMatrix<Matrix_>, typename Matrix_::Scalar>
{
    typedef ScalarMatrix<Matrix_> ReturnType;
};

template<typename Matrix_>
struct ScalarBinaryOpTraits<typename Matrix_::Scalar,
                            ScalarMatrix<Matrix_>,
                            internal::scalar_product_op<typename Matrix_::Scalar, ScalarMatrix<Matrix_>>>
    : internal::scalar_product_traits<typename Matrix_::Scalar, ScalarMatrix<Matrix_>>
{
    typedef ScalarMatrix<Matrix_> ReturnType;
};

/*
template<typename Matrix_, typename OtherMatrix_>
struct ScalarBinaryOpTraits<
       ScalarMatrix<Matrix_>,
       ScalarMatrix<OtherMatrix_>,
       internal::scalar_product_op<
         ScalarMatrix<Matrix_>,
         ScalarMatrix<OtherMatrix_>
       >
> : internal::scalar_product_traits<ScalarMatrix<Matrix_> , ScalarMatrix<OtherMatrix_> >
{
  typedef ScalarMatrix<Product< Matrix_, OtherMatrix_ > > ReturnType;
};
*/

} // namespace Eigen

#endif /* libslic3r_Geometry_EigenNesting_hpp_ */
