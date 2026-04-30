#define NOMINMAX

#include <catch2/catch_all.hpp>
#include <test_utils.hpp>

#include <fstream>

#include <libslic3r/Geometry/ScalarWrapper.hpp>

using namespace Eigen;
using namespace Catch::Matchers;

/*
 * For docs on using catch2;
 *
 * https://catch2-temp.readthedocs.io/en/latest/test-cases-and-sections.html
 */

// Strangely, Eigen doesn't define convenience types for 1x1 matrixes.
typedef Eigen::Matrix<float, 1, 1>  Matrix1f;
typedef Eigen::Matrix<double, 1, 1> Matrix1d;

/******************************************************************************
 * IS_VALID() checks if an eigen expression is valid.
 *
 * The idea here was to use SFINAE to check if an expression was valid and
 * didn't trigger compiler errors. Sadly this doesn't work, with some eigen
 * expressions having valid types with compiler errors on eval(), some having
 * valid eval() with compiler errors on indexing, and most of these errors
 * seem to be non-SFINAE-safe, requiring the checks to be commented out. */
template<typename T, typename = void>
constexpr bool is_valid = false;

template<typename T>
constexpr bool is_valid<T, std::void_t<decltype(std::declval<T>().eval())>> = true;
// constexpr bool is_valid<T, std::void_t<decltype(std::declval<T>()(0))> > = true;
// constexpr bool is_valid<T, std::void_t<typename Eigen::internal::traits<T>::Scalar>> = true;
// constexpr bool is_valid<T, std::void_t<typename T::Scalar>> = true;

#define IS_VALID(expr) is_valid<decltype(expr)>
// #define IS_VALID(expr) Eigen::internal::has_ReturnType<decltype(expr)>::value

/******************************************************************************
 * IS_COND() checks if an expression is a CHECK'able boolean condition. */
template<typename T>
constexpr bool is_condition = std::is_convertible_v<T, bool>;

#define IS_COND(expr) is_condition<decltype(expr)>

/******************************************************************************
 * IS_TRUE() checks if an expression is "true", with special handling for
 * Arrays.
 *
 * Normal condition expressions are just converted into a bool. However, for
 * Arrays, a1 == a2 is not a condition but an array of scalar == scalar
 * results. For Arrays of normal scalars this is an array of bools, but for
 * nested arrays using == would give an array of arrays of bools. However, ==
 * with nested arrays doesn't work because eigen assumes == of scalars always
 * returns a single bool. */
template<typename T>
constexpr std::enable_if_t<is_condition<T>, bool> is_true(const T& v)
{ return bool(v); }

template<typename D>
constexpr bool is_true(const Eigen::ArrayBase<D>& a)
{
    // We apply is_true() to each value to convert nested Arrays into bools first.
    return a.unaryExpr([](auto v) { return is_true(v); }).all();
    // return a.all();
}

#define IS_TRUE(expr) is_true(expr)

/******************************************************************************
 * IS_EQ() checks if two expressions are equal.
 *
 * Because == doesn't return a bool for arrays the idea was to use a smarter
 * is_equ() templated function instead. However, CHECK() nicely displays the
 * matrix contents for == comparisons, so in the end it was better to just
 * use it on arrays and use a0.matrix() == a1.matrix() checks. */
template<typename T1, typename T2>
constexpr std::enable_if_t<is_condition<decltype(std::declval<T1>() == std::declval<T2>())>, bool> is_eq(const T1& v1, const T2& v2)
{ return v1 == v2; }

template<typename D1, typename D2>
constexpr bool is_eq(const Eigen::ArrayBase<D1>& a1, const Eigen::ArrayBase<D2>& a2)
{ return a1.matrix() == a2.matrix(); }

#define IS_EQ(a1, a2) (a1).matrix() == (a2).matrix()

// These check basic Eigen Array and Matrix operations to show what works and
// what doesn't, and to verify that any of our defined eigen extensions
// haven't broken their compilation or execution.
SCENARIO("Using Eigen for basic Array and Vector operations", "[Eigen][BasicArray][BasicVector]")
{
    // Note the commented out checks show various different operations
    // that don't work but we can't check without generating a compiler error.
    GIVEN("A simple Array3f a3f{1,2,3}")
    {
        Array3f a3f(1, 2, 3);
        THEN("array operations work correctly")
        {
            // CHECK_FALSE(IS_VALID(a3f + Array3d(1, 2, 3))); // cannot add arrays with different Scalars.
            CHECK(IS_VALID(a3f + Array3f(1, 2, 3)));       // can add arrays with the same Scalars.
            CHECK_FALSE(IS_COND(a3f == Array3f(0, 0, 0))); // Cannot check array == array because the result isn't a condition.
            CHECK(IS_EQ(a3f, Array3f(1, 2, 3)));           // Can check IS_EQ(array, array) instead.
            CHECK_FALSE(IS_EQ(a3f, Array3f(1, 2, 0)));
            CHECK(IS_TRUE(a3f == Array3f(1, 2, 3))); // Can check IS_TRUE(array == array) also.
            CHECK_FALSE(IS_TRUE(a3f == Array3f(0, 2, 3)));
            CHECK(a3f.matrix() == Vector3f(1, 2, 3)); // Can check array.matrix() == vector.
            CHECK_FALSE(a3f.matrix() == Vector3f(1, 0, 3));
            CHECK(IS_EQ(a3f + 1.0f, Array3f(2, 3, 4)));                   // add array + scalar.
            CHECK(IS_EQ(1 + a3f, Array3f(2, 3, 4)));                      // add scalar + array.
            CHECK(IS_EQ(a3f - 1.0f, Array3f(0, 1, 2)));                   // subtract array - scalar.
            CHECK(IS_EQ(1 - a3f, Array3f(0, -1, -2)));                    // subtract scalar - array.
            CHECK(IS_EQ(2.0f * a3f, Array3f(2, 4, 6)));                   // multiply scalar * array.
            CHECK(IS_EQ(a3f * 2, Array3f(2, 4, 6)));                      // multiply array * scalar.
            CHECK(IS_EQ(2.0f / a3f, Array3f(2.0 / 1, 2.0 / 2, 2.0 / 3))); // divide scalar / array.
            CHECK(IS_EQ(a3f / 2, Array3f(1 / 2.0, 2 / 2.0, 3 / 2.0)));    // divide array / scalar.
            CHECK(IS_EQ(a3f + Array3f(1, 2, 3), Array3f(2, 4, 6)));       // add arrays.
            CHECK(IS_EQ(a3f - Array3f(1, 2, 3), Array3f(0, 0, 0)));       // subtract arrays.
            CHECK(IS_EQ(a3f * Array3f(1, 2, 3), Array3f(1, 4, 9)));       // multiply arrays.
            CHECK(IS_EQ(a3f / Array3f(1, 2, 3), Array3f(1, 1, 1)));       // divide arrays.
        };
    };
    GIVEN("A simple Vector3f v3f{1,2,3}")
    {
        Vector3f v3f(1, 2, 3);
        THEN("vector operations work correctly.")
        {
            // CHECK_FALSE(IS_VALID(v3f + Vector3d(1, 2, 3))); // cannot add vectors with different Scalars.
            CHECK(IS_VALID(v3f + Vector3f(1, 2, 3))); // can add vectors with the same Scalars.
            CHECK(IS_COND(v3f == Vector3f(0, 0, 0))); // Can check vector == vector because the result is a condition.
            CHECK(IS_EQ(v3f, Vector3f(1, 2, 3)));     // Can check IS_EQ(vector, vector).
            CHECK_FALSE(IS_EQ(v3f, Vector3f(1, 2, 0)));
            CHECK(IS_TRUE(v3f == Vector3f(1, 2, 3))); // Can check IS_TRUE(vector == vector) also.
            CHECK_FALSE(IS_TRUE(v3f == Vector3f(0, 2, 3)));
            CHECK(IS_EQ(v3f.array(), Array3f(1, 2, 3))); // Can check IS_EQ(vector.array() == array).
            CHECK_FALSE(IS_EQ(v3f.array(), Array3f(1, 0, 3)));
            // CHECK_FALSE(IS_VALID(v3f + 1.0f));    // add vector + scalar.
            // CHECK_FALSE(IS_VALID(1 + v3f));       // add scalar + vector.
            // CHECK_FALSE(IS_VALID(v3f - 1.0f));    // subtract vector - scalar.
            // CHECK_FALSE(IS_VALID(1 - v3f));     // subtract scalar - vector.
            CHECK(2.0f * v3f == Vector3f(2, 4, 6)); // multiply scalar * vector.
            CHECK(v3f * 2 == Vector3f(2, 4, 6));    // multiply vector * scalar.
            // CHECK_FALSE(IS_VALID(2.0f / v3f)); // divide scalar / vector.
            CHECK(v3f / 2 == Vector3f(1 / 2.0, 2 / 2.0, 3 / 2.0));                // divide vector / scalar.
            CHECK(v3f + Vector3f(1, 2, 3) == Vector3f(2, 4, 6));                  // add vectors.
            CHECK(v3f - Vector3f(1, 2, 3) == Vector3f(0, 0, 0));                  // subtract vectors.
            CHECK(RowVector3f(1, 2, 3) * v3f == Matrix1f(1 * 1 + 2 * 2 + 3 * 3)); // multiply vectors.
            CHECK(RowVector3f(1, 2, 3) * v3f == Matrix1f(1 * 1 + 2 * 2 + 3 * 3));
            CHECK(v3f * RowVector3f(1, 2, 3) == Matrix3f({{1, 2, 3}, {2, 4, 6}, {3, 6, 9}}));
            // CHECK_FALSE(IS_INVALID(v3f*Vector3f(1, 2, 3)));
            // CHECK_FALSE(IS_VALID(v3f / Vector3f(1, 2, 3)));       // divide vectors.
        };
    };
};

// Note defining these ScalarBinaryOpTraits entries is required to make
// operations like `Array<Array<int>> + int` or `Array<Array<int>> *
// Array<int>` work, but if they are defined they break normal operations
// like `Array<int> + Array<int>` because they become ambiguous; is one of
// those arrays actually a Scalar, giving a `Array<Array<int>>` result?
// So we can't define them.

/*
template <typename Scalar_, int Rows_, int Cols_, int Options_, int MaxRows_, int MaxCols_, typename BinaryOp_>
struct ScalarBinaryOpTraits<Array<Scalar_, Rows_, Cols_, Options_, MaxRows_, MaxCols_>, Scalar_, BinaryOp_>
{
    typedef Array<Scalar_, Rows_, Cols_, Options_, MaxRows_, MaxCols_>  ReturnType;
};
template <typename Scalar_, int Rows_, int Cols_, int Options_, int MaxRows_, int MaxCols_, typename BinaryOp_>
struct ScalarBinaryOpTraits<Scalar_, Array<Scalar_, Rows_, Cols_, Options_, MaxRows_, MaxCols_>, BinaryOp_>
{
    typedef Array<Scalar_, Rows_, Cols_, Options_, MaxRows_, MaxCols_>  ReturnType;
};
*/

// Eigen natively supports using Arrays as scalars inside Matrix, Vector, and
// Array, but it has some limitations. We include tests to test this and
// highlight the limitations.
SCENARIO("Using Eigen for nesting Arrays inside an Array", "[Eigen][NestedArray]")
{
    GIVEN("A nested Array3A3f a3a3f{a0,a1,a2}")
    {
        typedef Array<Array3f, 3, 1> Array3A3f;
        Array3f                      a0(1, 2, 3), a1(2, 3, 4), a2(3, 4, 5), a3(4, 5, 6);
        Array3A3f                    a3a3f(a0, a1, a2);
        THEN("Correct operations are valid.")
        {
            CHECK(IS_VALID(a3a3f + a3a3f));
            CHECK(IS_VALID((a3a3f + a3a3f)(0)));
            CHECK(IS_VALID(a3a3f * a3a3f));
            CHECK(IS_VALID((a3a3f * a3a3f)(0)));
            CHECK(IS_VALID(a3a3f / a3a3f));
            CHECK(IS_VALID((a3a3f / a3a3f)(0)));
        }
        THEN("Incorrect operations are invalid")
        {
            CHECK(IS_VALID(a3a3f == a3a3f));      // Strangely we can eval a nested array == operation...
            CHECK(IS_VALID((a3a3f == a3a3f)(0))); // .. but we can't index the result.
            CHECK(IS_VALID(a3a3f != a3a3f));      // Also we can eval a nested array != operation...
            CHECK(IS_VALID((a3a3f != a3a3f)(0))); // .. but we can't index the result.
            // CHECK(IS_VALID((a3a3f + 1.0f)));
            // CHECK(IS_VALID((1 + a3a3f)));
            // CHECK(IS_VALID((a3a3f - 1.0f)));
            // CHECK(IS_VALID((1 - a3a3f)));
            // CHECK(IS_VALID((a3a3f * 2.0f)));
            // CHECK(IS_VALID((2 * a3a3f)));
            // CHECK(IS_VALID((a3a3f / 2.0f)));
            // CHECK(IS_VALID((2 / a3a3f)));
        }
        WHEN("adding a3a3f + a3a3f")
        {
            auto ans = a3a3f + a3a3f;
            THEN("the result is correct")
            {
                // CHECK(IS_TRUE(a3a3f + a3a3f == 2 * a3a3f));
                // CHECK(IS_EQ(a3a3f + a3a3f, Array3A3f(2 * a0, 2 * a1, 2 * a2)));
                CHECK(IS_EQ(ans(0), 2 * a0));
                CHECK(IS_EQ(ans(1), 2 * a1));
                CHECK(IS_EQ(ans(2), 2 * a2));
            };
        };
        WHEN("subtracting a3a3f - a3a3f")
        {
            auto ans  = a3a3f - a3a3f;
            auto zero = Array3A3f::Zero();
            THEN("the result is all zeros")
            {
                CHECK(ans(0).matrix() == zero(0).matrix());
                CHECK(ans(1).matrix() == zero(1).matrix());
                CHECK(ans(2).matrix() == zero(2).matrix());
            };
        };
        WHEN("multiplying a3a3f * a3a3f")
        {
            auto ans = a3a3f * a3a3f;
            THEN("the result is correct")
            {
                CHECK(ans(0).matrix() == (a0 * a0).matrix());
                CHECK(ans(1).matrix() == (a1 * a1).matrix());
                CHECK(ans(2).matrix() == (a2 * a2).matrix());
            };
        };
        WHEN("dividing a3a3f / a3a3f.")
        {
            auto ans = a3a3f / a3a3f;
            THEN("the result is correct")
            {
                CHECK(ans(0).matrix() == (a0 / a0).matrix());
                CHECK(ans(1).matrix() == (a1 / a1).matrix());
                CHECK(ans(2).matrix() == (a2 / a2).matrix());
            };
        };
        WHEN("adding a3a3f + a3")
        {
            auto ans = a3a3f + a3;
            THEN("the result is corrects.")
            {
                CHECK(ans(0).matrix() == (a0 + a3).matrix());
                CHECK(ans(1).matrix() == (a1 + a3).matrix());
                CHECK(ans(2).matrix() == (a2 + a3).matrix());
            };
        };
        WHEN("adding a3 + a3a3f")
        {
            auto ans = a3 + a3a3f;
            THEN("the result is correct")
            {
                CHECK(ans(0).matrix() == (a3 + a0).matrix());
                CHECK(ans(1).matrix() == (a3 + a1).matrix());
                CHECK(ans(2).matrix() == (a3 + a2).matrix());
            };
        };
        WHEN("multiplying a3a3f * a3")
        {
            auto ans = a3a3f * a3;
            THEN("the result is correct")
            {
                CHECK(ans(0).matrix() == (a0 * a3).matrix());
                CHECK(ans(1).matrix() == (a1 * a3).matrix());
                CHECK(ans(2).matrix() == (a2 * a3).matrix());
            };
        };
        WHEN("multiplying a3 * a3a3f")
        {
            auto ans = a3 * a3a3f;
            THEN("the result is correct")
            {
                CHECK(ans(0).matrix() == (a3 * a0).matrix());
                CHECK(ans(1).matrix() == (a3 * a1).matrix());
                CHECK(ans(2).matrix() == (a3 * a2).matrix());
            };
        };
        /*
        // The following doesn't work because a3a3f + 7.0f doesn't compile
        // without the additional ScalarBinaryOpTraits<> defined.
        WHEN("adding a3a3f + 7.0f")
        {
            auto ans = a3a3f + 7.0f;
            THEN("Check a3a3f + float results.")
            {
                CHECK(ans(0).matrix() == (a0 + 7.0f).matrix());
                CHECK(ans(1).matrix() == (a1 + 7.0f).matrix());
                CHECK(ans(2).matrix() == (a2 + 7.0f).matrix());
            };
        };
        */
    };
};

SCENARIO("Using Eigen's Nesting an Array inside a Matrix", "[Eigen][NestedArray]")
{
    GIVEN("A Nested Array Matrix2A3f m2a3f{{a0,a1},{a2,a3}}")
    {
        typedef Matrix<Array3f, 2, 2> Matrix2A3f;
        Array3f                       a0(0, 0, 0), a1(0, 1, 1), a2(1, 0, 2), a3(1, 1, 3);
        Matrix2A3f                    m2a3f{{a0, a1}, {a2, a3}};
        Matrix2f                      m2f{{0, 1}, {2, 3}};
        WHEN("adding m2a3f + m2a3f")
        {
            auto ans = m2a3f + m2a3f;
            THEN("the result is 2 * m2a3f")
            {
                CHECK(ans(0, 0).matrix() == (2 * a0).matrix());
                CHECK(ans(0, 1).matrix() == (2 * a1).matrix());
                CHECK(ans(1, 0).matrix() == (2 * a2).matrix());
                CHECK(ans(1, 1).matrix() == (2 * a3).matrix());
            };
        };
        WHEN("subtracting m2a3f - m2a3f")
        {
            auto ans  = m2a3f - m2a3f;
            auto zero = Matrix2A3f::Zero();
            THEN("the result is Matrix2A3f::Zero().")
            {
                CHECK(ans(0, 0).matrix() == zero(0, 0).matrix());
                CHECK(ans(0, 1).matrix() == zero(0, 1).matrix());
                CHECK(ans(1, 0).matrix() == zero(1, 0).matrix());
                CHECK(ans(1, 1).matrix() == zero(1, 1).matrix());
            };
        };
        WHEN("multiplying m2a3f * m2a3f.")
        {
            auto ans = m2a3f * m2a3f;
            THEN("the result is correct.")
            {
                CHECK(ans(0, 0).matrix() == (a0 * a0 + a1 * a2).matrix());
                CHECK(ans(0, 1).matrix() == (a0 * a1 + a1 * a3).matrix());
                CHECK(ans(1, 0).matrix() == (a2 * a0 + a3 * a2).matrix());
                CHECK(ans(1, 1).matrix() == (a2 * a1 + a3 * a3).matrix());
            };
        };
        // Matrix2<Array3f> / Matrix2<Array3f> doesn't work because you can't divide matrixes.
        // Matrix2<Array3f> + Array3f doesn't work because you can't add a matrix to a scalar.
        // Array3f + Matrix2<Array3f> doesn't work because you can't add a scalar to a matrix.
        WHEN("multiplying m2a3f * a3")
        {
            auto ans = m2a3f * a3;
            THEN("the result is correct")
            {
                CHECK(ans(0, 0).matrix() == (a0 * a3).matrix());
                CHECK(ans(0, 1).matrix() == (a1 * a3).matrix());
                CHECK(ans(1, 0).matrix() == (a2 * a3).matrix());
                CHECK(ans(1, 1).matrix() == (a3 * a3).matrix());
            };
        };
        WHEN("multiplying a3 * m2a3f")
        {
            auto ans = a3 * m2a3f;
            THEN("the result is correct")
            {
                CHECK(ans(0, 0).matrix() == (a3 * a0).matrix());
                CHECK(ans(0, 1).matrix() == (a3 * a1).matrix());
                CHECK(ans(1, 0).matrix() == (a3 * a2).matrix());
                CHECK(ans(1, 1).matrix() == (a3 * a3).matrix());
            };
        };
        WHEN("dividing m2a3f / a3")
        {
            auto ans = m2a3f / a3;
            THEN("the result is correct")
            {
                CHECK(ans(0, 0).matrix() == (a0 / a3).matrix());
                CHECK(ans(0, 1).matrix() == (a1 / a3).matrix());
                CHECK(ans(1, 0).matrix() == (a2 / a3).matrix());
                CHECK(ans(1, 1).matrix() == (a3 / a3).matrix());
            };
        };
        // Array3f / Matrix2<Array3f> doesn't work because you can't divide a scalar by a matrix.
        /* This doesn't work because m2a3f * 7.0f doesn't compile, even though Array3f * 0.7f works.
        WHEN("Matrix2<Array3f> is multiplied by a float.")
        {
            auto ans = m2a3f * 7.0f;
            THEN("Check m2a3f * 7.0f results.")
            {
                CHECK(ans(0, 0).matrix() == (a0 * 7.0f).matrix());
                CHECK(ans(0, 1).matrix() == (a1 * 7.0f).matrix());
                CHECK(ans(1, 0).matrix() == (a2 * 7.0f).matrix());
                CHECK(ans(1, 1).matrix() == (a3 * 7.0f).matrix());
            };
        };
        */
        /* This doesn't work because m2 * m2a3f doesn't compile, even though float*Array3f works.
        WHEN("Matrix2f is multiplied by Matrix2<Array3f>.")
        {
            auto ans = m2f * m2a3f;
            THEN("Check m2f * m2a3f results.")
            {
                CHECK(ans(0, 0).matrix() == (m2f(0, 0) * a0 + m2f(0, 1) * a2).matrix());
                CHECK(ans(0, 1).matrix() == (m2f(0, 0) * a1 + m2f(0, 1) * a3).matrix());
                CHECK(ans(1, 0).matrix() == (m2f(1, 0) * a0 + m2f(1, 1) * a2).matrix());
                CHECK(ans(1, 1).matrix() == (m2f(1, 0) * a1 + m2f(1, 1) * a3).matrix());
            };
        };
        */
    };
};

SCENARIO("ScalarWrapper can wrap arrays and matrixs", "[ScalarWrapper]")
{
    GIVEN("Values v1=Vector3f(2,3,4) and a1=Array3f(3,4,5)")
    {
        typedef ScalarWrapper<Vector3f> SVector3f;
        typedef ScalarWrapper<Array3f>  SArray3f;
        Vector3f                        v1(2, 3, 4);
        Array3f                         a1(3, 4, 5);
        WHEN("sv1 and sa1 are auto-typed wrappers around lvalues v1 and a1")
        {
            auto sv1 = ScalarWrapper(v1);
            auto sa1 = ScalarWrapper(a1);
            THEN("derived(), matrix(), and array() return the right values")
            {
                CHECK(sv1.derived() == Vector3f(2, 3, 4));
                CHECK((sa1.derived() == Array3f(3, 4, 5)).all());
                CHECK(sv1.matrix() == Vector3f(2, 3, 4));
                CHECK(sa1.matrix() == Vector3f(3, 4, 5));
                CHECK((sa1.array() == Array3f(3, 4, 5)).all());
                CHECK((sv1.array() == Array3f(2, 3, 4)).all());
            };
            AND_WHEN("mutating sv1 derived(), matrix(), and array() elements with +=1")
            {
                sv1.derived()(0) += 1;
                sv1.matrix()(1) += 1;
                sv1.array()(2) += 1;
                THEN("sv1 elements have been incremented")
                {
                    CHECK(sv1.derived() == Vector3f(3, 4, 5));
                    AND_THEN("v1 elements have been incremented") { CHECK(v1 == Vector3f(3, 4, 5)); };
                };
            };
            AND_WHEN("mutating with sv1 *= 2")
            {
                sv1 *= 2;
                THEN("v1 is mutated")
                {
                    CHECK(sv1.derived() == v1);
                    CHECK(v1 == Vector3f(4, 6, 8));
                };
            };
            AND_WHEN("mutating with v1*=2")
            {
                v1 *= 2;
                THEN("sv1 has been doubled")
                {
                    CHECK(sv1.derived() == v1);
                    CHECK(sv1.derived() == Vector3f(4, 6, 8));
                };
            };
            AND_WHEN("sv1 is set to another sv2 wrapper around lvalue v2={5,6,7}")
            {
                Vector3f v2(5, 6, 7);
                auto     sv2 = ScalarWrapper(v2);
                sv1          = sv2;
                THEN("v1 now equals v2")
                {
                    CHECK(v1 == v2);
                    CHECK(v1 == Vector3f(5, 6, 7));
                };
                AND_WHEN("mutating with v2*=2")
                {
                    v2 *= 2;
                    THEN("sv2 is doubled")
                    {
                        CHECK(sv2.derived() == v2);
                        CHECK(sv2.derived() == Vector3f(10, 12, 14));
                        AND_THEN("v1 is not changed like v2")
                        {
                            CHECK(v1 != v2);
                            CHECK(v1 == Vector3f(5, 6, 7));
                        };
                    };
                };
            };
        };
        WHEN("Intializing with auto wrappers around rvalues sv1=ScalarWrapper(Vector3f(2, 3, 4)) and sa1=ScalarWrapper(Array3f(3, 4, 5))")
        {
            auto sv1 = ScalarWrapper(Vector3f(2, 3, 4));
            auto sa1 = ScalarWrapper(Array3f(3, 4, 5));
            THEN("sv1 and sa1 have the right values")
            {
                CHECK(sv1.derived() == Vector3f(2, 3, 4));
                CHECK((sa1.derived() == Array3f(3, 4, 5)).all());
            };
            AND_WHEN("mutating sv1 derived(), matrix(), and array() elements with +=1")
            {
                sv1.derived()(0) += 1;
                sv1.matrix()(1) += 1;
                sv1.array()(2) += 1;
                THEN("sv1 elements have been incremented") { CHECK(sv1.derived() == Vector3f(3, 4, 5)); };
            };
            AND_WHEN("mutating with sv1 *= 2")
            {
                sv1 *= 2;
                THEN("sv1 is doubled") { CHECK(sv1.derived() == Vector3f(4, 6, 8)); };
            };
        };
        WHEN("Initializing with explicit-typed lvalues initialized with sv1{2,3,4} and sa1{3,4,5}")
        {
            SVector3f sv1{2, 3, 4};
            SArray3f  sa1{3, 4, 5};
            THEN("sv1 and sa1 have the right values")
            {
                CHECK(sv1.derived() == Vector3f(2, 3, 4));
                CHECK((sa1.derived() == Array3f(3, 4, 5)).all());
            };
            AND_WHEN("mutating sv1 derived(), matrix(), and array() elements with +=1")
            {
                sv1.derived()(0) += 1;
                sv1.matrix()(1) += 1;
                sv1.array()(2) += 1;
                THEN("sv1 elements have been incremented") { CHECK(sv1.derived() == Vector3f(3, 4, 5)); };
            };
            AND_WHEN("mutating with sv1 *= 2")
            {
                sv1 *= 2;
                THEN("sv1 is doubled") { CHECK(sv1.derived() == Vector3f(4, 6, 8)); };
            };
        };
    };
};

SCENARIO("Test ScalarWrapper<Vector3f> inside a Matrix", "[ScalarWrapper]")
{
    GIVEN("Classes are defined and instances are initialized.")
    {
        typedef ScalarWrapper<Vector3f> ScalarV3f;
        typedef Matrix<ScalarV3f, 2, 2> Matrix2V3f;
        Vector3f                        v0(0, 0, 0), v1(2, 3, 4), v2(3, 4, 5), v3(4, 5, 6);
        // Matrix2V3f                      m2v3f{{v0, v1}, {v2, v3}};
        // Matrix2V3f                      m2v3f({{v0, v1}, {v2, v3}});
        // Matrix2V3f                      m2v3f({v0, v1, v2, v3});
        // Matrix2V3f                      m2v3f(v0, v1, v2, v3);
        Matrix2V3f m2v3f{{scalar(v0), scalar(v1)}, {scalar(v2), scalar(v3)}};
        m2v3f(0, 0) = v0;
        m2v3f(0, 1) = v1;
        m2v3f(1, 0) = v2;
        m2v3f(1, 1) = v3;
        Matrix2f m2f{{0, 1}, {2, 3}};
        // Matrix2V3f m1{{v0,v1},{v2,v3}};
        THEN("Check scalar() and matrix() operations.")
        {
            CHECK(v0 == Vector3f(0, 0, 0));
            CHECK(scalar(v0).matrix() == Vector3f(0, 0, 0));
            // CHECK(m2v3f(0, 0) == v0); need to uwrap before we can compare.
            CHECK(m2v3f(0, 0).matrix() == v0);
            CHECK(m2v3f(0, 1).matrix() == v1);
            CHECK(m2v3f(1, 0).matrix() == v2);
            CHECK(m2v3f(1, 1).matrix() == v3);
        };
        WHEN("Matrix2<ScalarWrapper<Vector3f>> instances are added and subtracted.")
        {
            THEN("Check m2v3f additions, subtractions, and N * m2v3f results.")
            {
                CHECK(m2v3f + m2v3f == 2.0f * m2v3f);
                CHECK(m2v3f + m2v3f == m2v3f * 2);
                CHECK(m2v3f + m2v3f + m2v3f == 3.0f * m2v3f);
                CHECK(2.0f * m2v3f + m2v3f == m2v3f * 3.0f);
                CHECK(m2v3f - m2v3f == Matrix2V3f::Zero());
            };
            THEN("Check m2v3f.array() + scalar(v3f) results.")
            {
                CHECK(((m2v3f.array() + scalar(v1)) == (scalar(v1) + m2v3f.array())).all());
                CHECK((m2v3f.array() + scalar(v0)).matrix() == m2v3f);
                CHECK((m2v3f.array() + scalar(v1)).matrix() ==
                      Matrix2V3f({{scalar(v0 + v1), scalar(v1 + v1)}, {scalar(v2 + v1), scalar(v3 + v1)}}));
                CHECK((scalar(v2) + m2v3f.array()).matrix() ==
                      Matrix2V3f({{scalar(v2 + v0), scalar(v2 + v1)}, {scalar(v2 + v2), scalar(v2 + v3)}}));
                CHECK((m2v3f.array() + scalar(v1))(0, 0).matrix() == v0 + v1);
                CHECK((m2v3f.array() + scalar(v1))(0, 1).matrix() == v1 + v1);
                CHECK((m2v3f.array() + scalar(v1))(1, 0).matrix() == v2 + v1);
                CHECK((m2v3f.array() + scalar(v1))(1, 1).matrix() == v3 + v1);
            };
        };
        WHEN("Instances multiplied.")
        {
            THEN("Check m2f * scalar(v3f) gives m2v3f.")
            {
                CHECK(Matrix2V3f::Zero() == 0.0f * m2v3f);
                CHECK((m2f * scalar(v1)) == (scalar(v1) * m2f));
                CHECK((m2f * scalar(v0)) == Matrix2V3f::Zero());
                // CHECK(m2f * scalar(v1) == Matrix2V3f({{v1*0.0f, v1*1.0f}, {v1*2.0f, v1*3.0f}}));
                CHECK(m2f * scalar(v1) == Matrix2V3f({{scalar(v1 * 0.0f), scalar(v1 * 1.0f)}, {scalar(v1 * 2.0f), scalar(v1 * 3.0f)}}));
                CHECK((m2f * scalar(v1))(0, 0).matrix() == v0);
                CHECK((m2f * scalar(v1))(0, 1).matrix() == v1);
                CHECK((m2f * scalar(v1))(1, 0).matrix() == 2 * v1);
                CHECK((m2f * scalar(v1))(1, 1).matrix() == 3 * v1);
            };
            THEN("Check m2v3f * N gives m2v3f.")
            {
                CHECK((m2v3f * 2.0f) == (2.0f * m2v3f));
                CHECK((m2v3f * 0.0f) == Matrix2V3f::Zero());
                CHECK((m2v3f * 2.0f)(0, 0).matrix() == 2 * v0);
                CHECK((m2v3f * 2.0f)(0, 1).matrix() == 2 * v1);
                CHECK((m2v3f * 2.0f)(1, 0).matrix() == 2 * v2);
                CHECK((m2v3f * 2.0f)(1, 1).matrix() == 2 * v3);
            };
            THEN("Check m2f * m2v3f gives m2v3f.")
            {
                CHECK((m2f * m2v3f) != (m2v3f * m2f));
                CHECK((m2f * m2v3f)(0, 0) == m2f(0, 0) * m2v3f(0, 0) + m2f(0, 1) * m2v3f(1, 0));
                CHECK((m2f * m2v3f)(0, 1) == m2f(0, 0) * m2v3f(0, 1) + m2f(0, 1) * m2v3f(1, 1));
                CHECK((m2f * m2v3f)(1, 0) == m2f(1, 0) * m2v3f(0, 0) + m2f(1, 1) * m2v3f(1, 0));
                CHECK((m2f * m2v3f)(1, 1) == m2f(1, 0) * m2v3f(0, 1) + m2f(1, 1) * m2v3f(1, 1));
            };
        };
    };
}

// The following is to setup benchmarking of different ways to do a tri-cubic
// interpolation of a 3D field.

// Get the compile-time size of a VXY field.
constexpr int N_VXY(int N) { return (N == Eigen::Dynamic) ? Eigen::Dynamic : N * N; }
// Get the compile-time size of a VXYZ field.
constexpr int N_VXYZ(int N) { return (N == Eigen::Dynamic) ? Eigen::Dynamic : N * N * N; }

// 1D Field implementation types.
template<typename Scalar, int N = Eigen::Dynamic>
using Field1D_VX = Eigen::Vector<Scalar, N>;

template<typename Scalar, int N = Eigen::Dynamic, int S = Eigen::Dynamic>
using MapField1D_VX = Eigen::Map<Field1D_VX<Scalar, S>>;

// 2D Field implementation types.
template<typename Scalar, int N = Eigen::Dynamic>
using Field2D_VXY = Field1D_VX<Scalar, N_VXY(N)>;

template<typename Scalar, int N = Eigen::Dynamic>
using Field2D_MXY = Eigen::Matrix<Scalar, N, N>;

template<typename Scalar, int N = Eigen::Dynamic, int S = Eigen::Dynamic>
using MapField2D_MXY = Eigen::Map<Field2D_MXY<Scalar, S>, 0, Eigen::OuterStride<N>>;

template<typename Scalar, int N = Eigen::Dynamic>
using Field2D_VYVX = Field1D_VX<ScalarWrapper<Field1D_VX<Scalar, N>>, N>;

template<typename Scalar, int N = Eigen::Dynamic>
using MapField2D_VYVX = Field1D_VX<ScalarWrapper<MapField1D_VX<Scalar, N>>, N>;

// 3D Field implementation types.
template<typename Scalar, int N = Eigen::Dynamic>
using Field3D_VXYZ = Eigen::Vector<Scalar, N_VXYZ(N)>;

template<typename Scalar, int N = Eigen::Dynamic> // indexed with (x+y*n,z)
using Field3D_MXYZ = Eigen::Matrix<Scalar, N_VXY(N), N>;

template<typename Scalar, int N = Eigen::Dynamic> // indexed with (x,y+z*n)
using Field3D_MX_YZ = Eigen::Matrix<Scalar, N, N_VXY(N)>;

// TODO: this needs a NullaryExpr to work properly?
// template<typename Scalar, int N=Eigen::Dynamic, int S=Eigen::Dynamic>
// using MapField3D_MXY_Z = Eigen::Map<Field3D_MXYZ<Scalar, S>, 0, Eigen::OuterStride<N_VXY(N)>>;

template<typename Scalar, int N = Eigen::Dynamic>
using Field3D_VZMXY = Field1D_VX<ScalarWrapper<Field2D_MXY<Scalar, N>>, N>;

template<typename Scalar, int N = Eigen::Dynamic, int S = Eigen::Dynamic>
using MapField3D_VZMXY = Eigen::Vector<ScalarWrapper<MapField2D_MXY<Scalar, N, S>>, S>;

template<typename Scalar, int N = Eigen::Dynamic>
using Field3D_VZVYVX = Eigen::Vector<ScalarWrapper<Field2D_VYVX<Scalar, N>>, N>;

template<typename Scalar, int N = Eigen::Dynamic>
using MapField3D_VZVYVX = Eigen::Vector<ScalarWrapper<MapField2D_VYVX<Scalar, N>>, N>;

// Get the index for (x,y) in a size n VXY field.
constexpr int idx_VXY(int x, int y, int n) { return x + y * n; }
// Get the index for index xy in a size s VXY view at (x,y) of a size n VXY field.
constexpr int blk_VXY(int xy, int x, int y, int n, int s) { return idx_VXY(x + xy % s, y + xy / s, n); }
// Get the index for (x,y,z) in a size n VXYZ field.
constexpr int idx_VXYZ(int x, int y, int z, int n) { return x + (y + z * n) * n; }
// Get the index for index xyz in a size s VXYZ view at (x,y,z) of a size n VXYZ field.
constexpr int blk_VXYZ(int xyz, int x, int y, int z, int n, int s)
{ return idx_VXYZ(x + xyz % s, y + (xyz / s) % s, z + xyz / (s * s), n); }
// Indexer for size s VXY views at (x,y) in a size n VXY field.
struct blkidx_VXY
{
    const int _x, _y, _n, _s;
    constexpr blkidx_VXY(int x, int y, int n, int s) : _x{x}, _y{y}, _n{n}, _s{s} {}
    constexpr int size() const { return _s * _s; }
    constexpr int operator[](int xy) const { return blk_VXY(xy, _x, _y, _n, _s); }
};

// The following are for getting size 4 1D, 4x4 2D, and 4x4x4 3D blocks from
// 3D fields.
template<typename Scalar, int N = Eigen::Dynamic>
auto get4(Field3D_MXYZ<Scalar, N>& f, Eigen::Index x, Eigen::Index y, Eigen::Index z)
{
    int n = f.cols();
    return f.template block<4, 1>(idx_VXY(x, y, n), z);
    // return f(seqN(idx_VXY(x,y,n), fix<4>), z);
    // return MapField1D_VX<Scalar,4>(f.data() + idx_VXY(x,y,n));
    // return Field1D_VX::NullaryExpr(4,[&f, xyo=idx_VXY(x,y,n), zo=z](Index x) -> Scalar& { return f(xyo+x, zo); });
}

template<typename Scalar, int N = Eigen::Dynamic>
auto get44(Field3D_MXYZ<Scalar, N>& f, const Eigen::Index x, const Eigen::Index y, const Eigen::Index z)
{
    int n = f.cols();
    assert(f.rows() == n * n);
    return f.reshaped(n, n * n).template block<4, 4>(x, idx_VXY(y, z, n));
    // return f.reshaped(n, n*n)(seqN(x,fix<4>),seqN(idx_VXY(y,z,n), fix<4>));
    // return MapField2D_MXY<Scalar, N, 4>(f.data() + idx_VXYZ(x,y,z,n), 4, 4, OuterStride<N>(n));
    // return Field2D_MXY::NullaryExpr(4,4,[&f, xo=x, yo=y, zo=z](Index x, Index y) -> Scalar& { return f(idx_VXY(xo+x,yo+y), zo); });
}

template<typename Scalar, int N = Eigen::Dynamic>
auto get444(Field3D_MXYZ<Scalar, N>& f, Eigen::Index x, Eigen::Index y, Eigen::Index z)
{
    int n = f.cols();
    assert(f.rows() == n * n);
    auto a1 = (f.template middleCols<4>(z)).template reshaped<RowMajor>(n, 4 * n);          // (XY,Z), select Z, -> (Y, ZX)
    auto a2 = a1.template middleRows<4>(y).template reshaped<ColMajor>(fix<4 * 4>, n);      // (Y,ZX), select Y, -> (YZ, X)
    auto a3 = a2.template middleCols<4>(x).template reshaped<RowMajor>(fix<4>, fix<4 * 4>); // (YZ,X), select X, -> (Z, XY)
    return a3.transpose().template reshaped<ColMajor>(fix<4>, fix<4 * 4>);                  // (Z,XY) -> (XY,Z) -> (X, YZ)

    // return f(blkidx_VXY(x,y,n,4), seqN(z, fix<4>));
    // return Field3D_MXYZ<Scalar,4>::NullaryExpr(4*4,4,[&f, xo=x, yo=y, zo=z](Index xy, Index z){ return f(win_VXY(xy,xo,yo,n,4),zo); });
}

template<typename Scalar, int N = Eigen::Dynamic>
auto get44(Field3D_VZMXY<Scalar, N>& f, Eigen::Index x, Eigen::Index y, Eigen::Index z)
{
    int n = f.size();
    return f(z).derived().template block<4, 4>(x, y);
    // return f(z).derived()(seqN(x, fix<4>), seqN(y, fix<4>));
    // return MapField2D_MXY<Scalar, 4, N>(f(z).derived().data() + x + y*n, 4, 4, OuterStride<N>(n));
    // return Field2D_MXY::NullaryExpr(4,4,[&f, xo=x, yo=y, zo=z](Index x, Index y) -> Scalar& { return f(zo)(xo+x,yo+y); });
}

template<typename Scalar, int N = Eigen::Dynamic>
auto& get44(Field3D_VZVYVX<Scalar, N>& f, Eigen::Index x, Eigen::Index y, Eigen::Index z)
{
    return Field2D_VYVX<Scalar, 4>::NullaryExpr(
        [&f, xo = x, yo = y, zo = z](int y) { return ScalarWrapper(f(zo)(yo + y).derived().segment<4>(xo)); });
    // return Field2D_MXY::NullaryExpr(4,4,[&f, xo=x, yo=y, zo=z](Index x, Index y) -> Scalar& { return f(zo).derived()(yo+y),derived()(xo+x); });
}

template<typename Scalar>
struct cubic
{
    // This is the matrix for storing (multiple sets of) the 4 points to interoplate with.
    template<int Cols = Dynamic, int MaxCols = Cols>
    using FMatrix = Matrix<Scalar, 4, Cols, 0, 4, MaxCols>;
    // This is the matrix for storing(multiple sets of) the 4 cubic interpolation coefficients.
    template<int Rows = Dynamic, int MaxRows = Rows>
    using CMatrix = Matrix<Scalar, Rows, 4, 0, MaxRows, 4>;
    // This is the matrix for storing the intermediate and final interpolation results.
    template<typename FSource>
    using AMatrix = Matrix<typename FSource::Scalar, Dynamic, Dynamic, 0, 4, std::min(FSource::MaxColsAtCompileTime / 4, Dynamic)>;

    template<typename C, typename F>
    inline static auto cubicN(const Eigen::MatrixBase<C>& c, const Eigen::MatrixBase<F>& f)
    {
        assert(c.cols() == 4 && c.rows() >= 1);
        assert(f.rows() == 4 && f.cols() % (1 << (2 * (c.rows() - 1))) == 0);
        AMatrix<F> fn = c.row(0) * f;
        for (int r = 1; r < c.rows(); r++) {
            fn = c.row(r) * fn.reshaped(fix<4>, AutoSize);
        }
        return fn;
    }
};

template<typename Scalar, int N = Eigen::Dynamic>
auto getF1(Field3D_MXYZ<Scalar, N>& f, Eigen::Index x, Eigen::Index y, Eigen::Index z)
{ return get4(f, x, y, z); }

template<typename Scalar, int N = Eigen::Dynamic>
auto getF2(Field3D_MXYZ<Scalar, N>& f, Eigen::Index x, Eigen::Index y, Eigen::Index z)
{ return get44(f, x, y, z); }

template<typename Scalar, int N = Eigen::Dynamic>
auto getF3(Field3D_MXYZ<Scalar, N>& f, Eigen::Index x, Eigen::Index y, Eigen::Index z)
{
    typename cubic<Scalar>::template FMatrix<4 * 4> ans;
    for (auto i = 0; i < 4; i++)
        ans.template block<4, 4>(0, i * 4) = get44(f, x, y, z + i);
    // ans << get44(f, x, y, z), get44(f, x, y, z + 1), get44(f, x, y, z + 2), get44(f, x, y, z + 3);
    // assert(ans.rows() == 4);
    // assert(ans.cols() == 4 * 4);
    return ans;
}

template<typename Scalar, int N = Eigen::Dynamic>
auto getF3stdvec(Field3D_MXYZ<Scalar, N>& f, Eigen::Index x, Eigen::Index y, Eigen::Index z)
{
    typename std::vector<typename cubic<Scalar>::template FMatrix<4>> fs{get44(f, x, y, z), get44(f, x, y, z + 1), get44(f, x, y, z + 2),
                                                                         get44(f, x, y, z + 3)};
    return fs;
}

template<typename C, typename F>
inline static auto cubic3stdvec(const Eigen::MatrixBase<C>& c, const F& f)
{
    using Scalar = typename C::Scalar;

    assert(c.cols() == 4 && c.rows() == 3);
    Vector<Scalar, 4> fz;
    auto              c2 = c.topRows(2);
    fz(0)                = cubic<Scalar>::cubicN(c2, f[0])(0);
    fz(1)                = cubic<Scalar>::cubicN(c2, f[1])(0);
    fz(2)                = cubic<Scalar>::cubicN(c2, f[2])(0);
    fz(3)                = cubic<Scalar>::cubicN(c2, f[3])(0);
    return cubic<Scalar>::cubicN(c.row(2), fz);
}

// Helper functions for getting random x,y,z in range 0-(n-3)] for iteration i.
int          lcg(int i) { return int((1103515245u * uint32_t(i) + 12345u) & 0x7fffffff); }
Eigen::Index getx(int i, int n) { return lcg(i) % (n - 3); }
Eigen::Index gety(int i, int n) { return lcg(i + 1) % (n - 3); }
Eigen::Index getz(int i, int n) { return lcg(i + 2) % (n - 3); }

// Helper function for building SECTION/BENCHMARK/ETC titles.
template<typename... Args>
std::string title(Args&&... args)
{
    std::stringstream ss;
    (ss << ... << std::forward<Args>(args));
    return ss.str();
}

SCENARIO("Benchmark ScalarWrapper", "[ScalarWrapper][!benchmark]")
{
    GIVEN("Classes are defined and instances are initialized.")
    {
        constexpr Eigen::Index N = 300;
        Eigen::Index           n = N;
        // These are the interpolation coefficients.
        THEN("benchmark reshape()")
        {
            Field3D_MXYZ<float> f3f;
            f3f.setRandom(n * n, n);
            Matrix4f A = Matrix4f::Random();
            BENCHMARK(title("reshape() to shift axies and index for block3 with incremental expressions"), i)
            {
                auto x  = getx(i, n);
                auto y  = gety(i, n);
                auto z  = getz(i, n);
                auto a1 = f3f.template reshaped<ColMajor>(n, n * n);
                auto a2 = a1(seqN(x, fix<4>), blkidx_VXY(y, z, n, 4));
                auto a3 = A.row(0) * a2;
                auto a4 = a3.template reshaped<ColMajor>(fix<4>, fix<4>);
                auto a5 = A.row(1) * a4;
                auto a6 = a5.reshaped();
                auto a7 = A.row(2) * a6;
                return a7.eval();
            };
            BENCHMARK(title("reshape() to shift axies and index for block3 as a single line expression"), i)
            {
                auto x = getx(i, n);
                auto y = gety(i, n);
                auto z = getz(i, n);
                return (A.row(2) *
                        ((A.row(1) * ((A.row(0) * ((f3f.template reshaped<ColMajor>(n, n * n))(seqN(x, fix<4>), blkidx_VXY(y, z, n, 4))))
                                          .template reshaped<ColMajor>(fix<4>, fix<4>)))
                             .reshaped()))
                    .eval();
            };
            BENCHMARK(title("Map to shift axies and index for block3"), i)
            {
                auto x = getx(i, n);
                auto y = gety(i, n);
                auto z = getz(i, n);
                return (A.row(2) * ((A.row(1) * ((A.row(0) * (MatrixXf::Map(f3f.data(), n, n * n)(seqN(x, fix<4>), blkidx_VXY(y, z, n, 4))))
                                                     .template reshaped<ColMajor>(fix<4>, fix<4>))))
                                       .reshaped())
                    .eval();
            };
        };
        THEN("benchmark cubicN()")
        {
            std::stringstream         suffix;
            cubic<double>::CMatrix<3> c3d;
            cubic<float>::CMatrix<3>  c3f;
            c3d.setRandom(3, 4);
            c3f = c3d.cast<float>();
            Field3D_MXYZ<double> f3d;
            Field3D_MXYZ<float>  f3f;
            f3d.setRandom(n * n, n);
            f3f = f3d.cast<float>();
            REQUIRE(f3f.cols() == n);
            REQUIRE(f3f.rows() == n * n);
            REQUIRE(f3d.cols() == n);
            REQUIRE(f3d.rows() == n * n);
            REQUIRE(c3f.cols() == 4);
            REQUIRE(c3f.rows() == 3);
            REQUIRE(c3d.cols() == 4);
            REQUIRE(c3d.rows() == 3);
            CHECK(getF3(f3d, 1, 2, 3) == get444(f3d, 1, 2, 3));
            BENCHMARK(title("Field3D_MXYZ<float>(", f3f.rows(), ", ", f3f.cols(), ") cubic1"), i)
            {
                auto x = getx(i, n);
                auto y = gety(i, n);
                auto z = getz(i, n);
                return cubic<float>::cubicN(c3f.row(0), getF1(f3f, x, y, z));
            };
            BENCHMARK(title("Field3D_MXYZ<double>(", f3f.rows(), ", ", f3f.cols(), ") cubic1"), i)
            {
                auto x = getx(i, n);
                auto y = gety(i, n);
                auto z = getz(i, n);
                return cubic<double>::cubicN(c3d.row(0), getF1(f3d, x, y, z));
            };
            BENCHMARK(title("Field3D_MXYZ<float> cubic2 size=", n, "x", n), i)
            {
                auto x = getx(i, n);
                auto y = gety(i, n);
                auto z = getz(i, n);
                return cubic<float>::cubicN(c3f.topRows(2), getF2(f3f, x, y, z));
            };

            BENCHMARK(title("Field3D_MXYZ<double> cubic2 size=", n, "x", n), i)
            {
                auto x = getx(i, n);
                auto y = gety(i, n);
                auto z = getz(i, n);
                return cubic<double>::cubicN(c3d.topRows(2), getF2(f3d, x, y, z));
            };
            BENCHMARK(title("cubicN(c3f, getF3(MXYZf, x, y, z) size=", n, "x", n), i)
            {
                auto x = getx(i, n);
                auto y = gety(i, n);
                auto z = getz(i, n);
                return cubic<float>::cubicN(c3f, getF3(f3f, x, y, z));
            };
            BENCHMARK(title("cubic3stdvec(c3f, getF3stdvec(MXYZf, x, y, z)) size=", n), i)
            {
                auto x = getx(i, n);
                auto y = gety(i, n);
                auto z = getz(i, n);
                return cubic3stdvec(c3f, getF3stdvec(f3f, x, y, z));
            };
            BENCHMARK(title("cubicN(c3d, MXYZd(blkidx_VXY(x,y,n,4), seqN(z, fix<4>)).reshaped(fix<4>,fix<4*4>)) size=", n), i)
            {
                auto x = getx(i, n);
                auto y = gety(i, n);
                auto z = getz(i, n);
                return cubic<double>::cubicN(c3d, f3d(blkidx_VXY(x, y, n, 4), seqN(z, fix<4>)).reshaped(fix<4>, fix<4 * 4>));
            };
            BENCHMARK(title("cubicN(c3d, MXYZd.reshaped(n,n*n)(seqN(x, fix<4>), blkidx_VXY(y,z,n,4))) size=", n), i)
            {
                auto x = getx(i, n);
                auto y = gety(i, n);
                auto z = getz(i, n);
                return cubic<double>::cubicN(c3d, f3d.reshaped(n, n * n)(seqN(x, fix<4>), blkidx_VXY(y, z, n, 4)));
            };
            BENCHMARK(title("cubicN(c3d, getF3(MXYZd, x, y, z)) size=", n), i)
            {
                auto x = getx(i, n);
                auto y = gety(i, n);
                auto z = getz(i, n);
                return cubic<double>::cubicN(c3d, getF3(f3d, x, y, z));
            };
            BENCHMARK(title("cubicN(c3d, get444(MXYZd,x,y,z)) size=", n), i)
            {
                auto x = getx(i, n);
                auto y = gety(i, n);
                auto z = getz(i, n);
                return cubic<double>::cubicN(c3d, get444(f3d, x, y, z));
            };
            BENCHMARK(title("cubic3stdvec(c3d, getF3stdvec(MXYZd, x, y, z)) size=", n), i)
            {
                auto x = getx(i, n);
                auto y = gety(i, n);
                auto z = getz(i, n);
                return cubic3stdvec(c3d, getF3stdvec(f3d, x, y, z));
            };
        };
    };
};
