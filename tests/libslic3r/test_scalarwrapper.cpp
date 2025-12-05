#define NOMINMAX

#include <catch2/catch.hpp>
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
SCENARIO("Test Eigen's Nesting an Array inside an Array", "[ScalarWrapper]")
{
    GIVEN("Nested Array types are defined and instances are initialized.")
    {
        typedef Array<Array3f, 3, 1> Array3A3f;
        Array3f                      a0(1, 2, 3), a1(2, 3, 4), a2(3, 4, 5), a3(4, 5, 6);
        Array3A3f                    a3a3f(a0, a1, a2);
        THEN("Check basic array operations")
        {
            // These check basic Array operations to show what works and what
            // doesn't, and to verify that any defined eigen extensions
            // haven't broken their compilation or execution.
            //
            // Note the commented out checks show various different Array
            // operations that don't work. This includes showing that `==`
            // for arrays generates a component-wise array of bool which
            // cannot be cast into a bool for the check. However, `==` for
            // Matrix does generate a single bool.
            //
            // CHECK(Array3f(0,0,0) == Array3f(0,0,0));  // cannot check if array == array.
            CHECK((Array3f(0, 0, 0) == Array3f(0, 0, 0)).all());                              // check (array == array).all().
            CHECK(Array3f(1, 2, 3).matrix() == Vector3f(1, 2, 3));                            // check array.matrix() == vector.
            CHECK((Array3f(1, 2, 3) + 1.0f).matrix() == Vector3f(2, 3, 4));                   // add array + scalar.
            CHECK((1 + Array3f(1, 2, 3)).matrix() == Vector3f(2, 3, 4));                      // add scalar + array.
            CHECK((Array3f(1, 2, 3) - 1.0f).matrix() == Vector3f(0, 1, 2));                   // subtract array - scalar.
            CHECK((1 - Array3f(1, 2, 3)).matrix() == Vector3f(0, -1, -2));                    // subtract scalar - array.
            CHECK((2.0f * Array3f(1, 2, 3)).matrix() == Vector3f(2, 4, 6));                   // multiply scalar * array.
            CHECK((Array3f(1, 2, 3) * 2).matrix() == Vector3f(2, 4, 6));                      // multiply array * scalar.
            CHECK((2.0f / Array3f(1, 2, 3)).matrix() == Vector3f(2.0 / 1, 2.0 / 2, 2.0 / 3)); // divide scalar / array.
            CHECK((Array3f(1, 2, 3) / 2).matrix() == Vector3f(1 / 2.0, 2 / 2.0, 3 / 2.0));    // divide array / scalar.
            CHECK((Array3f(1, 2, 3) + Array3f(1, 2, 3)).matrix() == Vector3f(2, 4, 6));       // add arrays.
            CHECK((Array3f(1, 2, 3) - Array3f(1, 2, 3)).matrix() == Vector3f(0, 0, 0));       // subtract arrays.
            CHECK((Array3f(1, 2, 3) * Array3f(1, 2, 3)).matrix() == Vector3f(1, 4, 9));       // multiply arrays.
            CHECK((Array3f(1, 2, 3) / Array3f(1, 2, 3)).matrix() == Vector3f(1, 1, 1));       // divide arrays.
            // CHECK((Array3f(1, 2, 3) + Array3d(1, 2, 3)).matrix() == Vector3d(1, 1, 1));       // add arrays with promotable scalars.
        };
        WHEN("Array<Array3f> instances are added.")
        {
            auto ans = a3a3f + a3a3f;
            THEN("Check a3a3f + a3a3f = 2 * a3a3f.")
            {
                CHECK(ans(0).matrix() == (2 * a0).matrix());
                CHECK(ans(1).matrix() == (2 * a1).matrix());
                CHECK(ans(2).matrix() == (2 * a2).matrix());
            };
        };
        WHEN("Array<Array3f> instances are subtracted.")
        {
            auto ans  = a3a3f - a3a3f;
            auto zero = Array3A3f::Zero();
            THEN("Check a3a3f - a3a3f = Array3A3f::Zero().")
            {
                CHECK(ans(0).matrix() == zero(0).matrix());
                CHECK(ans(1).matrix() == zero(1).matrix());
                CHECK(ans(2).matrix() == zero(2).matrix());
            };
        };
        WHEN("Array<Array3f> instances are multiplied.")
        {
            auto ans = a3a3f * a3a3f;
            THEN("Check a3a3f * a3a3f results.")
            {
                CHECK(ans(0).matrix() == (a0 * a0).matrix());
                CHECK(ans(1).matrix() == (a1 * a1).matrix());
                CHECK(ans(2).matrix() == (a2 * a2).matrix());
            };
        };
        WHEN("Array<Array3f> instances are divided.")
        {
            auto ans = a3a3f / a3a3f;
            THEN("Check a3a3f / a3a3f results.")
            {
                CHECK(ans(0).matrix() == (a0 / a0).matrix());
                CHECK(ans(1).matrix() == (a1 / a1).matrix());
                CHECK(ans(2).matrix() == (a2 / a2).matrix());
            };
        };
        WHEN("Array<Array3f> is added to Array3f.")
        {
            auto ans = a3a3f + a3;
            THEN("Check a3a3f + a3 results.")
            {
                CHECK(ans(0).matrix() == (a0 + a3).matrix());
                CHECK(ans(1).matrix() == (a1 + a3).matrix());
                CHECK(ans(2).matrix() == (a2 + a3).matrix());
            };
        };
        WHEN("Array3f is added to Array<Array3f>.")
        {
            auto ans = a3 + a3a3f;
            THEN("Check a3 + a3a3f results.")
            {
                CHECK(ans(0).matrix() == (a3 + a0).matrix());
                CHECK(ans(1).matrix() == (a3 + a1).matrix());
                CHECK(ans(2).matrix() == (a3 + a2).matrix());
            };
        };
        WHEN("Array<Array3f> is multiplied by Array3f.")
        {
            auto ans = a3a3f * a3;
            THEN("Check a3a3f * a3 results.")
            {
                CHECK(ans(0).matrix() == (a0 * a3).matrix());
                CHECK(ans(1).matrix() == (a1 * a3).matrix());
                CHECK(ans(2).matrix() == (a2 * a3).matrix());
            };
        };
        WHEN("Array3f is multiplied by Array<Array3f>.")
        {
            auto ans = a3 * a3a3f;
            THEN("Check a3 * a3a3f results.")
            {
                CHECK(ans(0).matrix() == (a3 * a0).matrix());
                CHECK(ans(1).matrix() == (a3 * a1).matrix());
                CHECK(ans(2).matrix() == (a3 * a2).matrix());
            };
        };
        /* The following doesn't work because a3a3f + 7.0f doesn't compile.
        WHEN("Array<Array3f> is added to float.")
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

SCENARIO("Test Eigen's Nesting an Array inside a Matrix", "[ScalarWrapper]")
{
    GIVEN("Nested Array types are defined and instances are initialized.")
    {
        typedef Matrix<Array3f, 2, 2> Matrix2A3f;
        Array3f                       a0(0, 0, 0), a1(0, 1, 1), a2(1, 0, 2), a3(1, 1, 3);
        Matrix2A3f                    m2a3f{{a0, a1}, {a2, a3}};
        Matrix2f                      m2f{{0, 1}, {2, 3}};
        THEN("Check basic matrix operations")
        {
            // These check basic Matrix operations to show what works and
            // what doesn't, and to verify that any defined eigen extensions
            // haven't broken their compilation or execution.
            //
            // Note the commented out checks show various different Matrix
            // operations that don't work. This includes showing that `Matrix
            // + Scalar` and `Scalar / Matrix` are not supported operations.
            //
            CHECK(Vector3f(1, 2, 3) == Vector3f(1, 2, 3));                // check vector == vector.
            CHECK((Vector3f(1, 2, 3).array() == Array3f(1, 2, 3)).all()); // check (vector.array() == array).all().
            // CHECK((Vector3f(1, 2, 3) + 1.0f) == Vector3f(2, 3, 4));       // add vector + scalar.
            // CHECK((1 + Vector3f(1, 2, 3)) == Vector3f(2, 3, 4));          // add scalar + vector.
            // CHECK((Vector3f(1, 2, 3) - 1.0f) == Vector3f(0, 1, 2));       // subtract vector - scalar.
            // CHECK((1 - Vector3f(1, 2, 3)) == Vector3f(0, -1, -2));        // subtract scalar - vector.
            CHECK((2.0f * Vector3f(1, 2, 3)) == Vector3f(2, 4, 6)); // multiply scalar * vector.
            CHECK((Vector3f(1, 2, 3) * 2) == Vector3f(2, 4, 6));    // multiply vector * scalar.
            // CHECK((2.0f / Vector3f(1, 2, 3)) == Vector3f(2.0 / 1, 2.0 / 2, 2.0 / 3));           // divide scalar / vector.
            CHECK((Vector3f(1, 2, 3) / 2) == Vector3f(1 / 2.0, 2 / 2.0, 3 / 2.0)); // divide vector / scalar.
            CHECK((Vector3f(1, 2, 3) + Vector3f(1, 2, 3)) == Vector3f(2, 4, 6));   // add vectors.
            CHECK((Vector3f(1, 2, 3) - Vector3f(1, 2, 3)) == Vector3f(0, 0, 0));   // subtract vectors.
            // CHECK((Vector3f(1, 2, 3) * Vector3f(1, 2, 3)) ==
            //       (Vector3f(1, 2, 3) * Vector3f(1, 2, 3))); // multiply vectors.
            // CHECK((Vector3f(1, 2, 3) / Vector3f(1, 2, 3)) ==
            //       (Vector3f(1, 2, 3) / Vector3f(1, 2, 3))); // divide vectors.
            CHECK((RowVector3f(1, 2, 3) * Vector3f(1, 2, 3)) == (RowVector3f(1, 2, 3) * Vector3f(1, 2, 3))); // multiply rowvector * vector.
            CHECK((Vector3f(1, 2, 3) * RowVector3f(1, 2, 3)) == (Vector3f(1, 2, 3) * RowVector3f(1, 2, 3))); // multiply vector * rowvector.
            // CHECK((Vector3f(1, 2, 3) + Vector3d(1, 2, 3)) == Vector3d(2, 4, 6)); // add vectors with promotable scalars.
        };
        WHEN("Matrix2<Array3f> instances are added.")
        {
            auto ans = m2a3f + m2a3f;
            THEN("Check m2a3f + m2a3f = 2 * m2a3f.")
            {
                CHECK(ans(0, 0).matrix() == (2 * a0).matrix());
                CHECK(ans(0, 1).matrix() == (2 * a1).matrix());
                CHECK(ans(1, 0).matrix() == (2 * a2).matrix());
                CHECK(ans(1, 1).matrix() == (2 * a3).matrix());
            };
        };
        WHEN("Matrix2<Array3f> instances are subtracted.")
        {
            auto ans  = m2a3f - m2a3f;
            auto zero = Matrix2A3f::Zero();
            THEN("Check m2a3f - m2a3f = Matrix2A3f::Zero().")
            {
                CHECK(ans(0, 0).matrix() == zero(0, 0).matrix());
                CHECK(ans(0, 1).matrix() == zero(0, 1).matrix());
                CHECK(ans(1, 0).matrix() == zero(1, 0).matrix());
                CHECK(ans(1, 1).matrix() == zero(1, 1).matrix());
            };
        };
        WHEN("Matrix2<Array3f> instances are multiplied.")
        {
            auto ans = m2a3f * m2a3f;
            THEN("Check m2a3f * m2a3f results.")
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
        WHEN("Matrix2<Array3f> is multiplied by Array3f.")
        {
            auto ans = m2a3f * a3;
            THEN("Check m2a3f * a3 results.")
            {
                CHECK(ans(0, 0).matrix() == (a0 * a3).matrix());
                CHECK(ans(0, 1).matrix() == (a1 * a3).matrix());
                CHECK(ans(1, 0).matrix() == (a2 * a3).matrix());
                CHECK(ans(1, 1).matrix() == (a3 * a3).matrix());
            };
        };
        WHEN("Array3f is multiplied by Matrix2<Array3f>.")
        {
            auto ans = a3 * m2a3f;
            THEN("Check a3 * m2a3f results.")
            {
                CHECK(ans(0, 0).matrix() == (a3 * a0).matrix());
                CHECK(ans(0, 1).matrix() == (a3 * a1).matrix());
                CHECK(ans(1, 0).matrix() == (a3 * a2).matrix());
                CHECK(ans(1, 1).matrix() == (a3 * a3).matrix());
            };
        };
        WHEN("Matrix2<Array3f> is divided by Array3f.")
        {
            auto ans = m2a3f / a3;
            THEN("Check m2a3f / a3 results.")
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
                CHECK((m2v3f.array() + scalar(v1)).matrix() == Matrix2V3f({{scalar(v0 + v1), scalar(v1 + v1)}, {scalar(v2 + v1), scalar(v3 + v1)}}));
                CHECK((scalar(v2) + m2v3f.array()).matrix() == Matrix2V3f({{scalar(v2 + v0), scalar(v2 + v1)}, {scalar(v2 + v2), scalar(v2 + v3)}}));
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
