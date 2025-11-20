#define NOMINMAX

#include <catch2/catch.hpp>
#include <test_utils.hpp>

#include <fstream>

#include <libslic3r/Geometry/EigenNesting.hpp>

using namespace Eigen;
using namespace Catch::Matchers;

/*
 * For docs on using catch2;
 *
 * https://catch2-temp.readthedocs.io/en/latest/test-cases-and-sections.html
 */

// Methods to wrap/unwrap using casting.
// template<typename Matrix_> ScalarWrapper<Matrix_>& _wrap(Matrix_& m) { return static_cast<ScalarWrapper<Matrix_>&>(m); }
// template<typename Matrix_> Matrix_&               _unwrap(ScalarWrapper<Matrix_>& m) { return static_cast<Matrix_&>(m); }

// define the wrapping class and conversion methods we want to test.
#define Wrapper ScalarWrapper
#define wrap(a) scalar(a)
// #define wrap(a) ScalarWrapper(a)
#define unwrap(a) (a).nestedExpression()

/*
template<typename MatrixExpr>
    const ArrayBase<ArrayWrapper<MatrixExpr>>& as_ArrayWrapper(const MatrixBase<MatrixExpr>& m)
{
  return ArrayWrapper(m);
}
*/

SCENARIO("Test Nesting a Vector3f inside a Matrix", "[EigenNesting]")
{
    GIVEN("Classes are defined and instances are initialized.")
    {
        typedef Wrapper<Vector3f>       ScalarV3f;
        typedef Matrix<ScalarV3f, 2, 2> Matrix2V3f;
        Vector3f                        v0(0, 0, 0), v1(0, 1, 1), v2(1, 0, 2), v3(1, 1, 3);
        Matrix2V3f                      m2v3f;
        m2v3f(0, 0) = v0;
        m2v3f(0, 1) = v1;
        m2v3f(1, 0) = v2;
        m2v3f(1, 1) = v3;
        Matrix2f m2f{{0, 1}, {2, 3}};
        // Matrix2V3f m1{{v0,v1},{v2,v3}};
        THEN("Check initialized values")
        {
            // Note the commented out checks show various different Array and
            // Matrix operations that don't work. This includes showing that
            // `==` for arrays generates a component-wise array of bool which
            // cannot be cast into a bool for the check. However, `==` for
            // Matrix does generate a single bool.
            CHECK(Vector3f(0, 0, 0) == Vector3f(0, 0, 0)); // can check if matrix == matrix
            // CHECK(Array3f(0,0,0) == Array3f(0,0,0));  // cannot check if array == array
            // CHECK((Vector3f(1,1,1) + 1.0f) == Vector3f(2,2,2));  // cannot add matrix + scalar
            CHECK((Array3f(1, 1, 1) + 1.0f).matrix() == Vector3f(2, 2, 2)); // can add array + scalar
            CHECK((2.0f * Vector3f(1, 1, 1)) == Vector3f(2, 2, 2));         // can multiply scalar and matrix.
            CHECK((2.0f * Array3f(1, 1, 1)).matrix() == Vector3f(2, 2, 2)); // can multiply scalar and array

            CHECK((Vector3f(1, 1, 1).array() + 1.0f).matrix() == Vector3f(2, 2, 2)); // can add matrix.array() + scalar
            CHECK((2.0f * Vector3f(1, 1, 1).array()).matrix() == Vector3f(2, 2, 2)); // can multiply scalar and matrix.array()
            CHECK(Matrix2f{{0.0f, 1.0f}, {2.0f, 3.0f}} == m2f);

            CHECK(v0 == Vector3f(0, 0, 0));
            CHECK(unwrap(wrap(v0)) == Vector3f(0, 0, 0));
            CHECK(unwrap(m2v3f(0, 0)) == v0);
            CHECK(unwrap(m2v3f(0, 1)) == v1);
            CHECK(unwrap(m2v3f(1, 0)) == v2);
            CHECK(unwrap(m2v3f(1, 1)) == v3);
        };
        WHEN("Instances added.")
        {
            THEN("Check m2<v3<f>> + m2<v3<f>> = N * m2<v3<f>> operations.")
            {
                CHECK(m2v3f + m2v3f == 2.0f * m2v3f);
                CHECK(m2v3f + m2v3f + m2v3f == 3.0f * m2v3f);
                CHECK(2.0f * m2v3f + m2v3f == m2v3f * 3.0f);
                CHECK(m2v3f - m2v3f == Matrix2V3f::Zero());
            };
            THEN("Check m2<v3<f>>.array() + scalar(v3<f>).")
            {
                // CHECK((m2v3f.array() + wrap(v0)) == (wrap(v0) + m2v3f.array()));
                CHECK((m2v3f.array() + wrap(v0)).matrix() == m2v3f);
                CHECK(unwrap((m2v3f.array() + wrap(v1))(0, 0)) == v0 + v1);
                CHECK(unwrap((m2v3f.array() + wrap(v1))(0, 1)) == v1 + v1);
                CHECK(unwrap((m2v3f.array() + wrap(v1))(1, 0)) == v2 + v1);
                CHECK(unwrap((m2v3f.array() + wrap(v1))(1, 1)) == v3 + v1);
            };
        };
        WHEN("Instances multiplied.")
        {
            THEN("Check m2<f> * scalar(v3<f>) gives m2<v3f>.")
            {
                CHECK(Matrix2V3f::Zero() == 0.0f * m2v3f);
                CHECK((m2f * wrap(v1)) == (wrap(v1) * m2f));
                CHECK((m2f * wrap(v0)) == Matrix2V3f::Zero());
                // CHECK(m2f * wrap(v1) == Matrix2V3f({{v1*0.0f, v1*1.0f}, {v1*2.0f, v1*3.0f}}));
                CHECK(m2f * wrap(v1) == Matrix2V3f({{wrap(v1 * 0.0f), wrap(v1 * 1.0f)}, {wrap(v1 * 2.0f), wrap(v1 * 3.0f)}}));
                CHECK(unwrap((m2f * wrap(v1))(0, 0)) == v0);
                CHECK(unwrap((m2f * wrap(v1))(0, 1)) == v1);
                CHECK(unwrap((m2f * wrap(v1))(1, 0)) == 2 * v1);
                CHECK(unwrap((m2f * wrap(v1))(1, 1)) == 3 * v1);
            };
            THEN("Check m2<v3<f>> * N gives m2<v3<f>>.")
            {
                CHECK((m2v3f * 2.0f) == (2.0f * m2v3f));
                CHECK((m2v3f * 0.0f) == Matrix2V3f::Zero());
                CHECK(unwrap((m2v3f * 2.0f)(0, 0)) == 2 * v0);
                CHECK(unwrap((m2v3f * 2.0f)(0, 1)) == 2 * v1);
                CHECK(unwrap((m2v3f * 2.0f)(1, 0)) == 2 * v2);
                CHECK(unwrap((m2v3f * 2.0f)(1, 1)) == 2 * v3);
            };
            THEN("Check m2<f> * m2<v3<f>> gives m2<v3<f>>.")
            {
                Matrix2V3f ans = m2v3f * m2f;
                CHECK((m2f * m2v3f) != (m2v3f * m2f));
                CHECK(unwrap((m2f * m2v3f)(0, 0)) == unwrap(m2f(0, 0) * m2v3f(0, 0) + m2f(0, 1) * m2v3f(1, 0)));
                CHECK(unwrap((m2f * m2v3f)(0, 1)) == unwrap(m2f(0, 0) * m2v3f(0, 1) + m2f(0, 1) * m2v3f(1, 1)));
                CHECK(unwrap((m2f * m2v3f)(1, 0)) == unwrap(m2f(1, 0) * m2v3f(0, 0) + m2f(1, 1) * m2v3f(1, 0)));
                CHECK(unwrap((m2f * m2v3f)(1, 1)) == unwrap(m2f(1, 0) * m2v3f(0, 1) + m2f(1, 1) * m2v3f(1, 1)));
            };
        };
    };
}
