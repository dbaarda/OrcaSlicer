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

SCENARIO("Test Nesting a Vector3f inside a Matrix", "[EigenNesting]")
{
    GIVEN("Classes are defined and instances are initialized.")
    {
        typedef ScalarMatrix<Vector3f>  ScalarV3f;
        typedef Matrix<ScalarV3f, 2, 2> Matrix2V3f;
        Vector3f                        v0(0, 0, 0), v1(0, 1, 1), v2(1, 0, 2), v3(1, 1, 3);
        Matrix2V3f                      m2v3f;
        m2v3f(0, 0)  = v0;
        m2v3f(0, 1)  = v1;
        m2v3f(1, 0)  = v2;
        m2v3f(1, 1)  = v3;
        Matrix2f m2f = (Matrix2f() << 0, 1, 2, 3).finished();
        // Matrix2V3f m1(v0,v1,v2,v3);
        THEN("Check initialized values")
        {
            CHECK(m2v3f(0, 0) == v0);
            CHECK(m2v3f(0, 1) == v1);
            CHECK(m2v3f(1, 0) == v2);
            CHECK(m2v3f(1, 1) == v3);
        };
        WHEN("Instances added.")
        {
            THEN("Check m2<v3<f>> + m2<v3<f>> = N * m2<v3<f>> operations.")
            {
                CHECK(m2v3f + m2v3f == 2 * m2v3f);
                CHECK(m2v3f + m2v3f + m2v3f == 3 * m2v3f);
                CHECK(2 * m2v3f + m2v3f == m2v3f * 3);
                CHECK(m2v3f - m2v3f == Matrix2V3f::Zero());
            };
            THEN("Check m2<v3<f>>.array() + scalar(v3<f>).")
            {
                // CHECK((m2v3f.array() + as_scalar(v0)) == (as_scalar(v0) + m2v3f.array()));
                CHECK((m2v3f.array() + as_scalar(v0)).matrix() == m2v3f);
                CHECK((m2v3f.array() + as_scalar(v1))(0, 0) == v0 + v1);
                CHECK((m2v3f.array() + as_scalar(v1))(0, 1) == v1 + v1);
                CHECK((m2v3f.array() + as_scalar(v1))(1, 0) == v2 + v1);
                CHECK((m2v3f.array() + as_scalar(v1))(1, 1) == v3 + v1);
            };
        };
        WHEN("Instances multiplied.")
        {
            THEN("Check m2<f> * scalar(v3<f>) gives m2<v3f>.")
            {
                CHECK((m2f * as_scalar(v1)) == (as_scalar(v1) * m2f));
                CHECK((m2f * as_scalar(v0)) == Matrix2V3f::Zero());
                CHECK((m2f * as_scalar(v1))(0, 0) == v0);
                CHECK((m2f * as_scalar(v1))(0, 1) == v1);
                CHECK((m2f * as_scalar(v1))(1, 0) == 2 * v1);
                CHECK((m2f * as_scalar(v1))(1, 1) == 3 * v1);
            };
            THEN("Check m2<v3<f>> * N gives m2<v3<f>>.")
            {
                CHECK((m2v3f * 2) == (2 * m2v3f));
                CHECK((m2v3f * 0) == Matrix2V3f::Zero());
                CHECK((m2v3f * 2)(0, 0) == 2 * v0);
                CHECK((m2v3f * 2)(0, 1) == 2 * v1);
                CHECK((m2v3f * 2)(1, 0) == 2 * v2);
                CHECK((m2v3f * 2)(1, 1) == 2 * v3);
            };
            THEN("Check m2<f> * m2<v3<f>> gives m2<v3<f>>."){
                // Matrix2V3f ans = m2v3f * m2f;
                // CHECK((m2f*m2v3f) != (m2v3f*m2f));
                // CHECK((m2f*m2v3f)(0,0) == m2f(0,0)*m2v3f(0,0) + m2f(0,1)*m2v3f(1,0));
                // CHECK((m2f*m2v3f)(0,1) == m2f(0,0)*m2v3f(0,1) + m2f(0,1)*m2v3f(1,1));
                // CHECK((m2f*m2v3f)(1,0) == m2f(1,0)*m2v3f(0,0) + m2f(1,1)*m2v3f(1,0));
                // CHECK((m2f*m2v3f)(1,1) == m2f(1,0)*m2v3f(0,1) + m2f(1,1)*m2v3f(1,1));
            };
        };
    };
}
