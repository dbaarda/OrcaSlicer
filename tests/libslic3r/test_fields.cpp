#define NOMINMAX

#include <catch2/catch_all.hpp>
#include <test_utils.hpp>

#include <fstream>

#include <libslic3r/Geometry/Fields.hpp>
#include <libslic3r/Point.hpp>
#include <libslic3r/SVG.hpp>
#include <libslic3r/SLA/AGGRaster.hpp>
#include <libslic3r/StreamUtils.hpp>

using namespace Slic3r;
using namespace Eigen;
using namespace Slic3r::Geometry;
using namespace Catch::Matchers;

/*
 * For docs on using catch2;
 *
 * https://catch2-temp.readthedocs.io/en/latest/test-cases-and-sections.html
 */

static Slic3r::sla::RasterGrayscaleAAGammaPower create_raster(const size_t w = 100, const size_t h = 100, const double pxw = 0.1)
{ return sla::RasterGrayscaleAAGammaPower(sla::Resolution(w, h), sla::PixelDim(pxw, pxw), sla::RasterBase::Trafo()); }

/*
    // drawing to a png file.
    auto rst = create_raster();
    rst.draw(expoly);
    std::fstream out(name + ".png", std::ios::out);
    out << rst.encode(sla::PNGRasterEncoder{});
    out.close();

    // drawing to an svg file.
    auto raster_bb =  BoundingBox(Point::new_scale(-50., -50.), Point::new_scale(50.,50.));
    SVG svg(name + ".svg", raster_bb);
    svg.draw_grid(raster_bb, "grey", scale_(0.05), scale_(5.0));
    svg.draw_outline(expoly, "red", "red", scale_(0.3));
    svg.Close();

// sample test.
TEST_CASE("Test Field2D", "[Fields]")
{
    SECTION("")
    {
      INFO("Comparing ext[" << i << "] against closest ref[" << j << "]");
      UNSCOPED_INFO("extracted ext[" << i << "]: " << ext[i]);
      REQUIRE(val == tgt);
      CHECK(val == tgt);
      CHECK_THAT(val, WithinRel(tgt, 0.05) || WithinAbs(tgt, 1.0));
      BENCHMARK("indexed", i) { return bench(i); };
    }
}

 */

typedef SVec3<float>   SVec3f;
typedef SVec3<double>  SVec3d;
typedef Field2D<float>  Field2Df;
typedef Field2D<double> Field2Dd;
typedef Field2D<SVec3f> Field2DV3f;
typedef Field2D<SVec3d> Field2DV3d;
typedef Field3D<float>  Field3Df;
typedef Field3D<double> Field3Dd;
typedef Field3D<SVec3f> Field3DV3f;
typedef Field3D<SVec3d> Field3DV3d;

SCENARIO("Test Field2D Smoothing", "[Fields]")
{
    GIVEN("Setup f2(15,15) instance zeroed with f2(7,7)=1000.")
    {
        Field2D<float> f2(15, 15);

        f2.setZero();
        f2(7, 7) = 1000.0;

        THEN("Check f2 setup center.")
        {
            REQUIRE(f2(7, 7) == 1000.0);
            REQUIRE(f2(6, 6) == 0.0);
            REQUIRE(f2(8, 8) == 0.0);
        };
        WHEN("Apply smoothing for 1 iteration")
        {
            _smoothF2(f2, 1);
            THEN("Check smoothing result")
            {
                // Note we put brackets around the comparison because we
                // don't need to see all the zero's in the test output.
                CHECK((f2.topRows(6) == Matrix<float, 6, 15>::Zero()));
                CHECK((f2.bottomRows(6) == Matrix<float, 6, 15>::Zero()));
                CHECK((f2.leftCols(6) == Matrix<float, 15, 6>::Zero()));
                CHECK((f2.rightCols(6) == Matrix<float, 15, 6>::Zero()));
                CHECK(f2.block(6, 6, 3, 3).isApprox(Matrix<float, 3, 3>::Constant(1000.0 / 9.0)));
            };
        };
        WHEN("Apply smoothing for default=6 iterations")
        {
            _smoothF2(f2);
            THEN("Check smoothing result")
            {
                CHECK((f2.topRows(1) == Matrix<float, 1, 15>::Zero()));
                CHECK((f2.bottomRows(1) == Matrix<float, 1, 15>::Zero()));
                CHECK((f2.leftCols(1) == Matrix<float, 15, 1>::Zero()));
                CHECK((f2.rightCols(1) == Matrix<float, 15, 1>::Zero()));
                // This is mostly so we can see the full block in the test output.
                CHECK(f2.block(1, 1, 13, 13) != Matrix<float, 13, 13>::Zero());
                CHECK(f2.maxCoeff() == f2(7, 7));
                CHECK(f2.block(1, 1, 13, 13).minCoeff() == f2(1, 1));
                CHECK(f2.block(1, 1, 13, 13).minCoeff() > 0.0);
                // Result should be symetric about (7,7) and smaller as you go out.
                for (Index i = 0; i < 7; i++) {
                    INFO("Checking for i=" << i);
                    CHECK((f2.row(i) == f2.row(14 - i)));
                    CHECK((f2.col(i) == f2.col(14 - i)));
                    CHECK((f2.col(i) == f2.row(i).transpose()));
                    CHECK((f2.row(i).segment(1, 13).array() < f2.row(i + 1).segment(1, 13).array()).all());
                    CHECK((f2.col(i).segment(1, 13).array() < f2.col(i + 1).segment(1, 13).array()).all());
                }
            };
        };
    };
}

SCENARIO("Test Field3D Smoothing", "[Fields]")
{
    GIVEN("Setup f3(15,15,15) instance zeroed with f3(7,7,7)=1000.")
    {
        Field3D<float> f3(15, 15, 15);

        f3.setZero();
        f3(7, 7, 7) = 1000.0;

        THEN("Check f3 setup center.")
        {
            REQUIRE(f3(7, 7, 7) == 1000.0);
            REQUIRE(f3(6, 6, 6) == 0.0);
            REQUIRE(f3(8, 8, 8) == 0.0);
        };
        WHEN("Apply smoothing for 1 iteration")
        {
            _smoothF3(f3, 1);
            THEN("Check smoothing result")
            {
                CHECK(f3(7, 7, 5) == 0.0);
                CHECK(f3(7, 7, 9) == 0.0);
                CHECK(f3(7, 5, 7) == 0.0);
                CHECK(f3(7, 9, 7) == 0.0);
                CHECK((f3(5) == Matrix<float, 15, 15>::Zero()));
                CHECK((f3(9) == Matrix<float, 15, 15>::Zero()));
                for (Index z = 6; z < 9; z++) {
                    CHECK(f3(z).block(6, 6, 3, 3).isApprox(Matrix<float, 3, 3>::Constant(1000.0 / 27.0)));
                }
            };
        };
        WHEN("Apply smoothing for default=6 iterations")
        {
            _smoothF3(f3);
            THEN("Check smoothing result")
            {
                CHECK((f3(0) == Matrix<float, 15, 15>::Zero()));
                CHECK((f3(14) == Matrix<float, 15, 15>::Zero()));
                for (Index z = 1; z < 14; z++) {
                    INFO("Checking for z=" << z);
                    CHECK((f3(z).topRows(1) == Matrix<float, 1, 15>::Zero()));
                    CHECK((f3(z).bottomRows(1) == Matrix<float, 1, 15>::Zero()));
                    CHECK((f3(z).leftCols(1) == Matrix<float, 15, 1>::Zero()));
                    CHECK((f3(z).rightCols(1) == Matrix<float, 15, 1>::Zero()));
                    // This is mostly so we can see the full block in the test output.
                    CHECK(f3(z).block(1, 1, 13, 13) != Matrix<float, 13, 13>::Zero());
                    CHECK(f3(z).maxCoeff() == f3(z)(7, 7));
                    CHECK(f3(z).block(1, 1, 13, 13).minCoeff() == f3(z)(1, 1));
                    CHECK(f3(z).block(1, 1, 13, 13).minCoeff() > 0.0);
                }
                for (Index z = 0; z < 7; z++) {
                    INFO("Checking for z=" << z);
                    CHECK(f3(z).isApprox(f3(14 - z)));
                    CHECK((f3(z).block(1, 1, 13, 13).array() < f3(z + 1).block(1, 1, 13, 13).array()).all());
                }
            };
        };
    };
}

TEST_CASE("Benchmark Field2D.cubic()", "[Fields]")
{
    Eigen::Index n     = 100;
    Eigen::Index steps = n - 3;
    double       dstep = double(n) / double(steps);
    SECTION("Field2D.cubic()")
    {
        std::stringstream title;
        Field2D<double>    f2d(n, n);
        Field2D<float>     f2f(n, n);
        f2d.setRandom();
        f2f = f2d.cast<float>();
        title << "Field2D<float> size=" << n << "x" << n;
        BENCHMARK(title.str(), i)
        {
            float d = (i % steps) * dstep;
            return cubic(f2f, Vec2f(d, d));
        };
        title = std::stringstream();
        title << "Field2D<double> size=" << n << "x" << n;
        BENCHMARK(title.str(), i)
        {
            double d = (i % steps) * dstep;
            return cubic(f2d, Vec2d(d, d));
        };
    };
};

TEST_CASE("Benchmark Field3D.cubic()", "[Fields]")
{
    // We use steps and dstep as a cheap way to calculate semi-random floats
    // to interpolate in the range 0 to n-1.
    Eigen::Index n     = 100;
    Eigen::Index steps = n - 3;
    double       dstep = double(n) / double(steps);
    SECTION("Field3D.cubic()")
    {
        std::stringstream title;
        Field3D<double>    f3d(n, n, n);
        Field3D<float>     f3f(n, n, n);
        f3d.setRandom();
        f3f = f3d.cast<float>();
        title << "Field3D<float> size=" << n << "x" << n << "x" << n;
        BENCHMARK(title.str(), i)
        {
            float d = (i % steps) * dstep;
            return cubic(f3f, Vec3f(d, d, d));
        };
        title = std::stringstream();
        title << "Field3D<double> size=" << n << "x" << n << "x" << n;
        BENCHMARK(title.str(), i)
        {
            double d = (i % steps) * dstep;
            return cubic(f3d, Vec3d(d, d, d));
        };
    };
};

TEST_CASE("Benchmark Field3D cubic_z()", "[Fields]")
{
    Eigen::Index n     = 100;
    Eigen::Index steps = n - 3;
    double       dstep = double(n) / double(steps);
    SECTION("Field3D.cubic_z()")
    {
        std::stringstream     title;
        Field3D<double>        f3xd(n, n, n), f3yd(n, n, n), f3zd(n, n, n);
        Field3D<float>         f3xf(n, n, n), f3yf(n, n, n), f3zf(n, n, n);
        Field3D<SVec3<double>> f3v3d(n, n, n);
        Field3D<SVec3<float>>  f3v3f(n, n, n);
        f3xd.setRandom();
        f3yd.setRandom();
        f3zd.setRandom();
        f3xf = f3xd.cast<float>();
        f3yf = f3yd.cast<float>();
        f3zf = f3zd.cast<float>();
        for (int z = 0; z < n; z++) {
            for (Eigen::Index y = 0; y < n; y++) {
                for (Eigen::Index x = 0; x < n; x++) {
                    f3v3d(x, y, z) = Vec3d(f3xd(x, y, z), f3yd(x, y, z), f3zd(x, y, z));
                    f3v3f(x, y, z) = Vec3f(f3xf(x, y, z), f3yf(x, y, z), f3zf(x, y, z));
                }
            }
        }
        title << "3 x Field3D<float> size=" << n << "x" << n << "x" << n;
        BENCHMARK(title.str(), i)
        {
            float         d    = (i % steps) * dstep;
            Field2D<float> f2xf = cubic_z(f3xf, d);
            Field2D<float> f2yf = cubic_z(f3yf, d);
            Field2D<float> f2zf = cubic_z(f3zf, d);
            return Vec3f(f2xf(1, 1), f2yf(1, 1), f2zf(1, 1));
        };
        title = std::stringstream();
        title << "3 x Field3D<double> size=" << n << "x" << n << "x" << n;
        BENCHMARK(title.str(), i)
        {
            double         d    = (i % steps) * dstep;
            Field2D<double> f2xd = cubic_z(f3xd, d);
            Field2D<double> f2yd = cubic_z(f3yd, d);
            Field2D<double> f2zd = cubic_z(f3zd, d);
            return Vec3d(f2xd(1, 1), f2yd(1, 1), f2zd(1, 1));
        };
        title = std::stringstream();
        title << "3 x Field3D<Vec3f> size=" << n << "x" << n << "x" << n;
        BENCHMARK(title.str(), i)
        {
            float f = (i % steps) * dstep;
            // Field2DV3f f2v3f = cubic_z(f3v3f, f);
            // return f2v3f;
        };
    };
};

TEST_CASE("Benchmark Field2D smooth()", "[Fields]")
{
    SECTION("Field2D<float>.smooth()")
    {
        for (Index n = 50; n <= 200; n *= 2) {
            std::stringstream title;
            title << "f2<float> size=" << n << "x" << n;
            Field2D<float> f2(n, n);
            f2.setRandom();
            BENCHMARK(title.str())
            {
                smooth(f2);
                return f2;
            };
            BENCHMARK(title.str() + "vect")
            {
                _smoothF2(f2);
                return f2;
            };
        }
    };
    SECTION("Field2D<double>.smooth()")
    {
        for (Index n = 50; n <= 200; n *= 2) {
            std::stringstream title;
            title << "f2<double> size=" << n << "x" << n;
            Field2D<double> f2(n, n);
            f2.setRandom();
            BENCHMARK(title.str())
            {
                smooth(f2);
                return f2;
            };
            BENCHMARK(title.str() + "vect")
            {
                _smoothF2(f2);
                return f2;
            };
        }
    };
    /*
    SECTION("Using T=int32_t")
    {
        for (Index n = 50; n <= 200; n *= 2) {
            std::stringstream title;
            title << "size=" << n << "x" << n;
            Field2D<int32_t> f2(n, n);
            f2.setRandom();
            BENCHMARK(title.str())
            {
                smooth(f2);
                return f2;
            };
        }
    };
    */
};

TEST_CASE("Benchmark Field3D smooth()", "[Fields]")
{
    SECTION("Field3D<float>.smooth()")
    {
        for (Index n = 50; n <= 200; n *= 2) {
            std::stringstream title;
            title << "f3<float> size=" << n << "x" << n << "x" << n;
            Field3D<float> f3(n, n, n);
            f3.setRandom();
            BENCHMARK(title.str())
            {
                smooth(f3);
                return f3;
            };
            BENCHMARK(title.str() + "vect")
            {
                _smoothF3(f3);
                return f3;
            };
        }
    };
    SECTION("Field3D<double>.smooth()")
    {
        for (Index n = 50; n <= 200; n *= 2) {
            std::stringstream title;
            title << "f3<double> size=" << n << "x" << n << "x" << n;
            Field3D<double> f3(n, n, n);
            f3.setRandom();
            BENCHMARK(title.str())
            {
                smooth(f3);
                return f3;
            };
            BENCHMARK(title.str() + "vect")
            {
                _smoothF3(f3);
                return f3;
            };
        }
    };
};
