#define NOMINMAX

#include <catch2/catch.hpp>
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
TEST_CASE("Test Field2", "[Fields]")
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

SCENARIO("Test Field2 Smoothing", "[Fields]")
{
    GIVEN("Setup f2(15,15) instance zeroed with f2(7,7)=1000.")
    {
        Field2<float> f2(15, 15);

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
            smooth(f2, 1);
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
            smooth(f2);
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

SCENARIO("Test Field3 Smoothing", "[Fields]")
{
    GIVEN("Setup f3(15,15,15) instance zeroed with f3(7,7,7)=1000.")
    {
        Field3<float> f3(15, 15, 15);

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
            smooth(f3, 1);
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
            smooth(f3);
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

TEST_CASE("Benchmark Field2 cubic()", "[Fields]")
{
    Eigen::Index n     = 100;
    Eigen::Index steps = n - 3;
    double       dstep = double(n) / double(steps);
    SECTION("Using T=float")
    {
        std::stringstream title;
        title << "f2<float> size=" << n << "x" << n;
        Field2<float> f2(n, n);
        f2.setRandom(n, n);
        BENCHMARK(title.str(), i)
        {
            float d = (i % steps) * dstep;
            return cubic(f2, Vec2f(d, d));
        };
    };
    SECTION("Using T=double")
    {
        std::stringstream title;
        title << "f2<double> size=" << n << "x" << n;
        Field2<double> f2(n, n);
        f2.setRandom(n, n);
        BENCHMARK(title.str(), i)
        {
            double d = (i % steps) * dstep;
            return cubic(f2, Vec2d(d, d));
        };
    };
};

TEST_CASE("Benchmark Field3 cubic()", "[Fields]")
{
    Eigen::Index n     = 100;
    Eigen::Index steps = n - 3;
    double       dstep = double(n) / double(steps);
    SECTION("Using T=float")
    {
        std::stringstream title;
        title << "f3<float> size=" << n << "x" << n << "x" << n;
        Field3<float> f3(n, n, n);
        for (Index z = 0; z < f3.size(); z++) {
            f3(z).setRandom(n, n);
        }
        BENCHMARK(title.str(), i)
        {
            float d = (i % steps) * dstep;
            return cubic(f3, Vec3f(d, d, d));
        };
    };
    SECTION("Using T=double")
    {
        std::stringstream title;
        title << "f3<double> size=" << n << "x" << n << "x" << n;
        Field3<double> f3(n, n, n);
        for (Index z = 0; z < f3.size(); z++) {
            f3(z).setRandom(n, n);
        }
        BENCHMARK(title.str(), i)
        {
            float d = (i % steps) * dstep;
            return cubic(f3, Vec3d(d, d, d));
        };
    };
    /*
    SECTION("Using T=Vec3f")
    {
        std::stringstream title;
        title << "f3<Vec3f> size=" << n << "x" << n << "x" << n;
        Field3<Vec3f> f3(n, n, n);
        for (Index z = 0; z < f3.lays(); z++) {
          for (Index y = 0; y < f3.rows(); y++) {
            for (Index x = 0; x < f3.cols(); x++) {
              f3(z,y,x) = Vec3f(x,y,z);
            }
          }
        }
        BENCHMARK(title.str(), i)
        {
          float d = (i % steps) * dstep;
          return cubic(f3, Vec3f(d, d, d));
        };
    };
     */
};

TEST_CASE("Benchmark Field3 cubic_z()", "[Fields]")
{
    Eigen::Index n     = 100;
    Eigen::Index steps = n - 3;
    double       dstep = double(n) / double(steps);
    SECTION("Using T=float")
    {
        std::stringstream title;
        title << "f3<float> size=" << n << "x" << n << "x" << n;
        Field3<float> f3(n, n, n);
        for (Index z = 0; z < f3.size(); z++) {
            f3(z).setRandom(n, n);
        }
        BENCHMARK(title.str(), i)
        {
            float d = (i % steps) * dstep;
            // TODO: make this work.
            //return cubic_z(f3, d);
        };
    };
    SECTION("Using T=double")
    {
        std::stringstream title;
        title << "f3<double> size=" << n << "x" << n << "x" << n;
        Field3<double> f3(n, n, n);
        for (Index z = 0; z < f3.size(); z++) {
            f3(z).setRandom(n, n);
        }
        BENCHMARK(title.str(), i)
        {
            double d = (i % steps) * dstep;
            // TODO: make this work.
            //return cubic_z(f3, d);
        };
    };
};

TEST_CASE("Benchmark Field2 smooth()", "[Fields]")
{
    SECTION("Using T=float")
    {
        for (Index n = 50; n <= 200; n *= 2) {
            std::stringstream title;
            title << "size=" << n << "x" << n;
            Field2<float> f2(n, n);
            f2.setRandom(n, n);
            BENCHMARK(title.str())
            {
                smooth(f2);
                return f2;
            };
        }
    };
    SECTION("Using T=double")
    {
        for (Index n = 50; n <= 200; n *= 2) {
            std::stringstream title;
            title << "size=" << n << "x" << n;
            Field2<double> f2(n, n);
            f2.setRandom(n, n);
            BENCHMARK(title.str())
            {
                smooth(f2);
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
            Field2<int32_t> f2(n, n);
            f2.setRandom(n, n);
            BENCHMARK(title.str())
            {
                smooth(f2);
                return f2;
            };
        }
    };
    */
};

TEST_CASE("Benchmark Field3 smooth()", "[Fields]")
{
    SECTION("Using T=float")
    {
        for (Index n = 50; n <= 200; n *= 2) {
            std::stringstream title;
            title << "size=" << n << "x" << n << "x" << n;
            Field3<float> f3(n, n, n);
            for (Index z = 0; z < f3.size(); z++) {
                f3(z).setRandom(n, n);
            }
            BENCHMARK(title.str())
            {
                smooth(f3);
                return f3;
            };
        }
    };
    SECTION("Using T=double")
    {
        for (Index n = 50; n <= 200; n *= 2) {
            std::stringstream title;
            title << "size=" << n << "x" << n << "x" << n;
            Field3<double> f3(n, n, n);
            for (Index z = 0; z < f3.size(); z++) {
                f3(z).setRandom(n, n);
            }
            BENCHMARK(title.str())
            {
                smooth(f3);
                return f3;
            };
        }
    };
};
