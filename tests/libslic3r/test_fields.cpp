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
    GIVEN("Setup f2(20,20) instance zeroed with f2(10,10)=1000.")
    {
        Field2<float> f2(20, 20);

        f2.setZero();
        f2(10, 10) = 1000.0;

        THEN("Check f2 instance.") { REQUIRE(f2(10, 10) == 1000.0); };
        WHEN("Apply smoothing for 1 iteration")
        {
            smooth(f2, 1);
            THEN("Check smoothing result")
            {
                CHECK(f2.topRows(9) == Matrix<float, 9, 20>::Zero());
                CHECK(f2.bottomRows(8) == Matrix<float, 8, 20>::Zero());
                CHECK(f2.leftCols(9) == Matrix<float, 20, 9>::Zero());
                CHECK(f2.rightCols(8) == Matrix<float, 20, 8>::Zero());
                CHECK(f2.block(9, 9, 3, 3) == Matrix<float, 3, 3>::Constant(1000.0 / 9.0));
            };
        };
        WHEN("Apply smoothing for default=6 iterations")
        {
            smooth(f2);
            THEN("Check smoothing result")
            {
                CHECK(f2.topRows(4) == Matrix<float, 4, 20>::Zero());
                CHECK(f2.bottomRows(3) == Matrix<float, 3, 20>::Zero());
                CHECK(f2.leftCols(4) == Matrix<float, 20, 4>::Zero());
                CHECK(f2.rightCols(3) == Matrix<float, 20, 3>::Zero());
                CHECK(f2.block(4, 4, 13, 13) != Matrix<float, 13, 13>::Zero());
                CHECK(f2.block(4, 4, 13, 13).minCoeff() > 0.0);
                CHECK(f2.maxCoeff() == f2(10, 10));
                CHECK(f2.rowwise().maxCoeff() == f2.col(10));
                CHECK(f2.colwise().maxCoeff() == f2.row(10));
                CHECK(f2.block(4, 4, 13, 13).rowwise().minCoeff() == f2.block(4, 4, 13, 1));
                CHECK(f2.block(4, 4, 13, 13).colwise().minCoeff() == f2.block(4, 4, 1, 13));
            };
        };
    };
}

TEST_CASE("Benchmark Field2 cubic()", "[Fields]")
{
    Eigen::Index n = 100;
    SECTION("Using T=float")
    {
        std::stringstream title;
        title << "f2<float> size=" << n << "x" << n;
        Field2<float> f2(n, n);
        f2.setRandom(n, n);
        BENCHMARK(title.str(), i) { return cubic(f2, Vec2f(std::fmodf(i * 123.123456f, 100.0f), std::fmodf(i * 321.124505f, 100.0f))); };
    };
    SECTION("Using T=double")
    {
        std::stringstream title;
        title << "f2<double> size=" << n << "x" << n;
        Field2<double> f2(n, n);
        f2.setRandom(n, n);
        BENCHMARK(title.str(), i) { return cubic(f2, Vec2d(std::fmod(i * 123.123456, 100.0), std::fmod(i * 321.124505, 100.0))); };
    };
};

TEST_CASE("Benchmark Field3 cubic()", "[Fields]")
{
    Eigen::Index n = 100;
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
            return cubic(f3, Vec3f(std::fmodf(i * 123.123456, 100.0), std::fmodf(i * 321.124505f, 100.0f),
                                   std::fmodf(i * 987.654333f, 100.0f)));
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
        { return cubic(f3, Vec3d(std::fmod(i * 123.123456, 100.0), std::fmod(i * 321.124505, 100.0), std::fmod(i * 987.654333, 100.0))); };
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
