#include <catch2/catch.hpp>

#include "libslic3r/Point.hpp"
// #include "libslic3r/BoundingBox.hpp"
// #include "libslic3r/Geometry.hpp"

#include <random>
#include "libslic3r/SVG.hpp"

using namespace Slic3r;
using namespace Catch::Matchers;

// See https://stackoverflow.com/a/79336425/5720535 for the pattern used here.
template<typename T, typename T2, typename En = void> struct CoordTraits;

// A template to encapsulte the coordinate types and standard operations on them.
template<typename _T, typename _T2, typename En = void> struct CoordTraitsBase
{
    using Coord                            = CoordTraits<_T, _T2, En>;
    using T                                = _T;
    using T2                               = _T2;
    using V                                = Vec<2, T>;
    using V2                               = Vec<2, T2>;
    using Vs                               = typename std::vector<V>;
    constexpr static double ScalingFactor  = 1e-6;
    constexpr static double ScalingFactor2 = Coord::ScalingFactor * Coord::ScalingFactor;
    // We limit Tmax to int32_t max because that's the "squared" limit we target.
    constexpr static T  Tmax  = std::numeric_limits<int32_t>::max();
    constexpr static T2 T2max = Coord::sqr(Tmax);

    constexpr static inline T  scale(const double v) { return static_cast<T>(v / Coord::ScalingFactor); }
    constexpr static inline T2 scale2(const double v) { return static_cast<T2>(v / Coord::ScalingFactor2); }

    constexpr static inline double unscale(const T v) { return v * Coord::ScalingFactor; }
    constexpr static inline double unscale2(const T2 v) { return v * Coord::ScalingFactor2; }

    constexpr static inline T to_T(const T2 v) { return static_cast<T>(v); }
    constexpr static inline V to_T(const V2& v) { return v.template cast<T>(); }

    constexpr static inline T2 to_T2(const T v) { return static_cast<T2>(v); }
    constexpr static inline V2 to_T2(const V& v) { return v.template cast<T2>(); }

    constexpr static inline T sqrt(const T2 v) { return Coord::to_T(static_cast<T2>(std::sqrt(v))); }
    constexpr static inline V sqrt(const V2& v) { return Coord::to_T(v.cwiseSqrt()); }

    constexpr static inline T2 sqr(const T v) { return Coord::to_T2(v) * Coord::to_T2(v); }
    constexpr static inline V2 sqr(const V& v) { return Coord::to_T2(v).cwiseAbs2(); }

    constexpr static inline T2 length2(const V& v) { return Coord::to_T2(v).squaredNorm(); }
    constexpr static inline T  length(const V& v) { return Coord::to_T(Coord::to_T2(v).norm()); };
};

template<typename _T, typename _T2, typename En> struct CoordTraits : public CoordTraitsBase<_T, _T2, En>
{};

// A Special variant using a scaled int32_t for T2 to avoid overflow.
template<typename _T> struct CoordTraits<_T, int32_t> : public CoordTraitsBase<_T, int32_t>
{
    using T                                = _T;
    using T2                               = int32_t;
    using V                                = Vec<2, T>;
    using V2                               = Vec<2, T2>;
    constexpr static int    T2Shift        = 16; // shift half the bits so sqr() doesn't overflow.
    constexpr static int    T2Mult         = 1 << T2Shift;
    constexpr static double ScalingFactor  = 1e-6;
    constexpr static double ScalingFactor2 = (ScalingFactor * T2Mult) * (ScalingFactor * T2Mult);

    constexpr static inline T to_T(const T2 v) { return static_cast<T>(v) << T2Shift; }
    constexpr static inline V to_T(const V2& v) { return V(v.coeff(0) << T2Shift, v.coeff(1) << T2Shift); }

    constexpr static inline T2 to_T2(const T v) { return static_cast<T2>(v >> T2Shift); }
    constexpr static inline V2 to_T2(const V& v) { return V2(v.coeff(0) >> T2Shift, v.coeff(1) >> T2Shift); }
};

// Get the result of min(v.norm(), len) but short-circuit evaluating
// the length-squared if len is clearly shorter.
//
// This includes multiple different implementations A through D for benchmarking.
template<typename Coord, class En = void> struct MinLength_X
{
    using T  = typename Coord::T;
    using T2 = typename Coord::T2;
    using V  = typename Coord::V;
    using V2 = typename Coord::V2;

    // Just directly get length.
    constexpr static inline T min_length_minlen_none(const V& v, const T len)
    {
        T l = Coord::length(v);
        return std::min(l, len);
    }

    // Just directly get length2 and check against sqr(len).
    constexpr static inline T min_length_len2_noVal(const V& v, const T len)
    {
        T2 l2 = Coord::length2(v);
        if (l2 >= Coord::sqr(len))
            return len;
        return Coord::sqrt(l2);
    }

    // Convert vector square to check our length2() overheads.
    constexpr static inline T min_length_vect2_noVal(const V& v, const T len)
    {
        T2 l2 = v.template cast<T2>().squaredNorm();
        if (l2 >= Coord::sqr(len))
            return len;
        return Coord::sqrt(l2);
    }

    // Convert coords and do all checks.
    constexpr static inline T min_length_coef2_all(const V& v, const T len)
    {
        T  x, y;
        T2 l2;
        if (((x = v.coeff(0)) > len) || (x < -len) || ((y = v.coeff(1)) > len) || (y < -len) ||
            ((l2 = Coord::sqr(x) + Coord::sqr(y)) >= Coord::sqr(len)))
            return len;
        return Coord::sqrt(l2);
    }

    // Convert coords and do all checks.
    constexpr static inline T min_length_len2_all(const V& v, const T len)
    {
        T  x, y;
        T2 l2;
        if (((x = v.coeff(0)) > len) || (x < -len) || ((y = v.coeff(1)) > len) || (y < -len) ||
            ((l2 = Coord::length2(v)) >= Coord::sqr(len)))
            return len;
        return Coord::sqrt(l2);
    }

    // Sqr coordinates and skip unsquared checks.
    constexpr static inline T min_length_len_all(const V& v, const T len)
    {
        T x, y, l;
        if (((x = v.coeff(0)) > len) || (x < -len) || ((y = v.coeff(1)) > len) || (y < -len) || ((l = Coord::length(v)) >= len))
            return len;
        return l;
    }

}; // Minlength_X

// Get the result of min(v.squaredNorm(), len2) but short-circuit evaluating
// the length-squared if len2 is clearly shorter.
//
// This includes multiple different implementations A through D for benchmarking.
template<typename Coord, class En = void> struct MinLength2_X
{
    using T  = typename Coord::T;
    using T2 = typename Coord::T2;
    using V  = typename Coord::V;
    using V2 = typename Coord::V2;

    // Just directly get length2.
    constexpr static inline T2 min_length2_minlen2_none(const V& v, const T2 len2)
    {
        T2 l2 = Coord::length2(v);
        return std::min(l2, len2);
    }

    // Convert coords and do all checks.
    constexpr static inline T2 min_length2_coefT2_all(const V& v, const T2 len2)
    {
        T2 x2, y2, l2;
        if (((x2 = Coord::to_T2(v.coeff(0))) > len2) || (x2 < -len2) || ((y2 = Coord::to_T2(v.coeff(1))) > len2) || (y2 < -len2) ||
            ((x2 *= x2) > len2) || ((y2 *= y2) > len2) || ((l2 = x2 + y2) >= len2))
            return len2;
        return l2;
    }

    // Convert V then do all checks.
    constexpr static inline T2 min_length2_vectT2_all(const V& v, const T2 len2)
    {
        V2 v2 = Coord::to_T2(v);
        T2 x2, y2, l2;
        if (((x2 = v2.coeff(0)) > len2) || (x2 < -len2) || ((y2 = v2.coeff(1)) > len2) || (y2 < -len2) || ((x2 *= x2) > len2) ||
            ((y2 *= y2) > len2) || ((l2 = x2 + y2) >= len2))
            return len2;
        return l2;
    }

    // Sqr coordinates and skip unsquared checks.
    constexpr static inline T2 min_length2_coefSqr_noVal(const V& v, const T2 len2)
    {
        T2 x2, y2, l2;
        if (((x2 = Coord::sqr(v.coeff(0))) > len2) || ((y2 = Coord::sqr(v.coeff(1))) > len2) || ((l2 = x2 + y2) >= len2))
            return len2;
        return l2;
    }

    // Sqr V and skip unsquared checks.
    constexpr static inline T2 min_length2_vectSqr_noVal(const V& v, const T2 len2)
    {
        V2 v2 = Coord::sqr(v);
        T2 l2;
        if ((v2.coeff(0) > len2) || (v2.coeff(1) > len2) || ((l2 = v2.sum()) >= len2))
            return len2;
        return l2;
    }

    // Convert coordinates and skip squared checks.
    constexpr static inline T2 min_length2_coefT2_noSqr(const V& v, const T2 len2)
    {
        T2 x2, y2, l2;
        if (((x2 = Coord::to_T2(v.coeff(0))) > len2) || (x2 < -len2) || ((y2 = Coord::to_T2(v.coeff(1))) > len2) || (y2 < -len2) ||
            ((l2 = x2 * x2 + y2 * y2) >= len2))
            return len2;
        return l2;
    }

    // Convert V and skip squared checks.
    constexpr static inline T2 min_length2_vectT2_noSqr(const V& v, const T2 len2)
    {
        V2 v2 = Coord::to_T2(v);
        T2 x2, y2, l2;
        if (((x2 = v2.coeff(0)) > len2) || (x2 < -len2) || ((y2 = v2.coeff(1)) > len2) || (y2 < -len2) || ((l2 = v2.squaredNorm()) >= len2))
            return len2;
        return l2;
    }

}; // MinLength2_X

template<typename T> using Data  = std::vector<Vec<2, T>>;
using Data32                     = Data<int32_t>;
using Data64                     = Data<int64_t>;
template<typename T> using Datas = std::array<Data<T>, 100>;
using Data32s                    = Datas<int32_t>;
using Data64s                    = Datas<int64_t>;

static Data32s gen_points(size_t point_cnt, unsigned seed = 7)
{
    Data32s                                data;
    std::mt19937_64                        rg(seed);
    std::uniform_int_distribution<int32_t> dist(-300 * 1e6, 300 * 1e6);
    for (int i = 0; i < data.size(); i++) {
        Data32 datai;
        datai.reserve(point_cnt);
        for (size_t j = 0; j < point_cnt; j++) {
            int32_t x = dist(rg), y = dist(rg);
            datai.emplace_back(x, y);
        }
        data[i] = datai;
    }
    return data;
}

static Data64s copy_points(const Data32s& data32)
{
    Data64s data;
    for (int i = 0; i < data.size(); i++) {
        Data64 datai;
        datai.reserve(data32[i].size());
        for (auto& v : data32[i]) {
            datai.emplace_back(v.template cast<int64_t>());
        }
        data[i] = datai;
    }
    return data;
}

static Data32s data32 = gen_points(10000);
static Data64s data64 = copy_points(data32);

template<typename T> Datas<T> get_data();
template<> Datas<int32_t>     get_data<int32_t>() { return data32; }
template<> Datas<int64_t>     get_data<int64_t>() { return data64; }

template<typename Coord> struct CoordTests
{
    using T  = typename Coord::T;
    using T2 = typename Coord::T2;
    using V  = typename Coord::V;
    using V2 = typename Coord::V2;
    using Vs = typename Coord::Vs;

    using Len  = MinLength_X<Coord>;
    using Len2 = MinLength2_X<Coord>;

    // static const Vs data = Setup::gen_points(1000);
    constexpr static T      t_max  = Coord::Tmax;
    constexpr static T2     t2_max = Coord::T2max;
    constexpr static double d_max  = Coord::unscale(Coord::Tmax);
    constexpr static double d2_max = Coord::unscale2(Coord::T2max);
    constexpr static T      t_1m   = Coord::scale(1000.0);
    constexpr static T2     t2_1m  = Coord::sqr(t_1m);

    static void checkcoords()
    {
        CHECK(Coord::scale(1.0) == static_cast<T>(1 / Coord::ScalingFactor));
        CHECK(Coord::scale2(1.0) == static_cast<T2>(1 / Coord::ScalingFactor2));
        CHECK(Coord::unscale(1) == Coord::ScalingFactor);
        CHECK(Coord::unscale2(1) == Coord::ScalingFactor2);

        CHECK(Coord::scale(Coord::unscale(1)) == 1);
        CHECK(Coord::scale2(Coord::unscale2(1)) == 1);
        CHECK(Coord::scale(Coord::unscale(t_1m)) == t_1m);
        CHECK(Coord::scale2(Coord::unscale2(t2_1m)) == t2_1m);
        CHECK(Coord::scale(Coord::unscale(t_max)) == t_max);
        CHECK(Coord::scale2(Coord::unscale2(t2_max)) == t2_max);

        CHECK(Coord::unscale(Coord::scale(1.0)) == 1.0);
        CHECK(Coord::unscale2(Coord::scale2(1.0)) == 1.0);
        CHECK(Coord::unscale(Coord::scale(1e3)) == 1e3);
        CHECK(Coord::unscale2(Coord::scale2(1e6)) == 1e6);
        CHECK(Coord::unscale(Coord::scale(d_max)) == d_max);
        CHECK(Coord::unscale2(Coord::scale2(d2_max)) == d2_max);

        CHECK(Coord::sqr(T(1)) == Coord::to_T2(T(1)) * Coord::to_T2(T(1)));
        CHECK(Coord::sqr(T(2)) == Coord::to_T2(T(2)) * Coord::to_T2(T(2)));
        CHECK(Coord::sqr(t_1m) == Coord::to_T2(t_1m) * Coord::to_T2(t_1m));
        CHECK(Coord::sqr(t_max) == Coord::to_T2(t_max) * Coord::to_T2(t_max));

        CHECK(Coord::sqrt(T2(1)) == Coord::to_T(T2(1)));
        CHECK(Coord::sqrt(T2(4)) == Coord::to_T(T2(2)));
        CHECK(Coord::sqrt(Coord::sqr(T(0x1230000))) == T(0x1230000));
        CHECK(Coord::sqrt(Coord::sqr(t_1m)) == t_1m);
        CHECK(Coord::sqrt(Coord::sqr(t_max)) == t_max);
    } // checkcoords

    // benchmark runner for min_length() tests.
    template<T (*F)(const V&, T)> static T runlen(const Vs& data)
    {
        T m = t_max;
        for (auto v : data) {
            m = F(v, m);
        }
        return m;
    };

    // benchmark runner for min_length2() tests.
    template<T2 (*F)(const V&, T2)> static T2 runlen2(const Vs& data)
    {
        T2 m = t2_max;
        for (auto v : data) {
            m = F(v, m);
        }
        return m;
    };

    static void runbench_length()
    {
        static const auto data = get_data<T>();
        T                 la;
        BENCHMARK("get_min_length_minlen_none", i) { return la = runlen<Len::min_length_minlen_none>(data[i % 100]); };
        T lb;
        BENCHMARK("get_min_length_len2_noVal", i) { return lb = runlen<Len::min_length_len2_noVal>(data[i % 100]); };
        T lc;
        BENCHMARK("get_min_length_vect2_noVal", i) { return lc = runlen<Len::min_length_vect2_noVal>(data[i % 100]); };
        T ld;
        BENCHMARK("get_min_length_coef2_all", i) { return ld = runlen<Len::min_length_coef2_all>(data[i % 100]); };
        T le;
        BENCHMARK("get_min_length_len2_all", i) { return le = runlen<Len::min_length_len2_all>(data[i % 100]); };
        T lf;
        BENCHMARK("get_min_length_len_all", i) { return lf = runlen<Len::min_length_len_all>(data[i % 100]); };
        CHECK(Coord::unscale(la) > 0.0);   // should be > 0.0mm
        CHECK(Coord::unscale(la) < 100.0); // should be less than 100mm
        CHECK(lb == la);
        CHECK(lc == la);
        CHECK(ld == la);
        CHECK(le == la);
        CHECK(lf == la);
    } // runbench_length()

    static void runbench_length2()
    {
        static const auto data = get_data<T>();
        T2                l2a;
        BENCHMARK("get_min_length2_minlen2_none", i) { return l2a = runlen2<Len2::min_length2_minlen2_none>(data[i % 100]); };
        T2 l2b;
        BENCHMARK("get_min_length2_coefT2_all", i) { return l2b = runlen2<Len2::min_length2_coefT2_all>(data[i % 100]); };
        T2 l2c;
        BENCHMARK("get_min_length2_vectT2_all", i) { return l2c = runlen2<Len2::min_length2_vectT2_all>(data[i % 100]); };
        T2 l2d;
        BENCHMARK("get_min_length2_coefSqr_noVal", i) { return l2d = runlen2<Len2::min_length2_coefSqr_noVal>(data[i % 100]); };
        T2 l2e;
        BENCHMARK("get_min_length2_vectSqr_noVal", i) { return l2e = runlen2<Len2::min_length2_vectSqr_noVal>(data[i % 100]); };
        T2 l2f;
        BENCHMARK("get_min_length2_coefT2_noSqr", i) { return l2f = runlen2<Len2::min_length2_coefT2_noSqr>(data[i % 100]); };
        T2 l2g;
        BENCHMARK("get_min_length2_vectT2_noSqr", i) { return l2g = runlen2<Len2::min_length2_vectT2_noSqr>(data[i % 100]); };
        CHECK(Coord::unscale2(l2a) > 0.0);           // should be > 0.0mm^2
        CHECK(Coord::unscale2(l2a) < 100.0 * 100.0); // should be less than 100mm^2
        CHECK(Coord::unscale(Coord::sqrt(l2a)) < 100.0);
        CHECK(l2b == l2a);
        CHECK(l2c == l2a);
        CHECK(l2d == l2a);
        CHECK(l2e == l2a);
        CHECK(l2f == l2a);
        CHECK(l2g == l2a);
    } // runbench_length2()

}; // RunBench

TEMPLATE_TEST_CASE("Benchmark minlength implementations using T=", "[MinLength2]", int32_t, int64_t)
{
    using T = TestType;

    SECTION("Using T2=int32_T", "This requires scaling to avoid overflow")
    {
        using T2 = int32_t;
        CoordTests<CoordTraits<T, T2>>::runbench_length();
    };
    SECTION("Using T2=int64_T", "This uses integers big enough to avoid overflow")
    {
        using T2 = int64_t;
        CoordTests<CoordTraits<T, T2>>::runbench_length();
    };
    SECTION("Using T2=double", "This uses doubles to avoid overflow")
    {
        using T2 = double;
        CoordTests<CoordTraits<T, T2>>::runbench_length();
    };
};

TEMPLATE_TEST_CASE("Benchmark minlength2 implementations T=", "[MinLength2]", int32_t, int64_t)
{
    using T = TestType;

    SECTION("Using T2=int32_T", "This requires scaling to avoid overflow")
    {
        using T2 = int32_t;
        CoordTests<CoordTraits<T, T2>>::runbench_length2();
    };
    SECTION("Using T2=int64_T", "This uses integers big enough to avoid overflow")
    {
        using T2 = int64_t;
        CoordTests<CoordTraits<T, T2>>::runbench_length2();
    };
    SECTION("Using T2=double", "This uses doubles to avoid overflow")
    {
        using T2 = double;
        CoordTests<CoordTraits<T, T2>>::runbench_length2();
    };
};

TEMPLATE_TEST_CASE("Test Coord implementation using T=", "[MinLength2]", int32_t, int64_t)
{
    using T = TestType;

    SECTION("Using T2=int32_T", "This requires scaling to avoid overflow")
    {
        using T2 = int32_t;
        CoordTests<CoordTraits<T, T2>>::checkcoords();
    };
    SECTION("Using T2=int64_T", "This uses integers big enough to avoid overflow")
    {
        using T2 = int64_t;
        CoordTests<CoordTraits<T, T2>>::checkcoords();
    };
    SECTION("Using T2=double", "This uses doubles to avoid overflow")
    {
        using T2 = double;
        CoordTests<CoordTraits<T, T2>>::checkcoords();
    };
};
