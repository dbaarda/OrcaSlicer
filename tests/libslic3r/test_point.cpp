#include <catch2/catch.hpp>

#include "libslic3r/Point.hpp"
// #include "libslic3r/BoundingBox.hpp"
// #include "libslic3r/Geometry.hpp"

#include <random>
#include "libslic3r/SVG.hpp"

using namespace Slic3r;
using namespace Catch::Matchers;

// See https://stackoverflow.com/a/79336425/5720535 for the pattern used here.
template<typename T, typename T2> struct CoordTraits;

// A template to encapsulte the coordinate types and standard operations on them.
template<typename _T, typename _T2, class En = void> struct CoordTraitsBase
{
    using Coord                            = CoordTraits<_T, _T2>;
    using T                                = _T;
    using T2                               = _T2;
    using V                                = Vec<2, T>;
    using V2                               = Vec<2, T2>;
    using Vs                               = typename std::vector<V>;
    constexpr static double ScalingFactor  = 1e-6;
    constexpr static double ScalingFactor2 = Coord::ScalingFactor * Coord::ScalingFactor;

    constexpr static inline T  scale(const double v) { return static_cast<T>(v / Coord::ScalingFactor); }
    constexpr static inline T2 scale2(const double v) { return static_cast<T2>(v / Coord::ScalingFactor2); }

    constexpr static inline double unscale(const T v) { return v * Coord::ScalingFactor; }
    constexpr static inline double unscale2(const T2 v) { return v * Coord::ScalingFactor2; }

    constexpr static inline T to_T(const T2 v) { return static_cast<T>(v); }

    constexpr static inline V to_T(const V2& v) { return v.template cast<T>(); }

    constexpr static inline T2 to_T2(const T v) { return static_cast<T2>(v); }

    constexpr static inline auto to_T2(const V& v) { return v.template cast<T2>(); }

    constexpr static inline T  sqrt(const T2 v) { return Coord::to_T(static_cast<T2>(std::sqrt(v))); }
    constexpr static inline T2 sqr(const T v)
    {
        T2 v2 = Coord::to_T2(v);
        return v2 * v2;
    }
    // constexpr static inline V2 sqr(const V& v) { return v.template cast<T2>.cwiseAbs2(); }
    constexpr static inline V2 sqr(const V& v) { return Coord::to_T2(v).cwiseAbs2(); }

    constexpr static inline T2 length2(const V& v) { return Coord::to_T2(v).squaredNorm(); }
  constexpr static inline T  length(const V& v) { return Coord::to_T(Coord::to_T2(v).norm()); };
};

template<typename _T, typename _T2> struct CoordTraits : public CoordTraitsBase<_T, _T2>
{};

// A specific instance using a scaled int32_t for T2 to avoid overflow.
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
    static inline V           to_T(const V2& v) { return v.template cast<T>() * T2Mult; }
    // constexpr static inline V to_T(const V2& v) { return V(to_T(v.coeff(0)), to_T(v.coeff(1))); }
    // static inline V to_T(const V2& v) { return v.template cast<T>().unaryExpr([](const T x){ return x<<16; }); }
    constexpr static inline T2 to_T2(const T v) { return static_cast<T2>(v >> T2Shift); }
    static inline V2           to_T2(const V& v) { return (v / T2Mult).template cast<T2>(); }
    // constexpr static inline V2 to_T2(const V& v) { return V2(to_T2(v.coeff(0)), to_T2(v.coeff(1))); }
    // static inline V2 to_T2(const V& v) { return v.unaryExpr([](const T x){ return x>>16; }).template cast<T2>(); }
    // constexpr static inline V2 sqr(const V& v) { return (v / T2Mult).template cast<T2>.cwiseAbs2(); }
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
        T2 l2;
        if ((l2 = Coord::length2(v)) >= Coord::sqr(len))
          return len;
      return Coord::sqrt(l2);
    }

    // Convert vector square to check our length2() overheads.
    constexpr static inline T min_length_vect2_noVal(const V& v, const T len)
    {
        T2 l2;
        if ((l2 = v.template cast<T2>().squaredNorm()) >= Coord::sqr(len))
          return len;
      return Coord::sqrt(l2);
    }

    // Convert coords and do all checks.
    constexpr static inline T min_length_coef2_all(const V& v, const T len)
    {
        T x, y;
        T2 l2;
        if (((x = v.coeff(0)) > len) || (x < -len) ||
            ((y = v.coeff(1)) > len) || (y < -len) ||
            ((l2 = Coord::sqr(x) + Coord::sqr(y)) >= Coord::sqr(len)))
            return len;
      return Coord::sqrt(l2);
    }

    // Convert coords and do all checks.
    constexpr static inline T min_length_len2_all(const V& v, const T len)
    {
        T x, y;
        T2 l2;
        if (((x = v.coeff(0)) > len) || (x < -len) ||
            ((y = v.coeff(1)) > len) || (y < -len) ||
            ((l2 = Coord::length2(v)) >= Coord::sqr(len)))
            return len;
      return Coord::sqrt(l2);
    }

    // Sqr coordinates and skip unsquared checks.
    constexpr static inline T min_length_len_all(const V& v, const T len)
    {
        T x, y, l;
        if (((x = v.coeff(0)) > len) || (x < -len) ||
            ((y = v.coeff(1)) > len) || (y < -len) ||
            ((l = Coord::length(v)) >= len))
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

template<typename T> using Data = std::vector<Vec<2, T>>;
using Data32                    = Data<int32_t>;
using Data64                    = Data<int64_t>;

static Data32 gen_points(size_t point_cnt)
{
    Data32 data;
    data.reserve(point_cnt);
    auto                                   seed = std::random_device{}();
    std::mt19937_64                        rg{seed};
    std::uniform_int_distribution<int32_t> dist(-1000 * 1e6, 1000 * 1e6);
    for (size_t i = 0; i < point_cnt; ++i) {
        int32_t x = dist(rg), y = dist(rg);
        data.emplace_back(x, y);
    }
    return data;
}
static Data64 copy_points(const Data32& data32)
{
    Data64 data;
    data.reserve(data32.size());
    for (auto& v : data32) {
        data.emplace_back(v.template cast<int64_t>());
    }
    return data;
}

static Data32 data32 = gen_points(1000);
static Data64 data64 = copy_points(data32);

template<typename T> Data<T> get_data();
template<> Data<int32_t>     get_data<int32_t>() { return data32; }
template<> Data<int64_t>     get_data<int64_t>() { return data64; }


template<typename Coord> struct CoordTests
{
    using T  = typename Coord::T;
    using T2 = typename Coord::T2;
    using V  = typename Coord::V;
    using V2 = typename Coord::V2;
    using Vs = typename Coord::Vs;

    using Len = MinLength_X<Coord>;
    using Len2 = MinLength2_X<Coord>;

    // static const Vs data = Setup::gen_points(1000);
    constexpr static T  t_max  = std::numeric_limits<T>::max();
    constexpr static T2 t2_max = std::numeric_limits<T2>::max();

    static void checkcoords()
    {
        CHECK(Coord::to_T(T2(4)) == T(4));
        CHECK(Coord::to_T2(T(4)) == T2(4));

        CHECK(Coord::scale(Coord::unscale(t_max)) == t_max);
        CHECK(Coord::scale2(Coord::unscale2(t2_max)) == t2_max);

        CHECK(Coord::scale(1.0) == static_cast<T>(1 / Coord::ScalingFactor));
        CHECK(Coord::scale2(1.0) == static_cast<T2>(1 / Coord::ScalingFactor2));

        CHECK(Coord::unscale(1) == Coord::ScalingFactor);
        CHECK(Coord::unscale2(1) == Coord::ScalingFactor2);

        CHECK(Coord::scale(t_max * Coord::ScalingFactor) == t_max);
        CHECK(Coord::unscale(t_max) == t_max * Coord::ScalingFactor);
        CHECK(Coord::scale2(t2_max * Coord::ScalingFactor2) == t2_max);
        CHECK(Coord::unscale2(t2_max) == t2_max * Coord::ScalingFactor2);

        CHECK(Coord::sqr(T(1)) == Coord::to_T2(T(1)) * Coord::to_T2(T(1)));
        CHECK(Coord::sqr(T(2)) == Coord::to_T2(T(2)) * Coord::to_T2(T(2)));
        CHECK(Coord::sqr(t_max) == Coord::to_T2(t_max) * Coord::to_T2(t_max));

        CHECK(Coord::sqrt(T2(1)) == Coord::to_T(T2(1)));
        CHECK(Coord::sqrt(T2(4)) == Coord::to_T(T2(2)));
        CHECK(Coord::sqrt(Coord::sqr(T(0x1230000))) == T(0x1230000));
        CHECK(Coord::sqrt(Coord::sqr(t_max)) == t_max);
    } // checkcoords

    // benchmark runner for min_length() tests.
    template<T (*F)(const V&, T)> static T runlen(const Vs& data)
    {
        T m = std::numeric_limits<T>::max();
        for (auto v : data) {
            m = F(v, m);
        }
        return m;
    };

    // benchmark runner for min_length2() tests.
    template<T2 (*F)(const V&, T2)> static T2 runlen2(const Vs& data)
    {
        T2 m = std::numeric_limits<T2>::max();
        for (auto v : data) {
            m = F(v, m);
        }
        return m;
    };

    static void runbench_length()
    {
        static const auto data = get_data<T>();
        T                la;
        BENCHMARK("get_min_length_minlen_none") { return la = runlen<Len::min_length_minlen_none>(data); };
        T lb;
        BENCHMARK("get_min_length_len2_noVal") { return lb = runlen<Len::min_length_len2_noVal>(data); };
        T lc;
        BENCHMARK("get_min_length_vect2_noVal") { return lc = runlen<Len::min_length_vect2_noVal>(data); };
        T ld;
        BENCHMARK("get_min_length_coef2_all") { return ld = runlen<Len::min_length_coef2_all>(data); };
        T le;
        BENCHMARK("get_min_length_len2_all") { return le = runlen<Len::min_length_len2_all>(data); };
        T lf;
        BENCHMARK("get_min_length_len_all") { return lf = runlen<Len::min_length_len_all>(data); };
        CHECK(Coord::unscale(la) > 0.0);           // should be > 0.0mm
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
        BENCHMARK("get_min_length2_minlen2_none") { return l2a = runlen2<Len2::min_length2_minlen2_none>(data); };
        T2 l2b;
        BENCHMARK("get_min_length2_coefT2_all") { return l2b = runlen2<Len2::min_length2_coefT2_all>(data); };
        T2 l2c;
        BENCHMARK("get_min_length2_vectT2_all") { return l2c = runlen2<Len2::min_length2_vectT2_all>(data); };
        T2 l2d;
        BENCHMARK("get_min_length2_coefSqr_noVal") { return l2d = runlen2<Len2::min_length2_coefSqr_noVal>(data); };
        T2 l2e;
        BENCHMARK("get_min_length2_vectSqr_noVal") { return l2e = runlen2<Len2::min_length2_vectSqr_noVal>(data); };
        T2 l2f;
        BENCHMARK("get_min_length2_coefT2_noSqr") { return l2f = runlen2<Len2::min_length2_coefT2_noSqr>(data); };
        T2 l2g;
        BENCHMARK("get_min_length2_vectT2_noSqr") { return l2g = runlen2<Len2::min_length2_vectT2_noSqr>(data); };
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

TEMPLATE_TEST_CASE("Benchmark minlength implementations", "[MinLength2]", int32_t, int64_t)
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

TEMPLATE_TEST_CASE("Benchmark minlength2 implementations", "[MinLength2]", int32_t, int64_t)
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

TEST_CASE("Test Coord implementation using T=int32_t", "[MinLength2]")
{
    using T = int32_t;

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
