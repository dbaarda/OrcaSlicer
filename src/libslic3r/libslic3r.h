#ifndef _libslic3r_h_
#define _libslic3r_h_

#include "libslic3r_version.h"
#define SLIC3R_APP_FULL_NAME "Orca Slicer"
#define GCODEVIEWER_APP_NAME "OrcaSlicer G-code Viewer"
#define GCODEVIEWER_APP_KEY "OrcaSlicerGcodeViewer"
#define GCODEVIEWER_BUILD_ID std::string("OrcaSlicer G-code Viewer-") + std::string(SLIC3R_VERSION) + std::string("-RC")

// this needs to be included early for MSVC (listing it in Build.PL is not enough)
#include <memory>
#include <array>
#include <algorithm>
#include <ostream>
#include <iostream>
#include <math.h>
#include <queue>
#include <sstream>
#include <cstdio>
#include <stdint.h>
#include <stdarg.h>
#include <vector>
#include <cassert>
#include <cmath>
#include <type_traits>
#include <optional>

#ifdef _WIN32
// On MSVC, std::deque degenerates to a list of pointers, which defeats its purpose of reducing allocator load and memory fragmentation.
// https://github.com/microsoft/STL/issues/147#issuecomment-1090148740
// Thus it is recommended to use boost::container::deque instead.
#include <boost/container/deque.hpp>
#endif // _WIN32

#include "Technologies.hpp"
#include "Semver.hpp"

// These a the default types used for scaled and unscaled coordinates and their scaling factors.
//
// Using scaled integer types for coordinates is generally faster for most operations than using floating point types. However the size and
// scaling used are critical for sufficient range and resolution, and to avoid overflows during calculations.
//
// It is relatively easy to get sufficient range and resolution with a relatively small integer, however, calculating lengths or areas
// requires squaring or multiplying lengths, which requires much larger integers to avoid overflowing. For example using int32_t with a
// scaling factor of 1e-6 comfortably fits within +-2.1m with nice 1nm resolution, but would overflow when squaring anything outside
// +-0.046mm (sqrt(2^31) * 1e-6). Even using 1e-3 scaling for a tolerable 1um resolution can only square lengths within +-46.34mm without
// overflowing.
//
// In theory it is possible to work around this with careful use of down-scaling or casting to other types to avoid overflows when
// performing these kinds of operations. However that quickly becomes very complicated and impractical.
#if 0
// Saves around 32% RAM after slicing step, 6.7% after G-code export (tested on PrusaSlicer 2.2.0 final).
//
// Using int32_t with scaling factor 1e-2 fits within +-2.1e+4m at 0.01mm resolution, and anything within +-463mm will not overflow when squared.
// This should be OK for printers upto 460mm wide. However the resolution is pretty coarse and may cause noticeable artifacts. Note that approximate
// comparisions use EPSILON=1e-4 which scales to SCALED_EPSILON=0, reducing the approximate comparisons to exact ones for scaled coordinates. We can
// use float for unscaled coordinates without loss of resolution when using them for scaled values.
using coord_t = int32_t;
using coordf_t = float;
static constexpr coordf_t SCALING_FACTOR = 1e-2;
#else
// FIXME At least FillRectilinear2 and std::boost Voronoi require coord_t to be 32bit.
//
//  Using int64_t with scaling factor 1e-6 fits within +-9.2e+9m at 1nm resolution, and anything within +-3.0m will not overflow when
//  squared. This should be OK for very large printers upto 3m wide. If you really do have a printer larger than 3m, using 1e-5 for 10nm
//  resolution should be fine upto +-30.3m. We use double for unscaled coordinates to avoid loss of resolution when using them for scaled
//  values.
using coord_t                            = int64_t;
using coordf_t                           = double;
static constexpr coordf_t SCALING_FACTOR = 1e-6;
#endif
static constexpr coordf_t INV_SCALING_FACTOR = 1.0 / SCALING_FACTOR;

// This is the default type used for unscaled coordinates.

namespace Slic3r {

// A meta-predicate which is true for integers wider than or equal to coord_t.
template<class I>
struct is_scaled_coord
{
    static const constexpr bool value = std::is_integral<I>::value &&
                                        std::numeric_limits<I>::digits >= std::numeric_limits<coord_t>::digits;
};

// Meta predicates for floating, 'scaled coord' and generic arithmetic types
// Can be used to restrict templates to work for only the specified set of types.
// parameter T is the type we want to restrict
// parameter O (Optional defaults to T) is the type that the whole expression
// will be evaluated to.
// e.g. template<class T> FloatingOnly<T, bool> is_nan(T val);
// The whole template will be defined only for floating point types and the
// return type will be bool.
// For more info how to use, see docs for std::enable_if
//
template<class T, class O = T>
using FloatingOnly = std::enable_if_t<std::is_floating_point<T>::value, O>;

template<class T, class O = T>
using ScaledCoordOnly = std::enable_if_t<is_scaled_coord<T>::value, O>;

template<class T, class O = T>
using IntegerOnly = std::enable_if_t<std::is_integral<T>::value, O>;

template<class T, class O = T>
using ArithmeticOnly = std::enable_if_t<std::is_arithmetic<T>::value, O>;

template<class T, class O = T>
using IteratorOnly = std::enable_if_t<!std::is_same_v<typename std::iterator_traits<T>::value_type, void>, O>;

template<class From, class To>
using IsConvertible = std::enable_if_t<std::is_convertible_v<From, To>>;

// This gets the promoted result type of adding two types. It correctly handles things like coersion of ConfigOptions into their value type
// for the promotions, unlike std::common_type_t which can down-convert them into a lower precision/range type to match the other type.
template<class T1, class T2>
using get_result_t = decltype(std::declval<T1>() + std::declval<T2>());

// Conversion from any convertable unscaled type into any arithmetic scaled type.
// The return type defaults to coord_t but can be explicitly specified.
template<typename Tout = coord_t, typename Tin, typename = IsConvertible<Tin, Tout>, typename = ArithmeticOnly<Tout>>
inline constexpr Tout scaled(const Tin& v) noexcept
{
    // For efficiency we want to down-convert the scaling factor to the cheapest type that can be multiplied by Tin to give a result with
    // sufficient range for the result with at least as much resolution as Tin. Tin might also need to be coerced into this type if it is
    // something like a ConfigOption. Note INV_SCALING_FACTOR can be converted exactly into an int, so the cheapest possible type is int32_t.
    using UnscaledType = get_result_t<Tin, int32_t>;
    return Tout(v * UnscaledType(INV_SCALING_FACTOR));
}

// Conversion from any convertable scaled type to floating point unscaled type.
// The return type defaults to coordf_t but can be explicitly specified.
template<typename Tout = coordf_t, typename Tin, typename = IsConvertible<Tin, Tout>, typename = FloatingOnly<Tout>>
inline constexpr Tout unscaled(const Tin& v) noexcept
{ return Tout(v) * Tout(SCALING_FACTOR); }

// These are older macro versions of these functions that used to return a double expression.
#define scale_(val) scaled<coordf_t>(val)
#define unscale_(val) unscaled(val)

// And what about this one?
template<typename T, typename Q>
inline constexpr T unscale(Q v)
{ return unscaled<T>(v); }

// Orca: maximum number of extruders is 64. For SEMM printers, it defines maximum filament number.
static constexpr size_t MAXIMUM_EXTRUDER_NUMBER = 64;

// Orca: how many filament slots syncing an AMS setup may create. This was derived from
// EnforcerBlockerType::ExtruderMax, but that cap now covers 32 paintable filaments, so the AMS
// limit is pinned here to keep sync behaving as it does for projects without mixed-color filaments.
static constexpr size_t MAXIMUM_AMS_SYNC_FILAMENT_NUMBER = 16;

// Orca: maximum line width is 5 times the nozzle diameter
static constexpr float MAX_LINE_WIDTH_MULTIPLIER = 5;

static constexpr double PI = 3.141592653589793238;
#define POLY_SIDE_COUNT 24 // for brim ear circle
static constexpr coordf_t RESOLUTION                     = 0.0125;
static constexpr coord_t SCALED_RESOLUTION               = scaled(RESOLUTION);
static constexpr coordf_t SPARSE_INFILL_RESOLUTION       = 0.04;
static constexpr coord_t SCALED_SPARSE_INFILL_RESOLUTION = scaled(SPARSE_INFILL_RESOLUTION);

static constexpr coordf_t SUPPORT_RESOLUTION       = 0.0375;
static constexpr coord_t SCALED_SUPPORT_RESOLUTION = scaled(SUPPORT_RESOLUTION);
// Maximum perimeter length for the loop to apply the small perimeter speed.
#define SMALL_PERIMETER_LENGTH(LENGTH) scaled((LENGTH) * 2 * PI)
static constexpr coordf_t INSET_OVERLAP_TOLERANCE = 0.4;
// 3mm ring around the top / bottom / bridging areas.
// FIXME This is quite a lot.
static constexpr coordf_t EXTERNAL_INFILL_MARGIN = 3.;
static constexpr coordf_t BRIDGE_INFILL_MARGIN   = 1.;
static constexpr coordf_t WIPE_TOWER_MARGIN      = 1.;

// FIXME This epsilon value is used for many non-related purposes:
//  For a threshold of a squared Euclidean distance,
//  for a trheshold in a difference of radians,
//  for a threshold of a cross product of two non-normalized vectors etc.
static constexpr coordf_t EPSILON       = 1e-4;
static constexpr coord_t SCALED_EPSILON = scaled(EPSILON);
// A convenient templated EPSILON that is scaled or not depending on the type.
template<typename T>
static constexpr T COORD_EPSILON = is_scaled_coord<T>::value ? SCALED_EPSILON : EPSILON;

#ifndef UNUSED
#define UNUSED(x) (void) (x)
#endif /* UNUSED */

// BBS: some global const config which user can not change, but developer can
static constexpr bool g_config_support_sharp_tails                = true;
static constexpr bool g_config_remove_small_overhangs             = true;
static constexpr float g_config_tree_support_collision_resolution = 0.2;

// Write slices as SVG images into out directory during the 2D processing of the slices.
// #define SLIC3R_DEBUG_SLICE_PROCESSING

extern Semver SEMVER;

// On MSVC, std::deque degenerates to a list of pointers, which defeats its purpose of reducing allocator load and memory fragmentation.
template<class T, class Allocator = std::allocator<T>>
using deque =
#ifdef _WIN32
    // Use boost implementation, which allocates blocks of 512 bytes instead of blocks of 8 bytes.
    boost::container::deque<T, Allocator>;
#else  // _WIN32
    std::deque<T, Allocator>;
#endif // _WIN32

enum Axis {
    X = 0,
    Y,
    Z,
    E,
    F,
    // BBS: add I, J, P axis
    I,
    J,
    P,
    NUM_AXES,
    // For the GCodeReader to mark a parsed axis, which is not in "XYZEF", it was parsed correctly.
    UNKNOWN_AXIS = NUM_AXES,
    NUM_AXES_WITH_UNKNOWN,
};

template<typename T, typename Alloc, typename Alloc2>
inline void append(std::vector<T, Alloc>& dest, const std::vector<T, Alloc2>& src)
{
    if (dest.empty())
        dest = src;
    else
        dest.insert(dest.end(), src.begin(), src.end());
}

template<typename T, typename Alloc>
inline void append(std::vector<T, Alloc>& dest, std::vector<T, Alloc>&& src)
{
    if (dest.empty())
        dest = std::move(src);
    else {
        dest.reserve(dest.size() + src.size());
        std::move(std::begin(src), std::end(src), std::back_inserter(dest));
    }
    src.clear();
    src.shrink_to_fit();
}

template<class T, class... Args> // Arbitrary allocator can be used
void clear_and_shrink(std::vector<T, Args...>& vec)
{
    // shrink_to_fit does not garantee the release of memory nor does it clear()
    std::vector<T, Args...> tmp;
    vec.swap(tmp);
    assert(vec.capacity() == 0);
}

// Append the source in reverse.
template<typename T>
inline void append_reversed(std::vector<T>& dest, const std::vector<T>& src)
{
    if (dest.empty())
        dest = src;
    else
        dest.insert(dest.end(), src.rbegin(), src.rend());
}

// Append the source in reverse.
template<typename T>
inline void append_reversed(std::vector<T>& dest, std::vector<T>&& src)
{
    if (dest.empty())
        dest = std::move(src);
    else {
        dest.reserve(dest.size() + src.size());
        std::move(std::rbegin(src), std::rend(src), std::back_inserter(dest));
    }
    src.clear();
    src.shrink_to_fit();
}

// Casting an std::vector<> from one type to another type without warnings about a loss of accuracy.
template<typename T_TO, typename T_FROM>
std::vector<T_TO> cast(const std::vector<T_FROM>& src)
{
    std::vector<T_TO> dst;
    dst.reserve(src.size());
    for (const T_FROM& a : src)
        dst.emplace_back((T_TO) a);
    return dst;
}

template<typename T>
inline void remove_nulls(std::vector<T*>& vec)
{
    vec.erase(std::remove_if(vec.begin(), vec.end(), [](const T* ptr) { return ptr == nullptr; }), vec.end());
}

template<typename T>
inline void sort_remove_duplicates(std::vector<T>& vec)
{
    std::sort(vec.begin(), vec.end());
    vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
}

// Older compilers do not provide a std::make_unique template. Provide a simple one.
template<typename T, typename... Args>
inline std::unique_ptr<T> make_unique(Args&&... args)
{ return std::unique_ptr<T>(new T(std::forward<Args>(args)...)); }

// Variant of std::lower_bound() with compare predicate, but without the key.
// This variant is very useful in case that the T type is large or it does not even have a public constructor.
template<class ForwardIt, class LowerThanKeyPredicate>
ForwardIt lower_bound_by_predicate(ForwardIt first, ForwardIt last, LowerThanKeyPredicate lower_than_key)
{
    ForwardIt it;
    typename std::iterator_traits<ForwardIt>::difference_type count, step;
    count = std::distance(first, last);

    while (count > 0) {
        it   = first;
        step = count / 2;
        std::advance(it, step);
        if (lower_than_key(*it)) {
            first = ++it;
            count -= step + 1;
        } else
            count = step;
    }
    return first;
}

// from https://en.cppreference.com/w/cpp/algorithm/lower_bound
template<class ForwardIt, class T, class Compare = std::less<>>
ForwardIt binary_find(ForwardIt first, ForwardIt last, const T& value, Compare comp = {})
{
    // Note: BOTH type T and the type after ForwardIt is dereferenced
    // must be implicitly convertible to BOTH Type1 and Type2, used in Compare.
    // This is stricter than lower_bound requirement (see above)

    first = std::lower_bound(first, last, value, comp);
    return first != last && !comp(value, *first) ? first : last;
}

// from https://en.cppreference.com/w/cpp/algorithm/lower_bound
template<class ForwardIt, class LowerThanKeyPredicate, class EqualToKeyPredicate>
ForwardIt binary_find_by_predicate(ForwardIt first, ForwardIt last, LowerThanKeyPredicate lower_thank_key, EqualToKeyPredicate equal_to_key)
{
    // Note: BOTH type T and the type after ForwardIt is dereferenced
    // must be implicitly convertible to BOTH Type1 and Type2, used in Compare.
    // This is stricter than lower_bound requirement (see above)

    first = lower_bound_by_predicate(first, last, lower_thank_key);
    return first != last && equal_to_key(*first) ? first : last;
}

template<typename ContainerType, typename ValueType>
inline bool contains(const ContainerType& c, const ValueType& v)
{ return std::find(c.begin(), c.end(), v) != c.end(); }
template<typename T>
inline bool contains(const std::initializer_list<T>& il, const T& v)
{ return std::find(il.begin(), il.end(), v) != il.end(); }

template<typename ContainerType, typename ValueType>
inline bool one_of(const ValueType& v, const ContainerType& c)
{ return contains(c, v); }
template<typename T>
inline bool one_of(const T& v, const std::initializer_list<T>& il)
{ return contains(il, v); }

template<typename T>
inline constexpr T sqr(T x)
{ return x * x; }

// Is value approximately zero?
template<typename Number, typename Precision = Number>
inline constexpr bool is_zero(const Number value, const Precision precision = COORD_EPSILON<Number>)
{
    // Note we use <= here so that precision=0 works.
    return std::abs(value) <= precision;
}

// Is value approximately equal to test?
template<typename Arg1, typename Arg2, typename Number = get_result_t<Arg1, Arg2>, typename Precision = Number>
inline constexpr bool is_approx(const Arg1 value, const Arg2 test, const Precision precision = COORD_EPSILON<Number>)
{ return is_zero(value - test, precision); }

// Is value approximately < test?
template<typename Arg1, typename Arg2, typename Number = get_result_t<Arg1, Arg2>, typename Precision = Number>
inline constexpr bool is_approx_lt(const Arg1 value, const Arg2 test, const Precision precision = COORD_EPSILON<Number>)
{
    // Note we use < here so that precision=0 works.
    return value < (test - precision);
}

// Is value approximately == test? For std::optional numbers.
// Will return true if both have no value, or false if only one has a value.
template<typename Arg1, typename Arg2, typename Number = get_result_t<Arg1, Arg2>, typename Precision = Number>
inline constexpr bool is_approx(const std::optional<Arg1>& value,
                                const std::optional<Arg2>& test,
                                const Precision precision = COORD_EPSILON<Number>)
{
    return (!value.has_value() && !test.has_value()) ||
           (value.has_value() && test.has_value() && is_approx<Number>(*value, *test, precision));
}

// Is value approximately < test? For std:optional numbers.
// Having any value will be considered greater than no value.
template<typename Arg1, typename Arg2, typename Number = get_result_t<Arg1, Arg2>, typename Precision = Number>
inline constexpr bool is_approx_lt(const std::optional<Arg1>& value,
                                   const std::optional<Arg2>& test,
                                   const Precision precision = COORD_EPSILON<Number>)
{
    return (!value.has_value() && test.has_value()) ||
           (value.has_value() && test.has_value() && is_approx_lt<Number>(*value, *test, precision));
}

// Is value approximately > test?
template<typename Arg1, typename Arg2, typename Number = get_result_t<Arg1, Arg2>, typename Precision = Number>
inline constexpr bool is_approx_gt(const Arg1 value, const Arg2 test, const Precision precision = COORD_EPSILON<Number>)
{ return is_approx_lt(test, value, precision); }

// Is value approximately <= test?
template<typename Arg1, typename Arg2, typename Number = get_result_t<Arg1, Arg2>, typename Precision = Number>
inline constexpr bool is_approx_le(const Arg1 value, const Arg2 test, const Precision precision = COORD_EPSILON<Number>)
{ return !is_approx_lt(test, value, precision); }

// Is value approximately >= test?
template<typename Arg1, typename Arg2, typename Number = get_result_t<Arg1, Arg2>, typename Precision = Number>
inline constexpr bool is_approx_ge(const Arg1 value, const Arg2 test, const Precision precision = COORD_EPSILON<Number>)
{ return !is_approx_lt(value, test, precision); }

// Linearly interpolate between a and b by ratio t.
template<typename T, typename Number>
inline constexpr T lerp(const T& a, const T& b, Number t)
{
    assert(is_approx_le(Number(0), t) && is_approx_le(t, Number(1)));
    return (Number(1) - t) * a + t * b;
}

template<class T, class I, class... Args> // Arbitrary allocator can be used
IntegerOnly<I, std::vector<T, Args...>> reserve_vector(I capacity)
{
    std::vector<T, Args...> ret;
    if (capacity > I(0))
        ret.reserve(size_t(capacity));

    return ret;
}

// Borrowed from C++20
template<class T>
using remove_cvref_t = std::remove_cv_t<std::remove_reference_t<T>>;

// A very simple range concept implementation with iterator-like objects.
// This should be replaced by std::ranges::subrange (C++20)
template<class It>
class Range
{
    It from, to;

public:
    // The class is ready for range based for loops.
    It begin() const { return from; }
    It end() const { return to; }

    // The iterator type can be obtained this way.
    using iterator   = It;
    using value_type = typename std::iterator_traits<It>::value_type;

    Range() = default;
    Range(It b, It e) : from(std::move(b)), to(std::move(e)) {}

    // Some useful container-like methods...
    inline size_t size() const { return std::distance(from, to); }
    inline bool empty() const { return from == to; }
};

template<class Cont>
auto range(Cont&& cont)
{ return Range{std::begin(cont), std::end(cont)}; }

template<class T, class = FloatingOnly<T>>
constexpr T NaN = std::numeric_limits<T>::quiet_NaN();

constexpr float NaNf  = NaN<float>;
constexpr double NaNd = NaN<double>;

// Rounding up.
// 1.5 is rounded to 2
// 1.49 is rounded to 1
// 0.5 is rounded to 1,
// 0.49 is rounded to 0
// -0.5 is rounded to 0,
// -0.51 is rounded to -1,
// -1.5 is rounded to -1.
// -1.51 is rounded to -2.
// If input is not a valid float (it is infinity NaN or if it does not fit)
// the float to int conversion produces a max int on Intel and +-max int on ARM.
template<typename I>
inline IntegerOnly<I, I> fast_round_up(double a)
{
    // Why does Java Math.round(0.49999999999999994) return 1?
    // https://stackoverflow.com/questions/9902968/why-does-math-round0-49999999999999994-return-1
    return a == 0.49999999999999994 ? I(0) : I(floor(a + 0.5));
}

template<class T>
using SamePair = std::pair<T, T>;

} // namespace Slic3r

#endif
