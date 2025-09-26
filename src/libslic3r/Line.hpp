#ifndef slic3r_Line_hpp_
#define slic3r_Line_hpp_

#include "libslic3r.h"
#include "Point.hpp"

#include <type_traits>

namespace Slic3r {

class BoundingBox;
class Line;
class Line3;
class Linef3;
class Polyline;
class ThickLine;
typedef std::vector<Line> Lines;
typedef std::vector<Line3> Lines3;
typedef std::vector<ThickLine> ThickLines;

Linef3 transform(const Linef3& line, const Transform3d& t);

namespace line_alg {

template<class L, class En = void> struct Traits {
    static constexpr int Dim = L::Dim;
    using Scalar = typename L::Scalar;
    using Point = Vec<Dim, Scalar>;

    static Point& get_a(L &l) { return l.a; }
    static Point& get_b(L &l) { return l.b; }
    static const Point& get_a(const L &l) { return l.a; }
    static const Point& get_b(const L &l) { return l.b; }
};

template<class L> const constexpr int Dim = Traits<remove_cvref_t<L>>::Dim;
template<class L> using Scalar = typename Traits<remove_cvref_t<L>>::Scalar;
template<class L> using Point = typename Traits<remove_cvref_t<L>>::Point;

template<class L> auto get_a(L &&l) { return Traits<remove_cvref_t<L>>::get_a(l); }
template<class L> auto get_b(L &&l) { return Traits<remove_cvref_t<L>>::get_b(l); }

// Get the closest point on a line to a point.
template<class L>
const Point<L> closest_point(const L &line, const Point<L> &point)
{
    if (get_a(line) == get_b(line)) return get_a(line);
    const Vec<Dim<L>, double>  v  = (get_b(line) - get_a(line)).template cast<double>();
    const Vec<Dim<L>, double>  va = (point  - get_a(line)).template cast<double>();
    // Consider the line extending the segment, parameterized as a + t (b - a).
    // We find projection of this point onto the line.
    // It falls where t = [(this-a) . (b-a)] / |b-a|^2
    const double t = va.dot(v) / v.squaredNorm();
    if (t <= 0.0) {
        // beyond the 'a' end of the segment
        return get_a(line);
    } else if (t >= 1.0) {
        // beyond the 'b' end of the segment
        return get_b(line);
    }
    return get_a(line) + (t * v).template cast<Scalar<L>>();
}

// Get the point projected onto the infinite line.
// The projected_point could be beyond the 'a' and 'b' ends of the segment.
template<class L>
const Point<L> projected_point(const L &line, const Point<L> &point)
{
    if (get_a(line) == get_b(line)) return get_a(line);
    const Vec<Dim<L>, double>  v  = (get_b(line) - get_a(line)).template cast<double>();
    const Vec<Dim<L>, double>  va = (point  - get_a(line)).template cast<double>();
    const double t = va.dot(v) / v.squaredNorm();
    return get_a(line) + (t * v).template cast<Scalar<L>>();
}


// Minimum distance squared and nearest point from a line to point.
template<class L>
double distance_to_squared(const L &line, const Point<L> &point, Point<L> *nearest_point)
{
   *nearest_point = closest_point(line, point);
   return dist2(point, *nearest_point);
}

// Minimum distance squared from a line to a point.
template<class L>
double distance_to_squared(const L &line, const Point<L> &point)
{
    return dist2(point, closest_point(line, point));
}

// Minimum distance from a line to a point.
template<class L>
double distance_to(const L &line, const Point<L> &point)
{
    return std::sqrt(distance_to_squared(line, point));
}

// Returns a squared distance to the projected point on the infinite line.
// The projected_point could be beyond the 'a' and 'b' ends of the segment.
template<class L>
double projected_distance_squared(const L &line, const Point<L> &point, Point<L> *nearest_point)
{
   *nearest_point = projected_point(line, point);
   return dist2(point, *nearest_point);
}

// Returns a squared distance to the projected point on the infinite line
// The projected point could be beyond the 'a' and 'b' ends of the segment.
template<class L>
double projected_distance_squared(const L &line, const Point<L> &point)
{
   // If you don't need the point, this is cheaper.
    if (get_a(line) == get_b(line)) return dist2(get_a(line), point);
    const Vec<Dim<L>, double>  v  = (get_b(line) - get_a(line)).template cast<double>();
    const Vec<Dim<L>, double>  va = (point  - get_a(line)).template cast<double>();
    return sqr(cross2(v, va)) / v.squaredNorm();
}

// Returns a distance to the projected point on the infinite line.
// The projected point could be beyond the 'a' and 'b' ends of the segment.
template<class L>
double projected_distance(const L &line, const Point<L> &point)
{
   // If you don't need the point, this is cheaper.
    if (get_a(line) == get_b(line)) return dist2(get_a(line), point);
    const Vec<Dim<L>, double>  v  = (get_b(line) - get_a(line)).template cast<double>();
    const Vec<Dim<L>, double>  va = (point  - get_a(line)).template cast<double>();
    return std::abs(cross2(v, va)) / v.norm();
}

template<class L> bool intersection(const L &l1, const L &l2, Vec<Dim<L>, Scalar<L>> *intersection_pt)
{
    using Floating      = typename std::conditional<std::is_floating_point<Scalar<L>>::value, Scalar<L>, double>::type;
    using VecType       = const Vec<Dim<L>, Floating>;
    const VecType v1    = (l1.b - l1.a).template cast<Floating>();
    const VecType v2    = (l2.b - l2.a).template cast<Floating>();
    Floating      denom = cross2(v1, v2);
    if (std::abs(denom) < EPSILON)
#if 0
        // Lines are collinear. Return true if they are coincident (overlappign).
        return ! (std::abs(nume_a) < EPSILON && std::abs(nume_b) < EPSILON);
#else
        return false;
#endif
    const VecType v12 = (l1.a - l2.a).template cast<Floating>();
    Floating nume_a = cross2(v2, v12);
    Floating nume_b = cross2(v1, v12);
    Floating t1     = nume_a / denom;
    Floating t2     = nume_b / denom;
    if (t1 >= 0 && t1 <= 1.0f && t2 >= 0 && t2 <= 1.0f) {
        // Get the intersection point.
        (*intersection_pt) = (l1.a.template cast<Floating>() + t1 * v1).template cast<Scalar<L>>();
        return true;
    }
    return false; // not intersecting
}

} // namespace line_alg

class Line
{
public:
    Line() {}
    Line(const Point& _a, const Point& _b) : a(_a), b(_b) {}
    explicit operator Lines() const { Lines lines; lines.emplace_back(*this); return lines; }
    void   scale(double factor) { this->a *= factor; this->b *= factor; }
    void   translate(const Point &v) { this->a += v; this->b += v; }
    void   translate(double x, double y) { this->translate(Point(x, y)); }
    void   rotate(double angle, const Point &center) { this->a.rotate(angle, center); this->b.rotate(angle, center); }
    void   reverse() { std::swap(this->a, this->b); }
    double length() const { return a.distance_to(b); }
    Point  midpoint() const { return (this->a + this->b) / 2; }
    bool   intersection_infinite(const Line &other, Point* point) const;
    bool   operator==(const Line &rhs) const { return this->a == rhs.a && this->b == rhs.b; }
    Point  closest_point(const Point &point) const { return line_alg::closest_point(*this, point); };
    double distance_to_squared(const Point &point, Point *closest_point) const { return line_alg::distance_to_squared(*this, point, closest_point); }
    double distance_to_squared(const Point &point) const { return line_alg::distance_to_squared(*this, point); }
    double distance_to(const Point &point) const { return line_alg::distance_to(*this, point); }
    // Returns a distance to the closest point on the infinite line through a and b.
    // Closest point (and returned squared distance to this point) could be beyond the 'a' and 'b' ends of the segment.
    Point  projected_point(const Point &point) const { return line_alg::projected_point(*this, point); };
    double projected_distance_squared(const Point &point, Point *closest_point) const { return line_alg::projected_distance_squared(*this, point, closest_point); }
    double projected_distance_squared(const Point &point) const { return line_alg::projected_distance_squared(*this, point); }
    double projected_distance(const Point &point) { return line_alg::projected_distance(*this, point); }
    bool   parallel_to(double angle) const;
    bool   parallel_to(const Line& line) const;
    bool   perpendicular_to(double angle) const;
    bool   perpendicular_to(const Line& line) const;
    double atan2_() const { return atan2(this->b(1) - this->a(1), this->b(0) - this->a(0)); }
    double orientation() const;
    double direction() const;
    Vector vector() const { return this->b - this->a; }
    Vector normal() const { return Vector((this->b(1) - this->a(1)), -(this->b(0) - this->a(0))); }
    bool   intersection(const Line& line, Point* intersection) const;
    // Clip a line with a bounding box. Returns false if the line is completely outside of the bounding box.
    bool   clip_with_bbox(const BoundingBox &bbox);
    // Extend the line from both sides by an offset.
    void   extend(double offset);

    Point a;
    Point b;

    static const constexpr int Dim = 2;
    using Scalar = Point::Scalar;
};

class ThickLine : public Line
{
public:
    ThickLine() : a_width(0), b_width(0) {}
    ThickLine(const Point& a, const Point& b) : Line(a, b), a_width(0), b_width(0) {}
    ThickLine(const Point& a, const Point& b, double wa, double wb) : Line(a, b), a_width(wa), b_width(wb) {}

    double a_width, b_width;
};

class CurledLine : public Line
{
public:
    CurledLine() : curled_height(0.0f) {}
    CurledLine(const Point& a, const Point& b) : Line(a, b), curled_height(0.0f) {}
    CurledLine(const Point& a, const Point& b, float curled_height) : Line(a, b), curled_height(curled_height) {}

    float curled_height;
};

using CurledLines = std::vector<CurledLine>;

class Line3
{
public:
    Line3() : a(Vec3crd::Zero()), b(Vec3crd::Zero()) {}
    Line3(const Vec3crd& _a, const Vec3crd& _b) : a(_a), b(_b) {}

    double  length() const { return (this->a - this->b).cast<double>().norm(); }
    Vec3crd vector() const { return this->b - this->a; }

    Vec3crd a;
    Vec3crd b;

    static const constexpr int Dim = 3;
    using Scalar = Vec3crd::Scalar;
};

class Linef
{
public:
    Linef() : a(Vec2d::Zero()), b(Vec2d::Zero()) {}
    Linef(const Vec2d& _a, const Vec2d& _b) : a(_a), b(_b) {}

    Vec2d a;
    Vec2d b;

    static const constexpr int Dim = 2;
    using Scalar = Vec2d::Scalar;
};
using Linesf = std::vector<Linef>;

class Linef3
{
public:
    Linef3() : a(Vec3d::Zero()), b(Vec3d::Zero()) {}
    Linef3(const Vec3d& _a, const Vec3d& _b) : a(_a), b(_b) {}

    Vec3d   intersect_plane(double z) const;
    void    scale(double factor) { this->a *= factor; this->b *= factor; }
    Vec3d   vector() const { return this->b - this->a; }
    Vec3d   unit_vector() const { return (length() == 0.0) ? Vec3d::Zero() : vector().normalized(); }
    double  length() const { return vector().norm(); }

    Vec3d a;
    Vec3d b;

    static const constexpr int Dim = 3;
    using Scalar = Vec3d::Scalar;
};

BoundingBox get_extents(const Lines &lines);

} // namespace Slic3r

// start Boost
#include <boost/polygon/polygon.hpp>
namespace boost { namespace polygon {
    template <>
    struct geometry_concept<Slic3r::Line> { typedef segment_concept type; };

    template <>
    struct segment_traits<Slic3r::Line> {
        typedef coord_t coordinate_type;
        typedef Slic3r::Point point_type;

        static inline point_type get(const Slic3r::Line& line, direction_1d dir) {
            return dir.to_int() ? line.b : line.a;
        }
    };
} }
// end Boost

#endif // slic3r_Line_hpp_
