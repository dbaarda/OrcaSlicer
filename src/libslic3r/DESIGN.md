# OrcaSlicer libslic3r Design

::: contents
**Contents**
:::

## Introduction

This is attempts to collect history/details/ideas related to libslic3r in
OrcaSlicer. The libslic3r bits are used to implement most of the non-gui
related parts of OrcaSlicer.

## Design Details

### Inspirations and Origins

See https://github.com/SoftFever/OrcaSlicer?tab=readme-ov-file#some-background
for details of the origins/history.

Looking at the code it appears to have evolved fairly significanty over time,
and looks nearly nothing like the original Slic3r code now. It appears to have
pulled in features from many different places and incorporated them without
wasting much effort on trying to make them consistent with the rest of the
code.

This means the coding style and formatting varies considerably, with even
individual files having varying tab/space indenting, trailing whitespace, and
inconsistent formatting. This inconsistency means the code provides no
guidance or example of any sort of preferences, further encouraging ad-hock
approaches and formatting in contributions. There is a `.clang-format` config
in the repo, but AFAIKT it has not been tuned to represent any kind of
preferred style and has not been actively used. See #10917 for a proposed fix.

It has multiple different dependencies, often with redundant feature-set
overlaps, with different features using different dependencies to do similar
kinds of things. Some of those dependencies include flexible hooks that
integrate them into libslic3r, allowing features to be implemented or
integrated using the API/style/approach of that particular dependency instead
of any kind of "native" libslic3r API.

There is also significant redundancy in the code itself, with many features
re-implementing the same things repeatedly in slightly different ways with
varying degrees of optimization and code cleanliness. This is sometimes
because the core code is missing those frequently required features, often it
is hard to find them to reuse them, sometimes the interface is clumsy making
them hard to reuse, or the new code being integrated just happened to
implement them too.

Note these observations are not criticism. The current state of the code is a
direct result and mostly required outcome of rapidly integrating features from
many different sources, which is what has made OrcaSlicer so good.
Incorporating and supporting multiple different libraries with redundant
feature-overlaps and different APIs makes it much easier to incorporate
features that were implemented using those different libraries.

However, it's often a good idea to try and clean up some accumulated tech-debt
and rationalionalise things a bit. It would be good to;

1. Reformat all the core code using an automated code formatting tool and
require future submits comply/use it.

1. Identify all the different library dependencies and where/how they are
used. Where there are overlaps in functionality provide recommendations on
which to use for future development.

1. Clean up redundency in the core code, identifying and providing core
methods to peform commonly duplicated functionality.

### Design Philosophy

#### Geometry

##### Spacial data

When iterating through grids of spatial data the most common order is through
x (columns), then y (rows), then z (layers), with that being the order of the
inner-most to outer-most loops/methods/etc. This means for the best memory
locality data should be arranged in vectors or packed in memory as layers (z)
of rows (y) of columns (x), and indexed as `grid[z][y][z]`. This ensures
better cache hit-rates when iterating through the whole grid. It also makes it
easy to pick-off individual x-rows with `grid[z][y]`, or x-y-layers with
`grid[z]`, to pass as arguments to methods that handle individual rows or
layers.

So in general you should use the `[z][y][x]` layout unless there is a strong
reason not to. If there is a strong reason not to you should add comments
explaining why.

##### Distance methods

There are lots of places where you need to get the distance between different
kinds of geometry. For different kinds of geometry there are different ways to
measure distance, and sometimes even multiple different ways to measure distance
for the same geometry.

Often you want to find the minimum distance for multile different instances,
and there can be fast heuristics for comparing distances that don't always
require computing the distance accurately. Sometimes, but less often, you want
to find the maxiumum distance.

Sometimes you only need to find the minimum distance. Sometimes you need to
find the closest instance and don't need the distance. Sometimes you need
both. Sometimes you need to find the closest point, which might not be a point
defining the instance but a derived point on a line or plane of the instance.

At a fundimental low level distance is always about the length of a vector, so
it makes sense to start there with some consistent methods.

// get the result of min(v.squaredNorm(), len2) but short-circuit evaluating
//the length-squared if len2 is clearly shorter.
coord2_t min_length2(Vec2crd v, coord2_t len2)
{
   coord2_t x2,y2,l2;
   if ((v.x > len2) || (v.x < -len2) ||
       (v.y > len2) || (v.y < -len2) ||
       ((x2 = sqr(v.x)) > len2) ||
       ((y2 = sqr(v.y)) > len2) ||
       ((l2 = x2 + y2) >= len2))
     return len2;
   return l2;
}

// Note min_length2() is cheaper, but if you really need the actual length...
coord2_t min_length(Vec2crd v, coord_t len)
{
   coord2_t l2;
   if ((v.x > len) || (v.x < -len) ||
       (v.y > len) || (v.y < -len) ||
       // note sqr() is cheaper than sqrt() so this heuristic is worth it if it hits.
       ((l2 = sqr(x) + sqr(y)) >= sqr(len))
     return len;
   return sqrt(l2);
}

// returns true if length2(v) < len2, and updates len2 to the new minimum length. Use
// like;

// for (v : vects) if (is_length2_lt(v, minlen)) min_v = v;
bool is_length2_lt(Vec2crd v, coord2_t& len2)
{
  coord2_t l2 = len2;
  return ((len2 = min_length2(v, len2)) < l2)
}

// returns the shorter of v and vmin, updating len2 to the min length2.
Vect min_vect_len2(Vect v, Vect vmin, coord2_t& len2)
{
  return is_length2_lt(v, len2) ? v : vmin;
}

// returns true if v is shorter than len2, updating vmin and len2 to the
// shortest vector and len2.
bool is_vect_lt(Vect2crd v, Vect2crd& vmin, coord2_t& len2)
{
  if (v
}

###### Existing Implementations and where used.

Eigen::

* norm, squaredNorm: returns Scalar type so will overflow for coord, so they need
cast<double> or cast<int64_t> before using them. Used in lots of places.

* projected distance: infinite-line projected point and projected point
distance. infinite line is define with a point on the line and a direction
vector. Doesn't look particularly optimized, but maybe uses vector ops? Not
used?

Point: lots



##### Coordinates

There are several different coordinate types, with the default `coord_t` being
an alias for `int64_t` (but the code includes disabled options to use the
previous `int32_t`) used as a "fixed point" integer where i*10e-6 is the
distance in mm. There are scale() and unscale() methods for scaling to/from
coord_t, usually to/from double. It's also possible to cast to/from other
types without re-scaling.

The old int32_t range of coord_t was sufficient to fairly easily represent
points and distances within a typical 3D printer volume, and fits an interval
of +-2.147m. However, this means the square of a distance with an absolute
magniture of more than 0.065536mm would overflow (2^16 * 1e-6). Note that
calculating a distance involves calculating the sum-of-the-squares of the
axis-distances, so just calculating a distance would overflow. For this reason
distance calculations still include casting the coord_t to double (without
scaling) before calculating the distance as a double. However, some places in
the code cast to `int64_t` to calculate distances. Using int64_t is probably
more efficient unless you need the distance as a double anyway.

The new int64_t for coord_t fits an interval of +-9.2 million Km, but it is
important to note that when squaring distances anything longer than 4.294m
(2^32*1e-6) will still overflow. That's probably sufficient for even large
print volumes, but might be a concern if this code is ever used for things
like 3D printing full scale houses or boats. Note the old int32_t option
includes commented out options and notes about using 1e-5 as the scaling
factor for large printers, but IMHO 1e-3 would be better. That would mean
int64_t distances up to 4.2Km (2^32*1e-3) could be squared without
overflowing, while still giving micro-meter resolution.

We should probably switch to using int64_t most of the time, but probably need
to also support switching to double for the cases where it is usefull.

Although the code currently uses int64_t there are still some advantages to
using int32_t for applications that don't need the extra range. It would be
nice if this really was a clean build option with the code correctly compiling
with either, because using int32_t would be faster and use less memory.
However, it seems large-scale printers and large numbers of build plates
support was actually the motivating factor for the switch, so it looks like
int64_t is going to be the default. IMHO it should either be made to
compile/work for both as a build option, or all the int32_t assumptions and
support should be ripped out.

It would also be nice to be able to pick a different scaling factor for
applications that require much larger scales. This is actually much easier, as
the code for scaling/unscaling is well centralized and universally used.

To correctly support coord_t as either int32_t or int64_t would need some kind
of abstraction-wrapper for handling distances safely. This would need to
include a suitably sized/scaled type for holding squared-distances. Even if
only int64_t is used for coord_t having this abstraction could still be useful
future expansions.

Note that I also considered the idea of a "fast distance squared" that used a
rescaled int32_t for squared distances with;

```c++
using coord2_t = int32_t; //squared distance type.

coord2_t sqr(coord_t l)
{
  coord2_t l2 = l >> 16;
  return l2 * l2;
}

coord_t sqrt(coord2_t l2)
{
  coord_t l = std::sqrt(l2)
  return l << 16;
}
```

This would give approxmate distances that are within +-0.065536mm ((1<<16) *
1e-6) which is probably sufficient for many applications. However, performance
testing compared to just using an int64_t shows at best marginal wins, and
sometimes it's slower, so its not worth it. However this technique could maybe
be useful for coord_t = int64_t if there was a need for calculating distances
longer than 4.294m.

Note C typedefs are just type aliases and are not strongly typed, so you can't
e.g. overload sqrt for different typdefs of the same base type. So to do
strong-typing you need to hide it inside a struct. Or maybe you can do some
kind of template magic. Boost includes
[strong_typedef.hpp](https://www.boost.org/doc/libs/1_63_0/libs/serialization/doc/strong_typedef.html)
that can be used for this. Alternatively there is [fluentcpp on strong
types](https://www.fluentcpp.com/posts/#strong-types) [available on
github](https://github.com/joboccara/NamedType) which looks like it has
stronger and/or more configurable controls. In theory it would be possible to
define strong types and operators so that any time you multiply a coord_t by
another coord_t the result is a coord2_t, and it's not possible to
add/subtract/multiply/divide a coord_t and a coord2_t without explicit
casting.

It would also be possible to copy/extend Eigen's coordinate abstractions to add
support for a Derived::Squared type to indicate the type to use for operations
that square the values, and use that for all the distance related operations.
Note the following Eigen docs cover using a plugin to add methods and the
example adds length(), distanceTo() squaredDistanceTo() etc;

https://libeigen.gitlab.io/eigen/docs-nightly/TopicCustomizing_Plugins.html


##### Features

There are several feature types of increasing complexity

In general the primitive 2D features are;

* Point - a single 2D point. Implemented as an Eigen::Vector.

* Line - a pair of 2D Points defining a finite line.

* BoundingBox - a coordinate aligned rectangular area defined by a min and max
corner points

* MultiPoint - an (ordered?) collection of 2D points.

* Polyline - a MultiPoint that defines a 2D sequence of connected lines.

* Polygon - a Polyline that is closed by an implied line between the last and
  first points. To be valid the circular line cannot cross itself. They have
  an "inside" and "outside". If the points are counter-clockwise the "inside"
  is inside the circle of lines, and if it's counter-clockwise it is outside
  the circle of lines.

* ExPolygon - A 2D area with an external contour Polygon that can have multiple
  Polygon holes inside it. The external contour should be counter-clockwise
  and the holes should be clockwise, defining an area "inside" the contor with
  the holes cut out of it. To be valid the Polygons must be valid, the holes
  must be fully contained within the contour, and the holes cannot overlap.


 For a 3D polygon it has an
"inside" and "outside" face, with the "outside" face being the one where the
points are counter-clockwise when looking at it.


### General Architecture


## Indexes


### Definitions


### References

https://github.com/SoftFever/OrcaSlicer
