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
preferred style and has not been actively used.

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
require computing the distance accurately.

Sometimes you only need to find the minimum distance. Sometimes you need
to find the closest instance and don't need the distance. Sometimes you need
both. Sometimes you need to find the closest point, which might not be a point
defining the instance but a point on a line or plane that the instance defines.

##### Coordinates

There are several different coordinate types, with the default `coord_t` being
an alias for `int32_t` (unless a compiler option toggles it to `int64_t`) used
as a "fixed point" integer where i*10e-6 is the distance in mm. There are
scale() and unscale() methods for scaling to/from coord_t, usually to/from
double. It's also possible to cast to/from other types without re-scaling.

The range of coord_t is sufficient to fairly easily represent points and
distances within a typical 3D printer volume, and fits an interval of
(-2147.48mm, +2147.48mm). However, this means the square of a distance with an
absolute magniture of more than 46.34mm will overflow. Note that calculating a
distance involves calculating the sum-of-the-squares of the axis-distances, so
just calculating a normal distance will probably overflow. For this reason
distance calculations typically cast the coord_t to double (without scaling)
before calculating the distance as a double. However, some places in the code
cast to `int64_t` to calculate distances, and that is probably more efficient.

We should probably switch to using int64_t most of the time, but probably need
to also support switching to double for the cases where it is usefull.

It might make sense to copy/extend Eigen's coordinate abstractions to add
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
