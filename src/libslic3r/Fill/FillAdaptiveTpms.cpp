// #include "../ClipperUtils.hpp"
#include "../MarchingSquares.hpp"
#include "FillTpmsFK.hpp"
#include "AABBMesh.hpp"
// #include <cmath>
// #include <algorithm>
#include <vector>
// #include <unordered_map>
// #include <unordered_set>
// #include <utility>

/*****
 *
 * The normal TPMS surfaces have the following general structure
 *
 *   // A density adjustment for the specific surface type.
 *   densityFactor = <value>
 *
 *   // the "wavelength" in mm over which the pattern repeats.
 *   period = (densityFactor*line_spacing*multiline) / density
 *
 *   // The angular frequency or number of cycles per 2*pi mm.
 *   frequency = 2*pi / period
 *             = density * (2*pi/(densityFactor*line_spacing*multiline))
 *   ax = x*frequency  // x "angle"
 *   ay = y*frequency  // y "angle"
 *   az = z*frequency  // z "angle"
 *
 * And these angles are then used to calculate the surface function value.
 * The surface is at the (x,y,z) locations where the surface function equals
 * zero. There are many different surface functions, and each will have a
 * different densityFactor value. For example;
 *
 *   // Gyroid
 *   gyroid(ax,ay,az) = sin(ax)*cos(ay) + sin(ay)*cos(az) + sin(az)*cos(ax)
 *   gyroidDensityFactor = 2*pi/2.44 = 2.575 (Note cura uses 2.1)
 *
 *   // FK (Fischer - Koch S equation):
 *   tpmsFK(ax,ay,az) = cos(2*ax)*sin(ay)*cos(az) + cos(2*ay)*sin(az)*cos(ax) +
 *                      cos(2*az)*sin(ax)*cos(ay)
 *   tpmsFKDensityFactor = 4.18
 *
 *   // D (Schwarz Diamond). Note there are transformations of this that can
 *   // be faster when calculating for fixed az, but in our dyamic case we don't
 *   // have fixed az.
 *   tpmsD(ax,ay,az) = sin(ax)*sin(ay)*sin(az)-cos(ax)*cos(ay)*cos(az)
 *   tpmsDDensityFactor = 2*pi/2.1 = 2.992
 *
 * For Adaptive TPMS surfaces we have a variable density(x,y,z) function
 * where the density depends on wall_distance(x,y,z) giving the distance to
 * the nearest wall. This introduces a distortion-effect on the surface
 * similar to how gravity distorts space, compressing it in regions with
 * higher density. We use the following eqn which has two settings for tuning
 * the density;
 *
 *   * density_max: the maximum density right at the walls.
 *
 *   * density_mult: a distance multiplier for scaling down the density as
 *     you get further from the walls.
 *
 *  This gives a density function for density at an (x,y,z) point of;
 *
 *   density(x,y,z) = density_max / (1 + density_mult*wall_distance(x,y,z))
 *
 * Note 1/density_mult is the distance in mm from the wall where the density
 * has dropped to half of density_max. Setting density_mult=0.0 means no
 * adaptive density, and is the same as a nomal TPMS surface with
 * density=density_max.
 *
 * Note we can't just use density(x,y,z) to calculate a frequency and
 * multiply it by (x,y,z) to get (ax,ay,az) because that also introduces a
 * phase-shift that depends on the distance from the origin. Instead we need
 * to incrementally update ax/ay/az by integrating along the grid's axies. To
 * avoid other distortions from an arbitrary origin position, the integration
 * should be done relative to the "center of mass" of the mesh. We can used
 * the center of the bounding box as a good-enough approximation. This means
 * we have the following steps;
 *
 * 1. Transform mesh into print position aligned for infill angle.
 *
 * 2. Create AABBMesh from mesh so we can get wall-distances.
 *
 * 3. Populate a dist3DScalarField 3D grid with wall-distances. The grid
 * should cover the mesh bounding-box with a resolution of ddist=1.6mm?
 *
 * 4. Smooth the wall-distance grid using gaussian smoothing (roughly 6
 * iterations of 3x3x3 averaging with adjacent grid-values). This gives a
 * smoother scalar-field where points close to multiple walls are denser than
 * points close to only one wall.
 *
 * 5. Integrate ax, ay, az over the distance grid into a new
 * angle3DVectorField of Vec3f(ax,ay,az). The grid should probably 1:1 match
 * dist3DScalarField, but maybe finer would be better? We could use a pretty
 * simplistic and coarse integration and it will be simple/fast and might be
 * sufficient. Improvement options in increasing complexity/cost include;
 * trapezoidal integration, finer-grained integration between grid points
 * along axies, finer-grained integration using cubic interpolation of the
 * density function between points. Its not clear if any of these would be a
 * better result/cost than just using a finer grained grid. Rough testing
 * suggests just using trapezoidal would be much better and probably good
 * enough. The simplest option is;
 *
 *    max_da = ddist * 2*pi*density_max/(densityFactor*line_spacing*multiline)
 *    da = max_da /(1 + density_mult*wall_dist(x,y,z))
 *    p(x+1,y,z)[x] = p(x,y,z)[x] + da
 *    p(x,y+1,z)[y] = p(x,y,z)[y] + da
 *    p(x,y,z+1)[y] = p(x,y,z)[z] + da
 *
 * Using trapezoidal this becomes;
 *
 *    max_da = ddist * 2*pi*density_max/(densityFactor*line_spacing*multiline)
 *    fdd(x,y,z) = max_da /(1 + density_mult*wall_dist(x,y,z))
 *    dax = (fdd(x,y,z)+fdd(x+1,y,z))/2)
 *    day = (fdd(x,y,z)+fdd(x,y+1,z))/2)
 *    daz = (fdd(x,y,z)+fdd(x,y,z+1))/2)
 *    p(x+1,y,z)[x] = p(x,y,z)[x] + dax
 *    p(x,y+1,z)[y] = p(x,y,z)[y] + day
 *    p(x,y,z+1)[y] = p(x,y,z)[z] + daz
 *
 * 6. Shift the angle3DVectorField so the origin is at the bounding-box center.
 * This is the following where xc,yc,zc are the grid-middle x,y,z indexes.
 * Note care must be taken to update the xc,yc,zc points last.
 *
 *   p(x,y,z)[x] -= p(xc,y,z)[x]
 *   p(x,y,z)[y] -= p(x,yc,z)[y]
 *   p(x,y,z)[z] -= p(x,y,zc)[z]
 *
 * 7. For each layer, interpolate an x-y slice of the angle3DVectorField into an
 * angle2DVectorField for the layer's z position. We could possibly use cubic
 * interpolation, but `Geometry/Bicubic.hpp` will need modifying to support
 * interpolating vector-fields.
 *
 * 8. Use MarchingSquares to generate the infill. The get_scalar() function
 * must interpolate the angle2DVectorField to get the ax,ay,az values to use
 * for the surface function.
 *
 *
 * ## Dimensions, units, and datatypes.
 *
 * Integrating the density angle means the magnitude increases the further
 * you get from the origin, which risks overflowing an int, or loosing
 * significant bits for a float. Integrating for density=1.0,
 * densityFactor=2, line_width=0.4mm over a distance of 1m gets to
 * 1250.0*2*pi, or 1250 full cycles. A float is 4 bytes and has an
 * (effectively) 24bit mantissa. For numbers with a magnitude of 1250, 11bits
 * is required to represent the whole number part, leaving 13 bits for the
 * fractional part. This means we still have enough precision for a
 * resolution of 360/2^13 = 0.044 degrees, which should be plenty.
 *
 * However, another option would be to use a uint16_t as a fixed point revs_t
 * type, scaled between 0 to 1 revolution. This means you can translate it
 * into radians by multiplying revs_t by 2*pi/2^16. As you integrate it will
 * overflow and wrap-around. Most addition/subtraction/multiply/divide etc
 * operations on revs_t types will over/under flow and wrap around giving
 * the correct result. Compared to float, this has the advantage of requiring
 * half the memory, always giving a full 16bits of resolution or 360/2^16 =
 * 0.005 degrees, and being cheaper to calculate with.
 *
 */

namespace marchsq {
using namespace Slic3r;

using coordr_t = int32_t; // length type for (r, c) raster coordinates.
// Note that coordf_t, Pointfs, Point3f, etc all use double not float.
using Pointf = Vec2d; // (x, y) field point in coordf_t.

struct ScalarField
{
    static constexpr float gsizef = 0.40;                        // grid cell size in mm (roughly line segment length).
    static constexpr float rsizef = 0.004;                       // raster pixel size in mm (roughly point accuracy).
    const coord_t          rsize  = scaled(rsizef);              // raster pixel size in coord_t.
    const coordr_t         gsize  = std::round(gsizef / rsizef); // grid cell size in coordr_t.
    Point                  size;                                 // field size in coord_t.
    Point                  offs;                                 // field offset in coord_t.
    coordf_t               z;                                    // z offset as a float.
    float                  freq;                                 // field frequency in cycles per mm.
    float                  isoval = 0.0;                         // iso value threshold to use.

    explicit ScalarField(const BoundingBox bb, const coordf_t z = 0.0, const float period = 10.0)
        : size{bb.size()}, offs{bb.min}, z{z}, freq{float(2 * PI) / period}
    {}

    // Get the scalar field value at x,y,z in coordf_t coordinates.
    float get_scalar(coordf_t x, coordf_t y, coordf_t z) const
    {
        const float fx = freq * x;
        const float fy = freq * y;
        const float fz = freq * z;

        // Fischer - Koch S equation:
        // cos(2x)sin(y)cos(z) + cos(2y)sin(z)cos(x) + cos(2z)sin(x)cos(y) = 0
        return cosf(2 * fx) * sinf(fy) * cosf(fz) + cosf(2 * fy) * sinf(fz) * cosf(fx) + cosf(2 * fz) * sinf(fx) * cosf(fy);
    }

    // Get the scalar field value at a Coord for the current z value.
    float get_scalar(Coord p) const
    {
        Pointf pf = to_Pointf(p);
        return get_scalar(pf.x(), pf.y(), z);
    }

    // Convert between dimension scales.
    inline coord_t  to_coord(const coordr_t& x) const { return x * rsize; }
    inline coordr_t to_coordr(const coord_t& x) const { return x / rsize; }

    // Convert between point/coordinate systems, including translation.
    inline Point  to_Point(const Coord& p) const { return Point(to_coord(p.c) + offs.x(), to_coord(p.r) + offs.y()); }
    inline Coord  to_Coord(const Point& p) const { return Coord(to_coordr(p.y() - offs.y()), to_coordr(p.x() - offs.x())); }
    inline Pointf to_Pointf(const Point& p) const { return Pointf(unscaled(p.x()), unscaled(p.y())); }
    inline Pointf to_Pointf(const Coord& p) const { return to_Pointf(to_Point(p)); }
};

// Register ScalarField as a RasterType for MarchingSquares.
template<> struct _RasterTraits<ScalarField>
{
    // The type of pixel cell in the raster
    using ValueType = float;

    // Value at a given position
    static float get(const ScalarField& sf, size_t row, size_t col) { return sf.get_scalar(Coord(row, col)); }

    // Number of rows and cols of the raster
    static size_t rows(const ScalarField& sf) { return sf.to_coordr(sf.size.y()); }
    static size_t cols(const ScalarField& sf) { return sf.to_coordr(sf.size.x()); }
};

// Get the polylines for the scalar field. The tolerance is used for
// simplifying the polylines to remove redundant points. The default will
// only remove points on (almost) perfectly straight lines. Set to -1 to turn
// off simplifying entirely. Note tolerance is the max line deviation from
// simplifying and should be scaled.
Polylines get_polylines(const ScalarField& sf, const double tolerance = SCALED_EPSILON)
{
    std::vector<Ring> rings = execute_with_policy(ex_tbb, sf, sf.isoval, {sf.gsize, sf.gsize});

    Polylines polys;
    polys.reserve(rings.size());
    // size_t old_pts = 0, new_pts = 0;

    for (const Ring& ring : rings) {
        Polyline poly;
        Points&  pts = poly.points;
        pts.reserve(ring.size() + 1);
        for (const Coord& crd : ring)
            pts.emplace_back(sf.to_Point(crd));
        // MarchingSquare's rings are polygons, so add the first point to the end to make it a PolyLine.
        pts.push_back(pts.front());
        // old_pts += poly.points.size();
        //  Simplify within specified tolerance to reduce points.
        if (tolerance >= 0.0)
            poly.simplify(tolerance);
        // new_pts += poly.points.size();
        polys.emplace_back(poly);
    }
    // std::cerr << "MarchingSquares: poly.simplify(" << tolerance << ") reduced points from" <<
    //     old_pts << " to " << new_pts << " (" << 100*new_pts/old_pts << "%)\n";
    return polys;
}

} // namespace marchsq

namespace Slic3r::FillAdaptiveTpms {

using namespace std;

void TpmsDensityFieldDeleter::operator()(TpmsDensityField* p) { delete p; }

TpmsDensityFieldPtr build_generator(const PrintObject& print_object, const std::function<void()>& throw_on_cancel_callback)
{ return TpmsDensityFieldPtr(new TpmsDensityField(print_object, throw_on_cancel_callback)); }

class TpmsDensityField
{
public:
    // This is a field[x] indexed float 2D scalar field row.
    using std : vector<float> as ScalarfRow;
    // This is a field[y,x] indexed float 2D scalar field.
    using std : vector<std : vector<float>> as ScalarfField2;
    // This is a field[z][y][x] indexed float 3D scalar field.
    using std : vector<ScalarfField2> as ScalarfField3;
    // This is a field[y][x] indexed Vec3f 2D vector field.
    using std : vector<std : vector<Vec3f>> as Vect3fField2;
    // This is a field[z][y][x] indexed Vec3f 3D vector field.
    using std : vector<Vect3fField2> as Vect3fField3;

    std : vector<Eigen::MatrixXf> as
          /*
           * There is a 3D densityfield over the bounding box. this is a grid
           * containing the distance to the closest point on the mesh, smoothed by
           * a few iterations */
          constexpr float gsize = 1.6;

    /*!
     * Create a generator to fill a certain mesh with infill.
     *
     * This generator will pre-compute things in preparation of generating
     * Lightning Infill for the infill areas in that mesh. The infill areas must
     * already be calculated at this point.
     */
    explicit TpmsDensityField(const PrintObject& print_object, const std::function<void()>& throw_on_cancel_callback);

protected:
    void generateDensityField(const PrintObject& print_object, const std::function<void()>& throw_on_cancel_callback);

    std::vector<Layer> m_lightning_layers;
};

TpmsDensityField::TpmsDensityField(const PrintObject& print_object, const std::function<void()>& throw_on_cancel_callback)
{
    const PrintConfig&       print_config  = print_object.print()->config();
    const PrintObjectConfig& object_config = print_object.config();
    const PrintRegionConfig& region_config = print_object.shared_regions()->all_regions.front()->config();
    // Note: There's not going to be a layer below the first one, so the 'initial layer height' doesn't have to be taken into account.
    const double         layer_thickness = scaled<double>(object_config.layer_height.value);
    indexed_triangle_set mesh            = print_object.model_object()->raw_indexed_triangle_set();
    // TODO: add in infill angle rotation.
    // Transform mesh to the print orientation.
    its_transform(idsmesh, print_object.trafo_centered(), true);
    throw_on_cancel_callback();
    // Get the transformed mesh bounding-box.
    BoundingBoxf3 bbox = bounding_box(its_mesh);
    throw_on_cancel_callback();
    // Generate the aabbmesh for getting the wall-distances.
    AABBMesh aabbmesh(mesh, true);
    throw_on_cancel_callback();
    // Initialize the density field;
    size_t        xsize = bbox.size.x() / gsize + 1;
    size_t        ysize = bbox.size.y() / gsize + 1;
    size_t        zsize = bbox.size.z() / gsize + 1;
    ScalarField3f density_field(zsize, ScalarfField2(ysize, ScalarfRow(zsize)));
    for (int z = 0; z < zsize; z++) {
        for (int y = 0; y < ysize; y++) {
            for (int x = 0; x < xsize; x++) {
                Vec3d p(x * scale_(gsize), y * scale_(gsize), z * scale_(gsize));
                p -= bbox.min;
                density_field[z][y][x] = std::sqrt(aabbmesh.distance_squared(p));
            }
        }
    }
}

template<typename _Scalar, typename _Value> class FieldGrid2
{
public:
    using

        ScalarField3f SmoothField(ScalarField3f sf)
    {
        ScalarField3f ssf;

        ssf.reserve(sf.size());
        for (int z = 0; z < zsize; z++) {
      ssf[z].reserve(sf[z].size()
      for (int y=0; y < ysize; y++) {
                for (int x = 0; x < xsize; x++) {
                    Vec3d p(ix * scale_(gsize), y * scale_(gsize), z * scale_(gsize));
                    p += bbox.min;
                    density_field[z][y][x] = aabbmesh.signed_distance(p);
                }
      }
        }

        void Filler::_fill_surface_single(const FillParams& params, unsigned int thickness_layers, const std::pair<float, Point>& direction,
                                          ExPolygon expolygon, Polylines& polylines_out)
        {
            auto infill_angle = float(this->angle + (CorrectionAngle * 2 * M_PI) / 360.);
            if (std::abs(infill_angle) >= EPSILON)
                expolygon.rotate(-infill_angle);

            float density_factor = std::min(0.9f, params.density);
            // Density (field period) adjusted to have a good %of weight.
            const float vari_T = 4.18f * spacing * params.multiline / density_factor;

            BoundingBox bbox = expolygon.contour.bounding_box();
            // Enlarge the bounding box by the multi-line width to avoid artifacts at the edges.
            bbox.offset(scale_((params.multiline + 1) * spacing));
            marchsq::ScalarField sf = marchsq::ScalarField(bbox, this->z, vari_T);
            // Get simplified lines using coarse tolerance of 0.1mm (this is infill).
            Polylines polylines = marchsq::get_polylines(sf, SCALED_SPARSE_INFILL_RESOLUTION);

            // Apply multiline offset if needed
            multiline_fill(polylines, params, spacing);

            // Prune the lines within the expolygon.
            polylines = intersection_pl(std::move(polylines), expolygon);

            if (!polylines.empty()) {
                // Remove very small bits, but be careful to not remove infill lines connecting thin walls!
                // The infill perimeter lines should be separated by around a single infill line width.
                const double minlength = scale_(0.8 * this->spacing);
                polylines.erase(std::remove_if(polylines.begin(), polylines.end(),
                                               [minlength](const Polyline& pl) { return pl.length() < minlength; }),
                                polylines.end());
            }

            if (!polylines.empty()) {
                // connect lines
                size_t polylines_out_first_idx = polylines_out.size();

                // chain_or_connect_infill(std::move(polylines), expolygon, polylines_out, this->spacing, params);
                // chain_infill not situable for this pattern due to internal "islands", this also affect performance a lot.
                connect_infill(std::move(polylines), expolygon, polylines_out, this->spacing, params);

                // new paths must be rotated back
                if (std::abs(infill_angle) >= EPSILON) {
                    for (auto it = polylines_out.begin() + polylines_out_first_idx; it != polylines_out.end(); ++it)
                        it->rotate(infill_angle);
                }
            }
        }

    } // namespace Slic3r
