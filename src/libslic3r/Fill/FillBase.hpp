#ifndef slic3r_FillBase_hpp_
#define slic3r_FillBase_hpp_

#include <assert.h>
#include <memory.h>
#include <float.h>
#include <stdint.h>
#include <stdexcept>

#include <type_traits>

#include "../libslic3r.h"
#include "../BoundingBox.hpp"
#include "../Exception.hpp"
#include "../Utils.hpp"
#include "../ExPolygon.hpp"
//BBS: necessary header for new function
#include "../PrintConfig.hpp"
#include "../Flow.hpp"
#include "../ExtrusionEntity.hpp"
#include "../ExtrusionEntityCollection.hpp"
#include "../ShortestPath.hpp"

namespace Slic3r {

class Surface;
enum InfillPattern : int;

namespace FillAdaptive {
    struct Octree;
};

// Infill shall never fail, therefore the error is classified as RuntimeError, not SlicingError.
class InfillFailedException : public Slic3r::RuntimeError {
public:
    InfillFailedException() : Slic3r::RuntimeError("Infill failed") {}
};

struct LockRegionParam
{
    LockRegionParam() {}
    std::map<float, ExPolygons> skin_density_params;
    std::map<float, ExPolygons> skeleton_density_params;
    std::map<Flow, ExPolygons>  skin_flow_params;
    std::map<Flow, ExPolygons>  skeleton_flow_params;
};

struct FillParams
{
    // Is this a solid fill?
    inline bool  is_solid() const { return is_approx_ge(density, 1.0f); }
    // Should the fill lines have a connection line around the inner perimeter?
    inline bool  dont_connect() const { return anchor_length_max < 0.05f; }

    // Get the line_spacing from flow_spacing and multiline.
    inline coordf_t get_line_spacing() const {
        assert(flow.spacing() && "params.flow must be set first.");
        return flow.spacing() * multiline;
    }
    // Set multiline required to give as close as possible to the specified line_spacing.
    // This requires that flow be set correctly first.
    inline void set_line_spacing(const coordf_t line_spacing) {
        assert(flow.spacing() && "params.flow must be set first.");
        assert(is_approx_le(flow.spacing(), line_spacing) && "line_spacing cannot be less than flow_spacing.");
        multiline = int(std::round(line_spacing / flow.spacing()));
    }
    // Get the fill_spacing from line_spacing and density.
    inline coordf_t get_fill_spacing() const {
        assert(flow.spacing() && "params.flow must be set first.");
        return get_line_spacing() / density;
    }
    // Set the density required to give the specified fill_spacing.
    // This requires that flow and multiline be set correctly first.
    inline void set_fill_spacing(const coordf_t fill_spacing) {
        assert(flow.spacing() && "params.flow must be set first.");
        assert((fill_spacing >= get_line_spacing()) && "fill_spacing cannot be less than line_spacing.");
        density = double(get_line_spacing() / fill_spacing);
    }

    ExtrusionRole   extrusion_role{ ExtrusionRole(0) };
    InfillPattern pattern{ ipRectilinear };

    // Flow to use.
    Flow        flow;

    // Number of flow-lines per fill-line.
    int         multiline       { 1 };

    // Fill density (ratio of line-spacing to fill-spacing in the range <0, 1>).
    float       density         { 0.f };

    // Unscaled layer height (which can be different to flow.height()).
    coordf_t    layer_height    { 0.f };

    // Unscaled length of an infill anchor along the perimeter.
    // 1000mm is roughly the maximum length line that fits into a 32bit coord_t.
    float       anchor_length       { 1000.f };
    float       anchor_length_max   { 1000.f };

    // Unscaled G-code resolution.
    double      resolution          { 0.0125 };

    bool        monotonic{ false }; // Use strictly left to right for better surface quality.
    bool        use_arachne{ false }; // Use arachne variable line widths.
    bool        dont_adjust{ true }; // Don't adjust flow and spacing to fill the space evenly.
    bool        dont_sort{ false }; // Don't sort the lines, just simply connect them.
    bool        can_reverse{true}; // The lines can be reversed.

    // For Honeycomb.
    // we were requested to complete each loop;
    // in this case we don't try to make more continuous paths
    // TODO(dbaarda): remove this if obsolete.
    //bool        complete      { false };

    // For Lateral lattice
    coordf_t    lateral_lattice_angle_1    { 0.f };
    coordf_t    lateral_lattice_angle_2    { 0.f };

    // For Lateral Honeycomb
    float       infill_overhang_angle    { 60 };

    // For Gyroid: when true, use the parameterized "optimized" variant.
    bool        gyroid_optimized { false };

    // BBS
    const           PrintRegionConfig* config{ nullptr };

    // Orca: forced print order of surface fill loops/fragments for center-based patterns
    // (Concentric, Archimedean Chords, Octagram Spiral). Default keeps shortest-path ordering.
    SurfaceFillOrder fill_order { SurfaceFillOrder::Default };

    float           horiz_move{0.0}; //move infill to get cross zag pattern
    bool            symmetric_infill_y_axis{false};
    coord_t         symmetric_y_axis{0};
    bool            locked_zag{false};
    float           infill_lock_depth{0.0};
    float           skin_infill_depth{0.0};
    CenterOfSurfacePattern center_of_surface_pattern{CenterOfSurfacePattern::Each_Surface};
};
static_assert(IsTriviallyCopyable<FillParams>::value, "FillParams class is not POD (and it should be - see constructor).");

class Fill
{
public:
    // Index of the layer.
    size_t      layer_id;
    // Z coordinate of the top print surface, in unscaled coordinates
    coordf_t    z;

    // Adjusted flow to use. From params.flow possibly adjusted by scaling it to fit the surface better.
    inline Flow        flow() const { return _flow; };
    // Unscaled adjusted flow-spacing (space for a single line). From the adjusted flow.spacing().
    inline coordf_t    flow_spacing() const { return _flow.spacing(); }
    // Unscaled adjusted line-spacing (space for a single multi-line). The adjusted flow-spacing with multiline applied.
    inline coordf_t    line_spacing() const { return _line_spacing; }
    // Unscaled adjusted fill-spacing (spacing between adjacent multilines). The adjusted line-spacing with density applied.
    inline coordf_t    fill_spacing() const { return _fill_spacing; }

    // Scaled and adjusted flow-spacing (space of a single line).
    coord_t     scaled_flow_spacing() const { return _flow.scaled_spacing(); }
    // Scaled and adjusted line-spacing (space of a single multiline),
    coord_t     scaled_line_spacing() const { return scaled(_line_spacing); }
    // Scaled and adjusted fill-spacing (spacing between adjacent multilines).
    coord_t     scaled_fill_spacing() const { return scaled(_fill_spacing); }

    // Unscaled infill / perimeter overlap.
    coordf_t    overlap;
    // in radians, ccw, 0 = East
    float       angle;

    // Orca: Fill direction is fixed absolute angle if SurfaceFillParams.fixed_angle or config.ironing_angle_fixed
    bool        fixed_angle{false};

    // In scaled coordinates. Used by the concentric infill pattern to clip the loops to create extrusion paths.
    coord_t     loop_clipping;
    // In scaled coordinates. Bounding box of the 2D projection of the object.
    BoundingBox bounding_box;

    // Octree builds on mesh for usage in the adaptive cubic infill
    FillAdaptive::Octree* adapt_fill_octree = nullptr;

    // PrintConfig and PrintObjectConfig are used by infills that use Arachne (Concentric and FillEnsuring).
    // Orca: also used by gap fill function.
    const PrintConfig       *print_config        = nullptr;
    const PrintObjectConfig *print_object_config = nullptr;

    // BBS: all no overlap expolygons in same layer
    ExPolygons  no_overlap_expolygons;
    bool dont_alternate_fill_direction = false;

    static float infill_anchor;
    static float infill_anchor_max;

public:
    virtual ~Fill() {}
    virtual Fill* clone() const = 0;

    static Fill* new_from_type(const InfillPattern type);
    static Fill* new_from_type(const std::string &type);
    static bool  use_bridge_flow(const InfillPattern type);

    void         set_bounding_box(const Slic3r::BoundingBox &bbox) { bounding_box = bbox; }
    BoundingBox  extended_object_bounding_box() const;
    // Use bridge flow for the fill?
    virtual bool use_bridge_flow() const { return false; }

    // Do not sort the fill lines to optimize the print head path?
    virtual bool no_sort() const { return false; }

    virtual bool is_self_crossing() = 0;

    // Return true if infill has a consistent pattern between layers.
    virtual bool has_consistent_pattern() const { return false; }

    // Perform the fill.
    virtual Polylines fill_surface(const Surface *surface, const FillParams &params);
    virtual ThickPolylines fill_surface_arachne(const Surface* surface, const FillParams& params);
    virtual void set_lock_region_param(const LockRegionParam &lock_param){};
    // BBS: this method is used to fill the ExtrusionEntityCollection.
    // It call fill_surface by default
    virtual void fill_surface_extrusion(const Surface *surface, const FillParams &params, ExtrusionEntitiesPtr &out);

protected:
    Fill() :
        layer_id(size_t(-1)),
        z(0.),
        // Infill / perimeter overlap.
        overlap(0.),
        // Initial angle is undefined.
        angle(FLT_MAX),
        loop_clipping(0),
        // The initial bounding box is empty, therefore undefined.
        bounding_box(Point(0, 0), Point(-1, -1))
        {}

    // The expolygon may be modified by the method to avoid a copy.
    virtual void    _fill_surface_single(
        const FillParams                & /* params */,
        unsigned int                      /* thickness_layers */,
        const std::pair<float, Point>   & /* direction */,
        ExPolygon                         /* expolygon */,
        Polylines                       & /* polylines_out */) {}

    // Used for concentric infill to generate ThickPolylines using Arachne.
    virtual void _fill_surface_single(const FillParams& params,
        unsigned int                   thickness_layers,
        const std::pair<float, Point>& direction,
        ExPolygon                      expolygon,
        ThickPolylines& thick_polylines_out) {}

    // Set `_flow`, `_line_spacing`, and `_fill_spacing` from params, optionally adjusted to fit `distance` if provided.
    // Note this will honour `params.dont_adjust` and will not adjust them if it is set.
    void _set_flow_and_spacing(const FillParams& params, const coordf_t distance=0.);
    // Set `_flow`, `_line_spacing`, and `_fill_spacing`, optionally adjusted to fit `distance` if provided.
    void _set_flow_and_spacing(const Flow& flow, const int multiline=1, const double density=1.0, const coordf_t distance=0.);
    // Adjust `_flow`, `_line_spacing`, and `_fill_spacing` to fit `distance`. Note this will scale them at most by roughly +-10%, or
    // for bridge flows will downscale by up to -20% to ensure flow diameter <= nozzle_diameter.
    void _adjust_flow_and_spacing(const coordf_t distance);
    // Scale `_flow`, `_line_spacing`, and `_fill_spacing` by `scale`.
    void _scale_flow_and_spacing(const double scale);
    // Set `_fill_spacing` directly.
    inline void _set_fill_spacing(const coordf_t fill_spacing)
    {
        assert(fill_spacing >= _line_spacing && "Cannot set fill_spacing < line_spacing.");
        _fill_spacing = fill_spacing;
    }
    // Set `_fill_spacing = _line_spacing / fill_density`.
    inline void _set_fill_density(const double fill_density)
    {
        assert(is_approx_le(fill_density, 1.0) && "Cannot set density>1.0.");
        _fill_spacing = _line_spacing / fill_density;
    }
    // Get the fill_density from `_line_spacing` and `_fill_spacing`.
    inline double _get_fill_density() const
    {
        return _line_spacing / _fill_spacing;
    }

    virtual float _layer_angle(size_t idx) const { return fixed_angle ? 0.f : (idx & 1) ? float(M_PI/2.) : 0.f; }

    virtual std::pair<float, Point> _infill_direction(const Surface *surface) const;

    // Orca: Dedicated function to calculate gap fill lines for the provided surface, according to the print object parameters
    // and append them to the out ExtrusionEntityCollection.
    void _create_gap_fill(const Surface* surface, const FillParams& params, ExtrusionEntityCollection* out);

private:
    // Calculate the scaling adjustment required for `spacing` to make multiple spaces fit within `distance`.
    // It returns the scale within `min_scale` to `max_scale` closest to 1.0 that leaves the smallest gap. If the distance is too small to fit even a single
    // space at the `min_scale`, it returns 1.0 for "don't bother to scale it".
    static double _adjust_spacing_scale(const coordf_t distance, const coordf_t spacing, const double min_scale=1.0f, const double max_scale=1.2f);

    // Adjusted flow to use. From `params.flow` possibly adjusted to fit a distance.
    Flow        _flow;
    // Unscaled adjusted line-spacing (space for a single multi-line). The adjusted flow-spacing with multiline applied.
    coordf_t    _line_spacing{0.};
    // Unscaled adjusted fill-spacing (spacing between adjacent multilines). The adjusted line-spacing with density applied.
    coordf_t    _fill_spacing{0.};

public:
    static void connect_infill(Polylines &&infill_ordered, const ExPolygon &boundary, Polylines &polylines_out, const FillParams &params);
    static void connect_infill(Polylines &&infill_ordered, const Polygons &boundary, const BoundingBox& bbox, Polylines &polylines_out, const FillParams &params);
    static void connect_infill(Polylines &&infill_ordered, const ConstPolygonPtrs &boundary, const BoundingBox &bbox, Polylines &polylines_out, const FillParams &params);

    static void chain_or_connect_infill(Polylines &&infill_ordered, const ExPolygon &boundary, Polylines &polylines_out, const FillParams &params);

    static void connect_base_support(Polylines &&infill_ordered, const ConstPolygonPtrs &boundary_src, const BoundingBox &bbox, Polylines &polylines_out, const FillParams &params);
    static void connect_base_support(Polylines &&infill_ordered, const Polygons &boundary_src, const BoundingBox &bbox, Polylines &polylines_out, const FillParams &params);

};
   //Fill  Multiline
   void multiline_fill(Polylines& polylines, const int n_lines, coordf_t flow_spacing);
} // namespace Slic3r

#endif // slic3r_FillBase_hpp_
