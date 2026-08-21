#include "../ClipperUtils.hpp"
#include "../ExPolygon.hpp"
#include "../Surface.hpp"
#include "../VariableWidth.hpp"
#include "Arachne/WallToolPaths.hpp"

#include "FillConcentric.hpp"
#include <libslic3r/ShortestPath.hpp>

namespace Slic3r {

void FillConcentric::_fill_surface_single(
    const FillParams                &params, 
    unsigned int                     thickness_layers,
    const std::pair<float, Point>   &direction, 
    ExPolygon                        expolygon,
    Polylines                       &polylines_out)
{
    // no rotation is supported for this infill pattern
    BoundingBox bounding_box = expolygon.contour.bounding_box();
    
    // Adjust flow and spacing for the box size.
    this->_set_flow_and_spacing(params, unscaled(bounding_box.size().x()));
    coord_t scaled_flow_spacing = this->scaled_flow_spacing(); // scaled adjusted flow spacing.
    coord_t scaled_line_spacing = this->scaled_line_spacing(); // scaled adjusted mult-line spacing.
    coord_t scaled_fill_spacing = this->scaled_fill_spacing(); // scaled adjusted fill spacing.

    // Contract surface polygon by half line width to avoid excesive overlap with perimeter.
    // Note polygons are pre-contracted by params.flow.spacing()/2 which we account for.
    ExPolygons contracted = offset_ex(expolygon, float((scaled_line_spacing - params.flow.scaled_spacing())/2));

    Polygons loops = to_polygons(contracted);

    ExPolygons last { std::move(contracted) };
    while (! last.empty()) {
        last = offset2_ex(last, -(scaled_fill_spacing + scaled_line_spacing/2), +scaled_line_spacing/2);
        append(loops, to_polygons(last));
    }

    // generate paths from the outermost to the innermost, to avoid
    // adhesion problems of the first central tiny loops
    loops = union_pt_chained_outside_in(loops);

    // Orca: an outward fill order prints the innermost loops first instead.
    if (params.fill_order == SurfaceFillOrder::Outward)
        std::reverse(loops.begin(), loops.end());
    
    // split paths using a nearest neighbor search
    size_t iPathFirst = polylines_out.size();
    Point last_pos(0, 0);
    for (const Polygon &loop : loops) {
        polylines_out.emplace_back(loop.split_at_index(last_pos.nearest_point_index(loop.points)));
        last_pos = polylines_out.back().last_point();
    }

    // Apply multiline offset if needed
    multiline_fill(polylines_out, params.multiline, scaled_flow_spacing);

    // clip the paths to prevent the extruder from getting exactly on the first point of the loop
    // Keep valid paths only.
    size_t j = iPathFirst;
    for (size_t i = iPathFirst; i < polylines_out.size(); ++ i) {
        polylines_out[i].clip_end(this->loop_clipping);
        if (polylines_out[i].is_valid()) {
            if (j < i)
                polylines_out[j] = std::move(polylines_out[i]);
            ++ j;
        }
    }
    if (j < polylines_out.size())
        polylines_out.erase(polylines_out.begin() + j, polylines_out.end());
    //TODO: return ExtrusionLoop objects to get better chained paths,
    // otherwise the outermost loop starts at the closest point to (0, 0).
    // We want the loops to be split inside the G-code generator to get optimum path planning.
}

void FillConcentric::_fill_surface_single(const FillParams& params,
    unsigned int                   thickness_layers,
    const std::pair<float, Point>& direction,
    ExPolygon                      expolygon,
    ThickPolylines& thick_polylines_out)
{
    assert(params.use_arachne);
    assert(this->print_config != nullptr && this->print_object_config != nullptr);

    // no rotation is supported for this infill pattern
    Point   bbox_size = expolygon.contour.bounding_box().size();
    coord_t scaled_flow_spacing = this->scaled_flow_spacing();

    if (params.is_solid() && !params.dont_adjust) {
        coord_t                loops_count = std::max(bbox_size.x(), bbox_size.y()) / scaled_flow_spacing + 1;
        Polygons               polygons = offset(expolygon, float(scaled_flow_spacing/2));

        double min_nozzle_diameter = *std::min_element(print_config->nozzle_diameter.values.begin(), print_config->nozzle_diameter.values.end());
        Arachne::WallToolPathsParams input_params;
        input_params.min_bead_width = 0.85 * min_nozzle_diameter;
        input_params.min_feature_size = 0.25 * min_nozzle_diameter;
        input_params.wall_transition_length = 1.0 * min_nozzle_diameter;
        input_params.wall_transition_angle = 10;
        input_params.wall_transition_filter_deviation = 0.25 * min_nozzle_diameter;
        input_params.wall_distribution_count = 1;

        Arachne::WallToolPaths wallToolPaths(polygons, scaled_flow_spacing, scaled_flow_spacing, loops_count, 0, params.layer_height, input_params);

        std::vector<Arachne::VariableWidthLines>    loops = wallToolPaths.getToolPaths();
        std::vector<const Arachne::ExtrusionLine*> all_extrusions;
        for (Arachne::VariableWidthLines& loop : loops) {
            if (loop.empty())
                continue;
            for (const Arachne::ExtrusionLine& wall : loop)
                all_extrusions.emplace_back(&wall);
        }

        // Orca: a forced fill order prints the loops in strictly monotonic depth order so
        // that surfaces broken up by holes or slots cannot hop outward and back inward.
        const bool forced_fill_order = params.fill_order != SurfaceFillOrder::Default;
        if (forced_fill_order) {
            const bool outward = params.fill_order == SurfaceFillOrder::Outward;
            std::stable_sort(all_extrusions.begin(), all_extrusions.end(),
                             [outward](const Arachne::ExtrusionLine *a, const Arachne::ExtrusionLine *b) {
                                 return outward ? a->inset_idx > b->inset_idx : a->inset_idx < b->inset_idx;
                             });
        }

        // Split paths using a nearest neighbor search.
        size_t firts_poly_idx = thick_polylines_out.size();
        Point  last_pos(0, 0);
        for (const Arachne::ExtrusionLine* extrusion : all_extrusions) {
            if (extrusion->empty())
                continue;

            ThickPolyline thick_polyline = Arachne::to_thick_polyline(*extrusion);
            if (extrusion->is_closed)
                thick_polyline.start_at_index(last_pos.nearest_point_index(thick_polyline.points));
            thick_polylines_out.emplace_back(std::move(thick_polyline));
            last_pos = thick_polylines_out.back().last_point();
        }

        // clip the paths to prevent the extruder from getting exactly on the first point of the loop
        // Keep valid paths only.
        size_t j = firts_poly_idx;
        for (size_t i = firts_poly_idx; i < thick_polylines_out.size(); ++i) {
            thick_polylines_out[i].clip_end(this->loop_clipping);
            if (thick_polylines_out[i].is_valid()) {
                if (j < i)
                    thick_polylines_out[j] = std::move(thick_polylines_out[i]);
                ++j;
            }
        }
        if (j < thick_polylines_out.size())
            thick_polylines_out.erase(thick_polylines_out.begin() + int(j), thick_polylines_out.end());

        if (!forced_fill_order)
            reorder_by_shortest_traverse(thick_polylines_out);
    } else {
        Polylines polylines;
        this->_fill_surface_single(params, thickness_layers, direction, expolygon, polylines);
        append(thick_polylines_out, to_thick_polylines(std::move(polylines), scaled_flow_spacing));
    }
}

} // namespace Slic3r
