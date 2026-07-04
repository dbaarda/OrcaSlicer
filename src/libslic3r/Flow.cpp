#include "Flow.hpp"
#include "I18N.hpp"
#include <cmath>
#include <assert.h>

#include <boost/algorithm/string/predicate.hpp>

// Mark string for localization and translate.
#define L(s) Slic3r::I18N::translate(s)

namespace Slic3r {

FlowErrorNegativeSpacing::FlowErrorNegativeSpacing() : FlowError("Flow::spacing() is negative. Did you set some extrusion width too small?")
{}

FlowErrorNegativeHeight::FlowErrorNegativeHeight() : FlowError("Flow::height() is negative. Did you set some line setting too small?") {}

FlowErrorHeightTooLarge::FlowErrorHeightTooLarge()
    : FlowError("Flow::height() > Flow::nozzle_diameter(). Did you set layer height too large?")
{}

Flow Flow::with_width(float width) const
{
    float scale = width / m_width;
    return m_bridge ? Flow(width, scale * m_height, scale * m_spacing, m_nozzle_diameter, m_bridge) :
                      Flow(width, m_height, rrect_spacing(width, m_height, overlap_factor()), m_nozzle_diameter, m_bridge);
}

Flow Flow::with_height(float height) const
{
    float scale = height / m_height;
    return m_bridge ? Flow(scale * m_width, height, scale * m_spacing, m_nozzle_diameter, m_bridge) :
                      Flow(rrect_width(m_spacing, height, overlap_factor()), height, m_spacing, m_nozzle_diameter, m_bridge);
}

Flow Flow::with_spacing(float spacing) const
{
    float scale = spacing / m_spacing;
    return m_bridge ? Flow(scale * m_width, scale * m_height, spacing, m_nozzle_diameter, m_bridge) :
                      Flow(rrect_width(spacing, m_height, overlap_factor()), m_height, spacing, m_nozzle_diameter, m_bridge);
}

Flow Flow::as_bridge(float dmr) const { return Flow(dmr > 0.0 ? dmr : diameter(), m_nozzle_diameter); }

// This static method returns a sane extrusion width default.
float Flow::auto_extrusion_width(FlowRole role, float nozzle_diameter, bool bridge)
{
    if (bridge)
        return nozzle_diameter;
    switch (role) {
    case frSupportMaterial:
    case frSupportMaterialInterface:
    case frSupportTransition: return nozzle_diameter;
    default: return 1.125f * nozzle_diameter;
    }
}

// Used by the Flow::extrusion_width() function to provide hints to the user on default extrusion width values,
// and to provide reasonable values to the PlaceholderParser.
static inline FlowRole opt_key_to_flow_role(const std::string& opt_key)
{
    if (opt_key == "inner_wall_line_width" ||
        // or all the defaults:
        opt_key == "line_width" || opt_key == "initial_layer_line_width")
        return frPerimeter;
    else if (opt_key == "outer_wall_line_width")
        return frExternalPerimeter;
    else if (opt_key == "sparse_infill_line_width")
        return frInfill;
    else if (opt_key == "internal_solid_infill_line_width")
        return frSolidInfill;
    else if (opt_key == "bridge_line_width")
        return frSolidInfill;
    else if (opt_key == "top_surface_line_width")
        return frTopSolidInfill;
    else if (opt_key == "support_line_width")
        return frSupportMaterial;
    else
        throw Slic3r::RuntimeError("opt_key_to_flow_role: invalid argument");
};

// Used to provide hints to the user on default extrusion width values, and to provide reasonable values to the PlaceholderParser.
double Flow::extrusion_width(const std::string& opt_key, const ConfigOptionResolver& config, const unsigned int first_printing_extruder)
{
    auto opt_nozzle_diameters   = config.option_throw<ConfigOptionFloats>("nozzle_diameter");
    const float nozzle_diameter = float(opt_nozzle_diameters->get_at(first_printing_extruder));

    double value = config.option_throw<ConfigOptionFloatOrPercent>(opt_key)->get_abs_value(nozzle_diameter);
    // for non-bridge widths, default back to the "line_width" setting.
    if (value <= 0.) {
        if (opt_key == "bridge_line_width")
            return nozzle_diameter;
        value = config.option_throw<ConfigOptionFloatOrPercent>("line_width")->get_abs_value(nozzle_diameter);
    }

    // If the value still is zero, calculate a sane default width.
    return (value <= 0.) ? auto_extrusion_width(opt_key_to_flow_role(opt_key), nozzle_diameter) : value;
}

// This constructor builds a Flow object from an extrusion width config setting
// and other context properties.
Flow Flow::new_from_config_width(FlowRole role, const ConfigOptionFloatOrPercent& width, float nozzle_diameter, float height, bool bridge)
{
    if (nozzle_diameter <= 0)
        throw Slic3r::InvalidArgument("Invalid nozzle_diameter supplied to new_from_config_width()");
    if (height <= 0)
        throw Slic3r::InvalidArgument("Invalid flow height supplied to new_from_config_width()");

    // If user left option to 0, use a sane default width.
    float line_width = float(width.value <= 0. ? auto_extrusion_width(role, nozzle_diameter, bridge) : width.get_abs_value(nozzle_diameter));
    return bridge ? Flow(line_width, nozzle_diameter) : Flow(line_width, height, nozzle_diameter);
}

size_t Flow::get_extruder(const PrintObjectConfig& object_config, const PrintRegionConfig& region_config, FlowRole role)
{
    switch (role) {
    case frExternalPerimeter: return region_config.outer_wall_filament_id;
    case frPerimeter: return region_config.inner_wall_filament_id;
    case frInfill: return region_config.sparse_infill_filament_id;
    case frSolidInfill: return region_config.internal_solid_filament_id;
    case frTopSolidInfill: return region_config.top_surface_filament_id;
    case frBottomSolidInfill: return region_config.bottom_surface_filament_id;
    case frSupportMaterial: return object_config.support_filament;
    case frSupportMaterialInterface: return object_config.support_interface_filament;
    case frSupportTransition: return object_config.support_filament;
    default: throw Slic3r::InvalidArgument("Unknown flow role");
    }
}

float Flow::get_nozzle_diameter(const PrintConfig& print_config, size_t extruder)
{
    // Get the configured nozzle_diameter for the extruder associated to the flow role requested.
    // Here this->extruder(role) - 1 may underflow to MAX_INT, but then the get_at() will follback to zero'th element, so everything is all right.
    return print_config.nozzle_diameter.get_at(extruder - 1);
}

float Flow::get_nozzle_diameter(const PrintConfig& print_config,
                                const PrintObjectConfig& object_config,
                                const PrintRegionConfig& region_config,
                                FlowRole role)
{ return get_nozzle_diameter(print_config, get_extruder(object_config, region_config, role)); }

// Get the base line width config option for a flow role without bridge or first_line modifiers from configs.
const auto& get_base_config_width(const PrintObjectConfig& object_config, const PrintRegionConfig& region_config, FlowRole role)
{
    switch (role) {
    case frExternalPerimeter: return region_config.outer_wall_line_width;
    case frPerimeter: return region_config.inner_wall_line_width;
    case frInfill:; return region_config.sparse_infill_line_width;
    case frSolidInfill: return region_config.internal_solid_infill_line_width;
    case frTopSolidInfill: return region_config.top_surface_line_width;
    // Should bottom surfaces use top_surface_line_width? Or have its own line_width setting?
    case frBottomSolidInfill: return region_config.internal_solid_infill_line_width;
    case frSupportMaterial:
    case frSupportMaterialInterface:
    case frSupportTransition: return object_config.support_line_width;
    default: throw Slic3r::InvalidArgument("Unknown flow role");
    }
}

// Get the line width config option for a flow role with bridge or first_line modifiers from configs.
const auto& get_config_width(const PrintConfig& print_config,
                             const PrintObjectConfig& object_config,
                             const PrintRegionConfig& region_config,
                             FlowRole role,
                             bool bridge,
                             bool first_layer)
{
    if (bridge) {
        return region_config.bridge_line_width;
    } else if (first_layer && print_config.initial_layer_line_width.value > 0) {
        // On the first layer use initial_layer_line_width for all roles if set.
        return print_config.initial_layer_line_width;
    }
    const auto& config_width = get_base_config_width(object_config, region_config, role);
    return config_width.value > 0.0 ? config_width : object_config.line_width;
}

// Get the line_width for a flow role from configs.
float Flow::get_line_width(const PrintConfig& print_config,
                           const PrintObjectConfig& object_config,
                           const PrintRegionConfig& region_config,
                           FlowRole role,
                           bool bridge,
                           bool first_layer,
                           float nozzle_diameter)
{
    if (nozzle_diameter <= 0.0)
        nozzle_diameter = get_nozzle_diameter(print_config, object_config, region_config, role);
    const auto& config_width = get_config_width(print_config, object_config, region_config, role, bridge, first_layer);
    // Return the config value, falling back to the auto defaults if it is zero.
    return (config_width.value > 0.0) ? config_width.get_abs_value(nozzle_diameter) : auto_extrusion_width(role, nozzle_diameter, bridge);
}

Flow Flow::new_from_role(const PrintConfig& print_config,
                         const PrintObjectConfig& object_config,
                         const PrintRegionConfig& region_config,
                         FlowRole role,
                         float height,
                         bool bridge,
                         bool first_layer)
{
    float nozzle_diameter = get_nozzle_diameter(print_config, object_config, region_config, role);
    float width           = get_line_width(print_config, object_config, region_config, role, bridge, first_layer, nozzle_diameter);
    return bridge ? Flow(width, nozzle_diameter) : Flow(width, height, nozzle_diameter);
}

} // namespace Slic3r
