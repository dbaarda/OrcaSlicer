#include "Flow.hpp"
#include "I18N.hpp"
#include "Print.hpp"
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

// This static method returns a sane extrusion width default.
float Flow::auto_extrusion_width(FlowRole role, float nozzle_diameter)
{
    switch (role) {
    case frSupportMaterial:
    case frSupportMaterialInterface:
    case frSupportTransition:
    case frTopSolidInfill: return nozzle_diameter;
    default:
    case frExternalPerimeter:
    case frPerimeter:
    case frSolidInfill:
    case frInfill: return 1.125f * nozzle_diameter;
    }
}

// Used by the Flow::extrusion_width() funtion to provide hints to the user on default extrusion width values,
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
    if (value == 0.) {
        if (opt_key == "bridge_line_width") {
            // For bridge_line_width, default to the nozzle_diameter.
            value = nozzle_diameter;
        } else {
            // for other widths, default back to the "line_width" setting.
            value = config.option_throw<ConfigOptionFloatOrPercent>("line_width")->get_abs_value(nozzle_diameter);
        }
    }

    // If the value still is zero, calculate a sane default width.
    return (value == 0.) ? auto_extrusion_width(opt_key_to_flow_role(opt_key), nozzle_diameter) : value;
}

// This constructor builds a Flow object from an extrusion width config setting
// and other context properties.
Flow Flow::new_from_config_width(FlowRole role, const ConfigOptionFloatOrPercent& width, float nozzle_diameter, float height)
{
    if (height <= 0)
        throw Slic3r::InvalidArgument("Invalid flow height supplied to new_from_config_width()");

    float w;
    if (!width.percent && width.value <= 0.) {
        // If user left option to 0, calculate a sane default width.
        w = auto_extrusion_width(role, nozzle_diameter);
    } else {
        // If user set a manual value, use it.
        w = float(width.get_abs_value(nozzle_diameter));
    }

    return Flow(w, height, rrect_spacing(w, height), nozzle_diameter, false);
}

Flow support_material_flow(const PrintObject* object, float layer_height)
{
    return Flow::new_from_config_width(frSupportMaterial,
                                       // The width parameter accepted by new_from_config_width is of type ConfigOptionFloatOrPercent, the
                                       // Flow class takes care of the percent to value substitution.
                                       (object->config().support_line_width.value > 0) ? object->config().support_line_width :
                                                                                         object->config().line_width,
                                       // if object->config().support_filament == 0 (which means to not trigger tool change, but use the
                                       // current extruder instead), get_at will return the 0th component.
                                       float(object->print()->config().nozzle_diameter.get_at(object->config().support_filament - 1)),
                                       (layer_height > 0.f) ? layer_height : float(object->config().layer_height.value));
}

Flow support_material_1st_layer_flow(const PrintObject* object, float layer_height)
{
    const PrintConfig& print_config = object->print()->config();
    const auto& width               = (print_config.initial_layer_line_width.value > 0) ? print_config.initial_layer_line_width :
                                                                                          object->config().support_line_width;
    return Flow::new_from_config_width(frSupportMaterial,
                                       // The width parameter accepted by new_from_config_width is of type ConfigOptionFloatOrPercent, the
                                       // Flow class takes care of the percent to value substitution.
                                       (width.value > 0) ? width : object->config().line_width,
                                       float(print_config.nozzle_diameter.get_at(object->config().support_filament - 1)),
                                       (layer_height > 0.f) ? layer_height : float(print_config.initial_layer_print_height.value));
}

Flow support_material_interface_flow(const PrintObject* object, float layer_height)
{
    return Flow::new_from_config_width(frSupportMaterialInterface,
                                       // The width parameter accepted by new_from_config_width is of type ConfigOptionFloatOrPercent, the
                                       // Flow class takes care of the percent to value substitution.
                                       (object->config().support_line_width > 0) ? object->config().support_line_width :
                                                                                   object->config().line_width,
                                       // if object->config().support_interface_filament == 0 (which means to not trigger tool change, but
                                       // use the current extruder instead), get_at will return the 0th component.
                                       float(object->print()->config().nozzle_diameter.get_at(object->config().support_interface_filament -
                                                                                              1)),
                                       (layer_height > 0.f) ? layer_height : float(object->config().layer_height.value));
}

} // namespace Slic3r
