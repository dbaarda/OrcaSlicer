#include "Exception.hpp"
#include "Print.hpp"

namespace Slic3r {

unsigned int PrintRegion::extruder(FlowRole role) const
{
    switch (role) {
    case frExternalPerimeter: return m_config.outer_wall_filament_id;
    case frPerimeter: return m_config.inner_wall_filament_id;
    case frInfill: return m_config.sparse_infill_filament_id;
    case frSolidInfill: return m_config.internal_solid_filament_id;
    case frTopSolidInfill: return m_config.top_surface_filament_id;
    case frBottomSolidInfill: return m_config.bottom_surface_filament_id;
    // Should external bridges use the same fillament as the bottom infill, or internal bridges?
    case frExternalBridge: return m_config.bottom_surface_filament_id;
    case frInternalBridge: return m_config.internal_solid_filament_id;
    case frSupportMaterial:
    case frSupportMaterialInterface:
    case frSupportTransition: throw Slic3r::InvalidArgument("Cannot get extruder for support roles without an object.");
    default: throw Slic3r::InvalidArgument("Unknown flow role");
    }
}

unsigned int PrintRegion::extruder(const PrintObject& object, FlowRole role) const
{
    const PrintObjectConfig& object_config = object.config();
    switch (role) {
    case frSupportMaterial: return object_config.support_filament;
    case frSupportMaterialInterface: return object_config.support_interface_filament;
    case frSupportTransition: return object_config.support_filament;
    default: return this->extruder(role);
    }
}

float PrintRegion::nozzle_diameter(const PrintConfig& print_config, FlowRole role) const
{
    // Get the configured nozzle_diameter for the extruder associated to the flow role requested.
    // Here this->extruder(role) - 1 may underflow to MAX_INT, but then the get_at() will follback to zero'th element, so everything is all right.
    return print_config.nozzle_diameter.get_at(this->extruder(role) - 1);
}

float PrintRegion::nozzle_diameter(const PrintObject &object, FlowRole role) const
{
    const PrintConfig& print_config        = object.print()->config();
    // Get the configured nozzle_diameter for the extruder associated to the flow role requested.
    // Here this->extruder(role) - 1 may underflow to MAX_INT, but then the get_at() will follback to zero'th element, so everything is all right.
    return print_config.nozzle_diameter.get_at(this->extruder(object, role) - 1);
}

Flow PrintRegion::flow(const PrintObject& object, FlowRole role, double layer_height, bool first_layer) const
{
    const PrintConfig& print_config        = object.print()->config();
    const PrintObjectConfig& object_config = object.config();
    ConfigOptionFloatOrPercent config_width;
    // On the first layer use initial_layer_line_width for all roles if set.
    if (first_layer && print_config.initial_layer_line_width.value > 0) {
        config_width = print_config.initial_layer_line_width;
    } else {
        switch (role) {
        case frExternalPerimeter: config_width = m_config.outer_wall_line_width; break;
        case frPerimeter: config_width = m_config.inner_wall_line_width; break;
        case frInfill:; config_width = m_config.sparse_infill_line_width; break;
        case frSolidInfill: config_width = m_config.internal_solid_infill_line_width; break;
        case frTopSolidInfill: config_width = m_config.top_surface_line_width; break;
        // Should bottom surfaces use top_surface_line_width? Or have its own line_width setting?
        case frBottomSolidInfill: config_width = m_config.internal_solid_infill_line_width; break;
        case frExternalBridge:
        case frInternalBridge: config_width = m_config.bridge_line_width; break;
        case frSupportMaterial:
        case frSupportMaterialInterface:
        case frSupportTransition: config_width = object_config.support_line_width; break;
        default: throw Slic3r::InvalidArgument("Unknown flow role");
        }
    }
    // Fallback to line_width if the role-specific width is zero.
    if (config_width.value == 0)
        config_width = object_config.line_width;
    return Flow::new_from_config_width(role, config_width, this->nozzle_diameter(object, role), float(layer_height));
}

coordf_t PrintRegion::nozzle_dmr_avg(const PrintConfig& print_config) const
{
    // Note this doesn't include bridge flows because they didn't exist when this was created, and it doesn't include support flows because it can't without a
    // PrintObject. Currently bridge flows use the same nozzle as frSolidInfill and frBottomSolidInfill anyway, and it's also unclear in the one place where
    // this is used in LayerRegion.cpp if they should be included or not.
    return (this->nozzle_diameter(print_config, frExternalPerimeter) +
            this->nozzle_diameter(print_config, frPerimeter) +
            this->nozzle_diameter(print_config, frInfill) +
            this->nozzle_diameter(print_config, frSolidInfill) +
            this->nozzle_diameter(print_config, frTopSolidInfill) +
            this->nozzle_diameter(print_config, frBottomSolidInfill)) / 6.0;
}

coordf_t PrintRegion::bridging_height_avg(const PrintConfig& print_config) const
{
    // Note this only uses the external bridge nozzle diameter because the only place this is used is for supports to calculate the correct support distance,
    // and we never support internal bridges.
    coordf_t nozzle_diameter = this->nozzle_diameter(print_config, frExternalBridge);
    return m_config.get_abs_value("bridge_line_width", nozzle_diameter);
}

void PrintRegion::collect_object_printing_extruders(const PrintConfig& print_config,
                                                    const PrintRegionConfig& region_config,
                                                    const bool has_brim,
                                                    std::vector<unsigned int>& object_extruders)
{
    // These checks reflect the same logic used in the GUI for enabling/disabling extruder selection fields.
    // BBS
    auto num_extruders    = (int) print_config.filament_diameter.size();
    auto emplace_extruder = [num_extruders, &object_extruders](int extruder_id) {
        int i = std::max(0, extruder_id - 1);
        object_extruders.emplace_back((i >= num_extruders) ? 0 : i);
    };
    if (region_config.wall_loops.value > 0 || has_brim) {
        emplace_extruder(region_config.outer_wall_filament_id);
        if (region_config.wall_loops.value > 1)
            emplace_extruder(region_config.inner_wall_filament_id);
    }
    if (region_config.sparse_infill_density.value > 0)
        emplace_extruder(region_config.sparse_infill_filament_id);
    if (region_config.sparse_infill_density.value > 0 || region_config.top_shell_layers.value > 0 ||
        region_config.bottom_shell_layers.value > 0)
        emplace_extruder(region_config.internal_solid_filament_id);
    if (region_config.top_shell_layers.value > 0)
        emplace_extruder(region_config.top_surface_filament_id);
    if (region_config.bottom_shell_layers.value > 0)
        emplace_extruder(region_config.bottom_surface_filament_id);
}

void PrintRegion::collect_object_printing_extruders(const Print& print, std::vector<unsigned int>& object_extruders) const
{
    // PrintRegion, if used by some PrintObject, shall have all the extruders set to an existing printer extruder.
    // If not, then there must be something wrong with the Print::apply() function.
#ifndef NDEBUG
    // BBS
    auto num_extruders = int(print.config().filament_diameter.size());
    assert(this->config().outer_wall_filament_id <= num_extruders);
    assert(this->config().inner_wall_filament_id <= num_extruders);
    assert(this->config().sparse_infill_filament_id <= num_extruders);
    assert(this->config().internal_solid_filament_id <= num_extruders);
    assert(this->config().top_surface_filament_id <= num_extruders);
    assert(this->config().bottom_surface_filament_id <= num_extruders);
#endif
    collect_object_printing_extruders(print.config(), this->config(), print.has_brim(), object_extruders);
}

} // namespace Slic3r
