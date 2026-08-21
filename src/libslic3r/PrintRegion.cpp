#include "Exception.hpp"
#include "Print.hpp"

namespace Slic3r {
/*
unsigned int PrintRegion::extruder(const PrintObject& object, FlowRole role) const
{
    return object.extruder(*this, role);
}

float PrintRegion::nozzle_diameter(const PrintObject &object, FlowRole role) const
{
    return object.nozzle_diameter(*this, role);
}

Flow PrintRegion::flow(const PrintObject& object, FlowRole role, double layer_height, bool first_layer, bool bridge) const
{
    return object.flow(*this, role, layer_height, bridge, first_layer);
}
*/
coordf_t PrintRegion::nozzle_dmr_avg(const PrintObject& object) const
{    
    // Note this doesn't include bridge or support extruders because it didn't in the past. Currently bridge flows use the same nozzle as
    // frSolidInfill and frBottomSolidInfill anyway, and it's also unclear in the one place where this is used in LayerRegion.cpp if supports should
    // be included or not.
    return (object.nozzle_diameter(*this, frExternalPerimeter) +
            object.nozzle_diameter(*this, frPerimeter) +
            object.nozzle_diameter(*this, frInfill) +
            object.nozzle_diameter(*this, frSolidInfill) +
            object.nozzle_diameter(*this, frTopSolidInfill) +
            object.nozzle_diameter(*this, frBottomSolidInfill)) / 6.0;
}

coordf_t PrintRegion::bridging_height_avg(const PrintObject& object) const
{
    // Note this returns only the bottom surface bridge flow's height. The only place this is used is
    // for supports to calculate the correct support distance, and we never support internal bridges.
    return object.flow(*this, frBottomSolidInfill, -1.0, true).height();
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
