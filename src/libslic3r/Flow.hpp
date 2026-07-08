#ifndef slic3r_Flow_hpp_
#define slic3r_Flow_hpp_

#include "libslic3r.h"
#include "Config.hpp"
#include "Exception.hpp"
#include "PrintConfig.hpp"

namespace Slic3r {

class PrintObject;

// Each FlowRole value identifies a distinct extruder.
enum FlowRole {
    frNone,
    frExternalPerimeter,
    frPerimeter,
    frInfill,
    frSolidInfill,
    frTopSolidInfill,
    frBottomSolidInfill,
    frSupportMaterial,
    frSupportMaterialInterface,
    frSupportTransition, // BBS
};

class FlowError : public Slic3r::InvalidArgument
{
public:
    FlowError(const std::string& what_arg) : Slic3r::InvalidArgument(what_arg) {}
    FlowError(const char* what_arg) : Slic3r::InvalidArgument(what_arg) {}
};

class FlowErrorNegativeSpacing : public FlowError
{
public:
    FlowErrorNegativeSpacing();
};

class FlowErrorNegativeHeight : public FlowError
{
public:
    FlowErrorNegativeHeight();
};

class FlowErrorHeightTooLarge : public FlowError
{
public:
    FlowErrorHeightTooLarge();
};

class FlowErrorMissingVariable : public FlowError
{
public:
    FlowErrorMissingVariable(const std::string& what_arg) : FlowError(what_arg) {}
};

class Flow
{
public:
    // These are constants and helper functions used for calculating areas and overlaps of lines.
    //
    // See https://manual.slic3r.org/advanced/flow-math for details.
    //
    // To extend this to work for narrow gapfill-lines with flows so small that width < height, we switch to a narrow-elipse cross-section
    // model where height stays at the layer height and width shrinks as the flow gets smaller.
    //
    // For both rounded-rectangle and narrow-elipse models, overlap_factor=0 means adjacent lines just touch with spacing=width, and
    // overlap_factor=1 is when the flow cross-section area is "squished" into a narrower spacing to completely fill the available
    // rectangular height*spacing area, equivalent to a 100% flow_ratio.
    //
    // For circular bridge-lines, they can be "squished" horizontally and vertically into a rectangular vertical-spacing *
    // horizontal-spacing available area. We define the horizontal and vertical overlap_factor=0 as when adjacent lines/layers just touch,
    // so spacing=diameter and height=diameter, and overlap_factor=1 as when the spacing area is a square with the same area as the circular
    // cross-section. Note this this means 100% flow ratio to completely fill the spacing*height area requires both horizontal and vertical
    // overlap_factor=1.
    //
    // Overlap distance for overlap_factor=1 of a rounded-rectangle line, overlap = height * RRECT_OVERLAP
    static constexpr float RRECT_OVERLAP = 1.0 - M_PI_4;
    // Overlap distance for overlap_factor=1 of a circular bridge line, overlap = diameter * BRIDGE_OVERLAP
    static constexpr float BRIDGE_OVERLAP = 1.0 - std::sqrt(M_PI_4);

    // Simple helper functions for circle areas.
    static constexpr float area2dmr(float a) { return std::sqrt(a / M_PI_4); }
    static constexpr float dmr2area(float d) { return d * d * M_PI_4; }

    // Get the rounded-rectangle normal-line spacing.
    static constexpr float rrect_spacing(float width, float height, float overlap_factor = 1.0)
    {
        // overlap_factor cannot be so high it gives negative spacing.
        assert((overlap_factor * RRECT_OVERLAP) < (width / std::min(width, height)));
        return width - overlap_factor * RRECT_OVERLAP * std::min(width, height);
    }

    // Get the rounded-rectangle normal-line width.
    static constexpr float rrect_width(float spacing, float height, float overlap_factor = 1.0)
    {
        float width = spacing + overlap_factor * RRECT_OVERLAP * height;
        // If the rounded-rectangle width < height, use the narrow-elipse model.
        return width < height ? spacing / (1 - overlap_factor * RRECT_OVERLAP) : width;
    }

    // Get the rounded-rectangle normal-line cross-section area.
    static constexpr float rrect_area(float width, float height) { return rrect_spacing(width, height) * height; }

    // Get the rounded-rectangle normal-line overlap_factor.
    static constexpr float rrect_overlap_factor(float width, float spacing, float height)
    { return (width - spacing) / (RRECT_OVERLAP * std::min(width, height)); }

    // Get the circular bridge-line spacing.
    static constexpr float bridge_spacing(float diameter, float overlap_factor = 0.0)
    {
        // Overlap cannot be so high it gives negative spacing.
        assert(overlap_factor * BRIDGE_OVERLAP < 1);
        return diameter - overlap_factor * BRIDGE_OVERLAP * diameter;
    }

    // Get the circular bridge-line width.
    static constexpr float bridge_width(float spacing, float overlap_factor = 0.0)
    { return spacing / (1 - overlap_factor * BRIDGE_OVERLAP); }

    // Get the circular bridge-line cross-section area.
    static constexpr float bridge_area(float diameter) { return dmr2area(diameter); }

    // Get the circular bridge-line horizontal overlap_factor.
    static constexpr float bridge_overlap_factor(float diameter, float spacing)
    { return (diameter - spacing) / (BRIDGE_OVERLAP * diameter); }

    // Get the circular bridge-line height.
    static constexpr float bridge_height(float diameter, float vertical_overlap_factor = 0.0)
    {
        // vertical_overlap_factor cannot be so high it gives negative height.
        assert(vertical_overlap_factor * BRIDGE_OVERLAP < 1);
        return diameter - vertical_overlap_factor * BRIDGE_OVERLAP * diameter;
    }

    // Get the circular bridge-line vertical overlap_factor.
    static constexpr float bridge_vertical_overlap_factor(float diameter, float height)
    { return (diameter - height) / (BRIDGE_OVERLAP * diameter); }

    // For normal-line flows with width >= height using a rounded-rectangle
    // cross-section extrusion model;
    //
    // * m_width is the full width of the rounded-rectangle cross-section.
    // * m_height is the vertical height of the rounded-rectangle and layer.
    // * m_spacing is the horizontal distance between adjacent possibly overlapping lines.
    //
    // For gapfill-line flows with width < height using a narrow-elipse
    // cross-section extrusion model;
    //
    // * m_width is the width of the narrow-elipse cross-section.
    // * m_height is the vertical height of the narrow-elipse and layer.
    // * m_spacing is the horizontal distance between adjacent possibly overlapping lines.
    //
    // For bridge-line flows using a circular cross-section extrusion model;
    //
    // * m_width is the diameter of the circular cross-section.
    // * m_height is the vertical distance between possibly overlapping layers of bridge lines.
    // * m_spacing is the horizontal distance between adjacent possibly overlapping bridge lines.
    //
    // Although currently bridge-lines are assumed to never vertically overlap so m_width = m_height = diameter, always use m_width or
    // diameter() for the diameter incase support for vertical overlapping is added in the future.
    Flow() = default;
    // Initialize a standard normal-line rounded-rectangle or narrow-elipse flow.
    Flow(float width, float height, float nozzle_diameter) : Flow(width, height, rrect_spacing(width, height), nozzle_diameter, false) {}
    // Initialize a standard bridge-line circular flow.
    Flow(float diameter, float nozzle_diameter) : Flow(diameter, diameter, diameter, nozzle_diameter, true) {}

    // Vertical spacing between extrusion layers.
    float height() const { return m_height; }
    // Width of the flow's extrusion model cross-section.
    float width() const { return m_width; }
    coord_t scaled_width() const { return coord_t(scale_(m_width)); }
    // Horizontal spacing between the extrusion centerlines.
    float spacing() const { return m_spacing; }
    coord_t scaled_spacing() const { return coord_t(scale_(m_spacing)); }
    // Nozzle diameter.
    float nozzle_diameter() const { return m_nozzle_diameter; }
    // Is it a circular bridge-line flow?
    bool bridge() const { return m_bridge; }
    // Is this a narrow-elipse gapfill-line flow?
    bool narrow() const { return m_width < m_height; }

    // Cross section area of the extrusion.
    double mm3_per_mm() const { return m_bridge ? bridge_area(m_width) : rrect_area(m_width, m_height); }
    // The diameter of an equivalent volume circular flow.
    float diameter() const { return m_bridge ? m_width : area2dmr(mm3_per_mm()); }
    // The horizontal overlap_factor.
    float overlap_factor() const
    { return m_bridge ? bridge_overlap_factor(m_width, m_spacing) : rrect_overlap_factor(m_width, m_spacing, m_height); }
    // The vertical overlap_factor. For normal lines we return 1 to match 100% flow ratio.
    float vertical_overlap_factor() const { return m_bridge ? bridge_vertical_overlap_factor(m_width, m_height) : 1.0; }
    // The flow_ratio of the line cross-section area to available height*spacing area.
    float flow_ratio() const { return mm3_per_mm() / (m_spacing * m_height); }

    // Note flows with only spacing different are considered equal. Is this correct?
    inline bool operator==(const Flow& rhs) const
    { return m_width == rhs.m_width && m_height == rhs.m_height && m_nozzle_diameter == rhs.m_nozzle_diameter && m_bridge == rhs.m_bridge; }

    inline bool operator!=(const Flow& rhs) const { return !operator==(rhs); }

    bool operator<(const Flow& rhs) const { return this->mm3_per_mm() < rhs.mm3_per_mm(); }

    // Create a modified flow with a different width while maintaining the same overlap_factors. For bridge-lines this scales spacing and
    // height. For normal-lines it preserves height and adjusts spacing.
    Flow with_width(float width) const;

    // Create a modified flow with a different height while maintaining overlap_factors. For bridge lines this scales width and spacing. For
    // normal lines it preserves spacing and adjusts width.
    Flow with_height(float height) const;

    // Create a modified flow with a different spacing while maintaining overlap_factors. For bridge lines this scales width and height. For
    // normal lines it preserves height and adjusts width.
    Flow with_spacing(float spacing) const;

    // Create a modified flow that is a bridge version of a normal line flow. Without specifying diameter it returns a flow with the same
    // cross-section area.
    Flow as_bridge(float diameter = -1.0) const;

    // Create a flow for a FlowRole, width option, nozzle_diameter, and height.
    static Flow new_from_config_width(
        FlowRole role, const ConfigOptionFloatOrPercent& width, float nozzle_diameter, float height, bool bridge = false);

    // Sane extrusion width default based on nozzle diameter.
    // The defaults were derived from manual Prusa MK3 profiles.
    static float auto_extrusion_width(FlowRole role, float nozzle_diameter, bool bridge = false);

    // Extrusion width from full config, taking into account the defaults (when set to zero) and ratios (percentages). Precise value depends
    // on layer index (1st layer vs. other layers vs. variable layer height), on active extruder etc. Therefore the value calculated by this
    // function shall be used as a hint only.
    static double extrusion_width(const std::string& opt_key,
                                  const ConfigOptionResolver& config,
                                  const unsigned int first_printing_extruder = 0);

    // Get the extruder_id for a flow role from configs.
    static size_t get_extruder(const PrintObjectConfig& object_config, const PrintRegionConfig& region_config, FlowRole role);

    // Get the nozzle diameter for an extruder from configs.
    static float get_nozzle_diameter(const PrintConfig& print_config, size_t extruder);

    // Get the nozzle diameter for a flow role from configs.
    static float get_nozzle_diameter(const PrintConfig& print_config,
                                     const PrintObjectConfig& object_config,
                                     const PrintRegionConfig& region_config,
                                     FlowRole role);

    // Get the line_width for a flow role from configs. If nozzle_diameter is not provided it will be fetched from configs.
    static float get_line_width(const PrintConfig& print_config,
                                const PrintObjectConfig& object_config,
                                const PrintRegionConfig& region_config,
                                FlowRole role,
                                bool bridge           = false,
                                bool first_layer      = false,
                                float nozzle_diameter = 0.0);

    static Flow new_from_role(const PrintConfig& print_config,
                              const PrintObjectConfig& object_config,
                              const PrintRegionConfig& region_config,
                              FlowRole role,
                              float height,
                              bool bridge      = false,
                              bool first_layer = false);

private:
    Flow(float width, float height, float spacing, float nozzle_diameter, bool bridge)
        : m_width(width), m_height(height), m_spacing(spacing), m_nozzle_diameter(nozzle_diameter), m_bridge(bridge)
    {
        if (m_height < 0) {
            throw FlowErrorNegativeHeight();
        }
        // nozzle_diameter<=0 is used for "don't care" so assume height is valid.
        if (m_nozzle_diameter > 0 && m_height > m_nozzle_diameter + EPSILON) {
            throw FlowErrorHeightTooLarge();
        }
        if (m_spacing < 0) {
            throw FlowErrorNegativeSpacing();
        }
        // Negative width should never happen, unlike negative spacing or negative height, which might happen for invalid settings values.
        assert(m_width > 0);
    }

    float m_width{0};
    float m_height{0};
    float m_spacing{0};
    float m_nozzle_diameter{0};
    bool m_bridge{false};
};

} // namespace Slic3r
#endif
