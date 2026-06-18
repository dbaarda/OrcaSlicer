#ifndef slic3r_Flow_hpp_
#define slic3r_Flow_hpp_

#include "libslic3r.h"
#include "Config.hpp"
#include "Exception.hpp"
#include "ExtrusionEntity.hpp"

namespace Slic3r {

class PrintObject;

// These are constants used for calculating areas and overlaps of lines.
//
// See https://manual.slic3r.org/advanced/flow-math for details.
//
// To extend this to work for narrow gapfill lines with flows so small that
// width < height, we switch to a narrow-elipse cross-section model where
// height stays at the layer height and width shrinks as the flow gets
// smaller.
//
// For both rounded-rectangle and narrow-elipse models, overlap_factor=0
// means adjacent lines just touch with spacing=width, and overlap-factor=1
// is when the flow cross-section area is "squished" into a narrower spacing
// to completely fill the available rectangular height*spacing area,
// equivalent to a 100% flow_ratio.
//
// For circular bridge lines, they can be "squished" horizontally and
// vertically into a rectangular vertical-spacing * horizontal-spacing
// available area. We define the horizontal and vertical overlap_factor=0 as
// when adjacent lines/layers just touch, so spacing=diameter and
// height=diameter, and overlap_factor=1 as when the spacing area is a square
// with the same area as the circular cross-section. Note this this means
// 100% flow ratio to completely fill the spacing*height area requires both
// horizontal and vertical overlap_factor=1.
//
// Area of a circle, area = diameter^2 * PI_4.
static constexpr double PI_4 = PI / 4;
// Length of a square's side with the same area as a circle, length = diameter * SQRTPI_2
static constexpr double SQRTPI_2 = std::sqrt(PI_4);
// Overlap distance for overlap-factor=1 of a rounded-rectangle line, overlap = height * RRLINE_OVERLAP
static constexpr double RRLINE_OVERLAP = 1. - PI_4;
// Overlap distance for overlap-factor=1 of a circular bridge line, overlap = diameter * BRIDGE_OVERLAP
static constexpr double BRIDGE_OVERLAP = 1. - SQRTPI_2;

// Simple helper functions.
static constexpr float area2dmr(float a) { return std::sqrt(a / PI_4); }
static constexpr float dmr2area(float d) { return d * d * PI_4; }

// Get the rounded-rectangle or narrow-elipse spacing from width, height, and
// overlap-factor.
static constexpr float rrect_spacing(float width, float height, float overlap = 1.0)
{
    // Overlap cannot be so high it gives negative spacing.
    assert(overlap * RRLINE_OVERLAP < width / std::min(width, height));
    return width - overlap * RRLINE_OVERLAP * std::min(width, height);
}

// Get the rounded-rectangle or narrow-elipse width from spacing, height, and
// overlap-factor.
static constexpr float rrect_width(float spacing, float height, float overlap = 1.0)
{
    float width = spacing + overlap * RRLINE_OVERLAP * height;
    // If the rounded-rectangle width < height, use the eliptical model.
    return width < height ? spacing / (1 - overlap * RRLINE_OVERLAP) : width;
}

// Get the rounded-rectangle or narrow-elipse cross-section area from width,
// and height.
static constexpr float rrect_area(float width, float height) { return rrect_spacing(width, height) * height; }

// Get the rounded-rectangle or narrow-elipse flow-ratio from width, spacing,
// and height.
static constexpr float rrect_flow(float width, float spacing, float height) { return rrect_spacing(width, height) / spacing; }

// Get the rounded-rectangle or narrow-elipse overlap-factor from width,
// spacing, and height.
static constexpr float rrect_overlap(float width, float spacing, float height)
{ return (width - spacing) / (RRLINE_OVERLAP * std::min(width, height)); }

// Get the circular-bridge spacing from diameter and horizontal overlap-factor.
static constexpr float bridge_spacing(float diameter, float overlap = 0.0)
{
    // Overlap cannot be so high it gives negative spacing.
    assert(overlap * BRIDGE_OVERLAP < 1);
    return diameter - overlap * BRIDGE_OVERLAP * diameter;
}

// Get the circular-bridge width from spacing and horizontal overlap-factor.
static constexpr float bridge_width(float spacing, float overlap = 0.0) { return spacing / (1 - BRIDGE_OVERLAP * overlap); }

// Get the circular-bridge area from diameter.
static constexpr float bridge_area(float diameter) { return dmr2area(diameter); }

// Get the circular-bridge flow-ratio from diameter, spacing, and height.
static constexpr float bridge_flow(float diameter, float spacing, float height) { return bridge_area(diameter) / (spacing * height); }

// Get the circular-bridge horizontal overlap-factor from diameter and spacing.
static constexpr float bridge_overlap(float diameter, float spacing) { return (diameter - spacing) / (BRIDGE_OVERLAP * diameter); }

// Get the circular-bridge height from diameter and vertical overlap-factor.
static constexpr float bridge_height(float diameter, float vertoverlap = 0.0)
{
    // Vertoverlap cannot be so high it gives negative height.
    assert(vertoverlap * BRIDGE_OVERLAP < 1);
    return diameter - vertoverlap * BRIDGE_OVERLAP * diameter;
}

// Get the circular-bridge vertical overlap-factor from diameter and height.
static constexpr float bridge_vertoverlap(float diameter, float height) { return (diameter - height) / (BRIDGE_OVERLAP * diameter); }

enum FlowRole {
    frExternalPerimeter,
    frPerimeter,
    frInfill,
    frSolidInfill,
    frTopSolidInfill,
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

class FlowErrorNegativeFlow : public FlowError
{
public:
    FlowErrorNegativeFlow();
};

class FlowErrorMissingVariable : public FlowError
{
public:
    FlowErrorMissingVariable(const std::string& what_arg) : FlowError(what_arg) {}
};

class Flow
{
public:
    // For normal non-bridging flows with width >= height using a
    // rounded-rectangle extrusion cross-section model;
    //
    // * m_width is the full width of the rounded-rectangle cross-section.
    // * m_height is the vertical height of the rounded-rectangle and layer.
    // * m_spacing is the horizontal distance between adjacent possibly
    //   overlapping lines.
    //
    // For non-bridging narrow gapfill flows with width < height using a
    // narrow-elipse cross-section model;
    //
    // * m_width is the width of the narrow-elipse cross-section.
    // * m_height is the vertical height of the narrow-elipse and layer.
    // * m_spacing is the horizontal distance between adjacent possibly
    //   overlapping lines.
    //
    // For bridging flows using a circular extrusion cross-section model;
    //
    // * m_width is the diameter of the circular cross section.
    // * m_height is the vertical distance between possibly overlapping
    //   layers of bridge lines.
    // * m_spacing is the horizontal distance between adjacent possibly
    //   overlapping bridge lines.
    //
    // Although currently bridge lines are assumed to never vertically
    // overlap so m_width = m_height = diameter, always use m_width or
    // diameter() for the diameter in case support for vertical overlapping
    // is added in the future.
    Flow() = default;
    // Initialize a standard rounded-rectangle or narrow elipse line flow.
    Flow(float width, float height, float nozzle_diameter) : Flow(width, height, rrect_spacing(width, height), nozzle_diameter, false) {}
    // Initialize a standard circular-bridge line flow.
    Flow(float diameter, float nozzle_diameter) : Flow(diameter, diameter, diameter, nozzle_diameter, true) {}

    // Vertical spacing between extrusion layers.
    float height() const { return m_height; }
    // Width of the model rounded-rectangle, narrow-elipse, or
    // circular-bridge flow cross-section.
    float width() const { return m_width; }
    coord_t scaled_width() const { return coord_t(scale_(m_width)); }
    // Horizontal spacing between the extrusion centerlines.
    float spacing() const { return m_spacing; }
    coord_t scaled_spacing() const { return coord_t(scale_(m_spacing)); }
    // Nozzle diameter.
    float nozzle_diameter() const { return m_nozzle_diameter; }
    // Is it a bridge?
    bool bridge() const { return m_bridge; }
    // Is this a narrow infill flow?
    bool narrow() const { return m_width < m_height; }

    // The diameter of an equivalent volume circular flow.
    float diameter() const { return m_bridge ? m_width : area2dmr(mm3_per_mm()); }
    // The width/spacing ratio.
    float density() const { return m_width / m_spacing; }
    // The horizontal overlap-factor.
    float overlap_factor() const { return m_bridge ? bridge_overlap(m_width, m_spacing) : rrect_overlap(m_width, m_spacing, m_height); }
    // The vertical overlap-factor. For normal lines we return 1 to match 100% flow ratio.
    float vertical_overlap_factor() const { return m_bridge ? bridge_vertoverlap(m_width, m_height) : 1.0; }
    // The fill ratio of the line cross-section area to available height*spacing area.
    float flow_ratio() const { return m_bridge ? bridge_flow(m_width, m_spacing, m_height) : rrect_flow(m_width, m_spacing, m_height); }
    // Cross section area of the extrusion.
    double mm3_per_mm() const { return m_bridge ? bridge_area(m_width) : rrect_area(m_width, m_height); }

    // Override and set spacing, width, or height without changing the other attributes.
    void set_height(float height)
    {
        assert(height > 0);
        m_height = height;
    }
    void set_width(float width)
    {
        assert(width > 0);
        m_width = width;
    }
    void set_spacing(float spacing)
    {
        assert(spacing > 0);
        m_spacing = spacing;
    }

    // Elephant foot compensation spacing to be used to detect narrow parts, where the elephant foot compensation cannot be applied.
    // To be used on frExternalPerimeter only.
    // Enable some perimeter squish (see INSET_OVERLAP_TOLERANCE).
    // Here an overlap of 0.2x external perimeter spacing is allowed for by the elephant foot compensation.
    coord_t scaled_elephant_foot_spacing() const { return coord_t(0.5f * float(this->scaled_width() + 0.6f * this->scaled_spacing())); }

    inline bool operator==(const Flow& rhs) const
    { return m_width == rhs.m_width && m_height == rhs.m_height && m_nozzle_diameter == rhs.m_nozzle_diameter && m_bridge == rhs.m_bridge; }

    inline bool operator!=(const Flow& rhs) const { return !operator==(rhs); }

    bool operator<(const Flow& rhs) const { return this->mm3_per_mm() < rhs.mm3_per_mm(); }

    // Create a modified flow with a different width while maintaining the
    // same overlap-ratios. For bridges this scales spacing and height. For
    // normal lines it preserves height and adjusts spacing.
    Flow with_width(float width) const
    {
        float scale = width / m_width;
        return m_bridge ? Flow(width, scale * m_height, scale * m_spacing, m_nozzle_diameter, m_bridge) :
                          Flow(width, m_height, rrect_spacing(width, m_height, overlap_factor()), m_nozzle_diameter, m_bridge);
    }

    // Create a modified flow with a different height while maintaining
    // overlap-factors. For bridge lines this scales width and spacing. For
    // normal lines it preserves spacing and adjusts width.
    Flow with_height(float height) const
    {
        float scale = height / m_height;
        return m_bridge ? Flow(scale * m_width, height, scale * m_spacing, m_nozzle_diameter, m_bridge) :
                          Flow(rrect_width(m_spacing, height, overlap_factor()), height, m_spacing, m_nozzle_diameter, m_bridge);
    }

    // Create a modified flow with a different spacing while maintaining
    // overlap-factors. For bridge lines this scales width and height. For
    // normal lines it preserves height and adjusts width.
    Flow with_spacing(float spacing) const
    {
        float scale = spacing / m_spacing;
        return m_bridge ? Flow(scale * m_width, scale * m_height, spacing, m_nozzle_diameter, m_bridge) :
                          Flow(rrect_width(spacing, m_height, overlap_factor()), m_height, spacing, m_nozzle_diameter, m_bridge);
    }

    // Create a modified flow by scaling the flow area by a ratio while
    // maintaining overlap-factors. This is the same as scaling width by
    // the appropriate amount.
    Flow with_flow_ratio(double ratio) const
    {
        float width = m_bridge ?
                          // For bridges, scale width by sqrt(ratio);
                          std::sqrt(ratio) * m_width :
                          // For normal lines, get the width for scaled spacing using
                          // overlap-factor=1 (the equivalent rectangular area width).
                          rrect_width(ratio * rrect_spacing(m_width, m_height), m_height);
        return this->with_width(width);
    }

    static Flow bridging_flow(float dmr, float nozzle_diameter) { return Flow(dmr, nozzle_diameter); }

    static Flow new_from_config_width(FlowRole role, const ConfigOptionFloatOrPercent& width, float nozzle_diameter, float height);

    // Spacing of normal layer lines using the rounded-rectangle or narrow-elipse model.
    static constexpr float rounded_rectangle_extrusion_spacing(float width, float height) { return rrect_spacing(width, height); }
    // Width of normal layer lines using the rounded-rectangle or narrow-elipse model.
    static constexpr float rounded_rectangle_extrusion_width_from_spacing(float spacing, float height)
    { return rrect_width(spacing, height); }
    // Spacing of round thread extrusions.
    static constexpr float bridge_extrusion_spacing(float dmr) { return dmr; }

    // Sane extrusion width defautl based on nozzle diameter.
    // The defaults were derived from manual Prusa MK3 profiles.
    static float auto_extrusion_width(FlowRole role, float nozzle_diameter);

    // Extrusion width from full config, taking into account the defaults (when set to zero) and ratios (percentages).
    // Precise value depends on layer index (1st layer vs. other layers vs. variable layer height),
    // on active extruder etc. Therefore the value calculated by this function shall be used as a hint only.
    static double extrusion_width(const std::string& opt_key,
                                  const ConfigOptionFloatOrPercent* opt,
                                  const ConfigOptionResolver& config,
                                  const unsigned int first_printing_extruder = 0);
    static double extrusion_width(const std::string& opt_key,
                                  const ConfigOptionResolver& config,
                                  const unsigned int first_printing_extruder = 0);

private:
    Flow(float width, float height, float spacing, float nozzle_diameter, bool bridge)
        : m_width(width), m_height(height), m_spacing(spacing), m_nozzle_diameter(nozzle_diameter), m_bridge(bridge)
    {}

    float m_width{0};
    float m_height{0};
    float m_spacing{0};
    float m_nozzle_diameter{0};
    bool m_bridge{false};
};

extern Flow support_material_flow(const PrintObject* object, float layer_height = 0.f);
extern Flow support_transition_flow(const PrintObject* object); // BBS
extern Flow support_material_1st_layer_flow(const PrintObject* object, float layer_height = 0.f);
extern Flow support_material_interface_flow(const PrintObject* object, float layer_height = 0.f);

} // namespace Slic3r

#endif
