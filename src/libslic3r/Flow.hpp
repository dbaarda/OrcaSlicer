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
// Note overlap-factor=1 for a rounded-rectangle line is when the adjacent
// lines overlap enough that the cross-section area is "squished" to
// completely fill the available rectangular height*spacing area, for a 100%
// flow_ratio.
//
// For circular bridge lines, they can be "squished" horizontally and
// vertically into a rectangulra vertical-spacing * horizontal-spacing
// available area. We define the horizontal and vertical overlap-factor=1 as
// when the spacing area is a square with the same area as the circular
// cross-section.
//
// Area of a circle, area = diameter^2 * PI_4.
static const double PI_4 = PI/4;
// Length of a square's side with the same area as a circle, length = diameter * SQRTPI_2
static const double SQRTPI_2 = std::sqrt(PI_4);
// Overlap distance for overlap-factor=1 of a rounded-rectangle line, overlap = height * RRLINE_OVERLAP
static const double RRLINE_OVERLAP = 1. - PI_4;  // area ratio of the difference between a square and its enclosed circle.
// Overlap distance for overlap-factor=1 of a circular bridge line, overlap = diameter * BRIDGE_OVERLAP
static const double BRIDGE_OVERLAP = 1. - SQRTPI_2;     // length ratio of difference between a circle's diameter and same-area square's side.

enum FlowRole {
    frExternalPerimeter,
    frPerimeter,
    frInfill,
    frSolidInfill,
    frTopSolidInfill,
    frExternalBridge,
    frInternalBridge,
    frSupportMaterial,
    frSupportMaterialInterface,
    frSupportTransition,  // BBS
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
    // For normal non-bridging flows using a rounded-rectangle extrusion
    // cross-section model;
    //
    // * m_width is the full width of the rounded-rectangle cross-section.
    // * m_height is the vertical height and distance between layers of the
    //   rounded-rectangle cross-section.
    // * m_spacing is the horizontal distance between adjacent and possibly
    //   overlapping lines.
    //
    // For bridging flows using a circular extrusion cross-section model;
    //
    // * m_width is the diameter of the circular cross section.
    // * m_height is the vertical distance between possibly overlapping
    //   layers of bridge lines.
    // * m_spacing is the horizontal distance between possibly overlapping
    //   layers of bridge lines.
    //
    // Although currently bridge lines are assumed to never vertically
    // overlap so m_width = m_height = diameter, always use m_width or
    // diameter() for the diameter in case support for vertical overlapping
    // is added in the future.
    Flow() = default;
    // Initialize a standard rounded-rectangle line flow.
    Flow(float width, float height, float nozzle_diameter) :
        Flow(width, height, rounded_rectangle_extrusion_spacing(width, height), nozzle_diameter, false) {}
    // Initialize a standard circular-bridge line flow.
    Flow(float diameter, float nozzle_diameter) :
        Flow(diameter, diameter, diameter, nozzle_diameter, true) {}

    // Vertical spacing between extrusion layers.
    float   height()          const { return m_height; }
    // Width of the model rounded-rectangle or circular-bridge flow cross-section.
    float   width()           const { return m_width; }
    coord_t scaled_width()    const { return coord_t(scale_(m_width)); }
    // Horizontal spacing between the extrusion centerlines.
    float   spacing()         const { return m_spacing; }
    coord_t scaled_spacing()  const { return coord_t(scale_(m_spacing)); }
    // Nozzle diameter.
    float   nozzle_diameter() const { return m_nozzle_diameter; }
    // Is it a bridge?
    bool    bridge()          const { return m_bridge; }

    // The diameter of an equivalent volume circular flow.
    float   diameter() const { return m_bridge ? m_width : std::sqrt(mm3_per_mm()/PI_4); }
    // The width/spacing ratio.
    float   density()         const { return m_width / m_spacing; }
    // The horizontal overlap-factor for both rounded-rectangle and circular-bridge lines.
    float   overlap_factor()  const { return (m_width - m_spacing) / (m_bridge ? m_width * BRIDGE_OVERLAP : m_height * RRLINE_OVERLAP); }
    // The vertical overlap-factor for circular-bridge lines. For rounded-rectangles we return 1 to match 100% flow ratio.
    float   vertical_overlap_factor() const { return (m_bridge ? (m_width - m_height)/(m_width * BRIDGE_OVERLAP) : 1.); }
    // The fill ratio of the line cross-section area to available height*spacing area.
    float   flow_ratio() const { return mm3_per_mm() / (m_spacing * m_height); }
    // Cross section area of the extrusion.
    double  mm3_per_mm()      const;

    // Override and set spacing, width, or height without changing the other attributes.
    void    set_height(float height) { assert(height <= m_width); m_height = height; }
    void    set_width(float width) { assert (m_height <= width); m_width = width; }
    void    set_spacing(float spacing) { assert(spacing > 0); m_spacing = spacing; }

    // Elephant foot compensation spacing to be used to detect narrow parts, where the elephant foot compensation cannot be applied.
    // To be used on frExternalPerimeter only.
    // Enable some perimeter squish (see INSET_OVERLAP_TOLERANCE).
    // Here an overlap of 0.2x external perimeter spacing is allowed for by the elephant foot compensation.
    coord_t scaled_elephant_foot_spacing() const { return coord_t(0.5f * float(this->scaled_width() + 0.6f * this->scaled_spacing())); }

    inline bool operator==(const Flow &rhs) const { return m_width == rhs.m_width && m_height == rhs.m_height && m_nozzle_diameter == rhs.m_nozzle_diameter && m_bridge == rhs.m_bridge; }

    inline bool operator!=(const Flow &rhs) const { return ! operator==(rhs); }

    bool operator <(const Flow &rhs) const {
        return this->mm3_per_mm() < rhs.mm3_per_mm();
    }

    // Create modified flows changing width, spacing, or height while
    // maintaining the same overlap-ratio. Note for rounded-rectangle lines
    // maintaining overlap-factor scales the overlap or gap distance linearly
    // with height, so the density (width/spacing) changes, but for circular
    // bridge lines maintaining overlap factor maintains the same density, so
    // the gap/overlap distance changes.
    Flow        with_width (float width)  const {
        assert(! m_bridge);
        return Flow(width, m_height, rounded_rectangle_extrusion_spacing(width, m_height), m_nozzle_diameter, m_bridge);
    }
    Flow        with_height(float height) const {
        if (m_bridge) {
          // maintaining the vertical and horizontal overlap-factor requires
          // scaling height, width, and spacing together.
          auto scale = height / m_height;
          return Flow(m_width*scale, height, m_spacing*scale, m_nozzle_diameter, m_bridge);
        } else {
          return Flow(m_width, height, m_width - overlap_factor()*RRLINE_OVERLAP*height, m_nozzle_diameter, m_bridge);
        }
    }
    Flow        with_spacing(float spacing) const;
    // Adjust the width / height of a rounded extrusion model to reach the prescribed cross section area while maintaining extrusion spacing.
    Flow        with_cross_section(float area) const;
    Flow        with_flow_ratio(double ratio) const { return this->with_cross_section(this->mm3_per_mm() * ratio); }

    static Flow bridging_flow(float dmr, float nozzle_diameter) { return Flow(dmr, nozzle_diameter); }

    static Flow new_from_config_width(FlowRole role, const ConfigOptionFloatOrPercent &width, float nozzle_diameter, float height);

    // Spacing of extrusions with rounded extrusion model.
    static constexpr float rounded_rectangle_extrusion_spacing(float width, float height) { return width - height * RRLINE_OVERLAP; }
    // Width of extrusions with rounded extrusion model.
    static constexpr float rounded_rectangle_extrusion_width_from_spacing(float spacing, float height) { return spacing + height * RRLINE_OVERLAP; }
    // Spacing of round thread extrusions.
    static constexpr float bridge_extrusion_spacing(float dmr) { return dmr; }

    // Sane extrusion width defautl based on nozzle diameter.
    // The defaults were derived from manual Prusa MK3 profiles.
    static float auto_extrusion_width(FlowRole role, float nozzle_diameter);

    // Extrusion width from full config, taking into account the defaults (when set to zero) and ratios (percentages).
    // Precise value depends on layer index (1st layer vs. other layers vs. variable layer height),
    // on active extruder etc. Therefore the value calculated by this function shall be used as a hint only.
    static double extrusion_width(const std::string &opt_key, const ConfigOptionFloatOrPercent *opt, const ConfigOptionResolver &config, const unsigned int first_printing_extruder = 0);
    static double extrusion_width(const std::string &opt_key, const ConfigOptionResolver &config, const unsigned int first_printing_extruder = 0);

private:
    Flow(float width, float height, float spacing, float nozzle_diameter, bool bridge) :
        m_width(width), m_height(height), m_spacing(spacing), m_nozzle_diameter(nozzle_diameter), m_bridge(bridge)
        {
            // Gap fill violates this condition.
            //assert(width >= height);
        }

    float       m_width { 0 };
    float       m_height { 0 };
    float       m_spacing { 0 };
    float       m_nozzle_diameter { 0 };
    bool        m_bridge { false };
};

extern Flow support_material_flow(const PrintObject* object, float layer_height = 0.f);
extern Flow support_transition_flow(const PrintObject *object); //BBS
extern Flow support_material_1st_layer_flow(const PrintObject *object, float layer_height = 0.f);
extern Flow support_material_interface_flow(const PrintObject *object, float layer_height = 0.f);

}

#endif
