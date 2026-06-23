#include <catch2/catch_all.hpp>

#include <numeric>
#include <sstream>

#include "test_data.hpp" // get access to init_print, etc

#include "libslic3r/Config.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/Config.hpp"
#include "libslic3r/GCodeReader.hpp"
#include "libslic3r/Flow.hpp"
#include "libslic3r/libslic3r.h"

using namespace Slic3r::Test;
using namespace Slic3r;

/// Test the expected behavior for auto-width, spacing, etc
SCENARIO("Flow:new_from_config_width() for non-bridges", "[Flow]") {
    GIVEN("Nozzle Diameter of 0.4, a desired width of 1mm and layer height of 0.5") {
        ConfigOptionFloatOrPercent width(1.0, false);
        float nozzle_diameter = 0.4f;
        float layer_height    = 0.4f;

        // Spacing for non-bridges is has some overlap
        THEN("External perimeter flow has spacing fixed to 1.125 * nozzle_diameter") {
            auto flow = Flow::new_from_config_width(frExternalPerimeter, ConfigOptionFloatOrPercent(0, false), nozzle_diameter,
                                                    layer_height);
            REQUIRE(flow.spacing() == Catch::Approx(1.125 * nozzle_diameter - layer_height * (1.0 - M_PI_4)));
        }

        THEN("Internal perimeter flow has spacing fixed to 1.125 * nozzle_diameter") {
            auto flow = Flow::new_from_config_width(frPerimeter, ConfigOptionFloatOrPercent(0, false), nozzle_diameter, layer_height);
            REQUIRE(flow.spacing() == Catch::Approx(1.125 * nozzle_diameter - layer_height * (1.0 - M_PI_4)));
        }
        THEN("Spacing for supplied width is 0.8927f") {
            auto flow = Flow::new_from_config_width(frExternalPerimeter, width, nozzle_diameter, layer_height);
            REQUIRE(flow.spacing() == Catch::Approx(width.value - layer_height * (1.0 - M_PI_4)));
            flow = Flow::new_from_config_width(frPerimeter, width, nozzle_diameter, layer_height);
            REQUIRE(flow.spacing() == Catch::Approx(width.value - layer_height * (1.0 - M_PI_4)));
        }
    }
    /// Check the min/max
    GIVEN("Nozzle Diameter of 0.25") {
        float nozzle_diameter = 0.25f;
        float layer_height    = 0.5f;
        WHEN("layer height is set to 0.2") {
            layer_height = 0.15f;
            THEN("Max width is set.") {
                auto flow = Flow::new_from_config_width(frPerimeter, ConfigOptionFloatOrPercent(0, false), nozzle_diameter, layer_height);
                REQUIRE(flow.width() == Catch::Approx(1.125 * nozzle_diameter));
            }
        }
        WHEN("Layer height is set to 0.25") {
            layer_height = 0.25f;
            THEN("Min width is set.") {
                auto flow = Flow::new_from_config_width(frPerimeter, ConfigOptionFloatOrPercent(0, false), nozzle_diameter, layer_height);
                REQUIRE(flow.width() == Catch::Approx(1.125 * nozzle_diameter));
            }
        }
    }
}

SCENARIO("Flow: Flow math for normal-lines", "[Flow]") {
    /// Check for width < height using the narrow-elipse flow model.
    GIVEN("Input width=0.5 is greater than height=0.2 for nozzle_diameter=0.4") {
        float width           = 0.5f;
        float height          = 0.2f;
        float nozzle_diameter = 0.4f;
        WHEN("flow=Flow(width, height, nozzle_diameter) is called") {
            auto flow = Flow(width, height, nozzle_diameter);
            THEN("flow attributes are correct for the rounded-rectangle extrusion model") {
                CHECK(flow.width() == width);
                CHECK(flow.height() == height);
                CHECK(flow.spacing() == Catch::Approx(width - height * (1.0 - M_PI_4)));
                CHECK(flow.nozzle_diameter() == nozzle_diameter);
                CHECK(flow.bridge() == false);
                CHECK(flow.narrow() == false);
                CHECK(flow.mm3_per_mm() == Catch::Approx(flow.spacing() * height));
                CHECK(flow.diameter() == Catch::Approx(std::sqrt(flow.mm3_per_mm() / M_PI_4)));
                CHECK(flow.density() == Catch::Approx(width / flow.spacing()));
                CHECK(flow.overlap_factor() == Catch::Approx(1.0));
                CHECK(flow.vertical_overlap_factor() == Catch::Approx(1.0));
                CHECK(flow.flow_ratio() == Catch::Approx(1.0));
            }
        }
    }

    /// Check for width == height boundary between switching flow models.
    GIVEN("Input width=0.2 equals height=0.2 for nozzle_diameter=0.4") {
        float width           = 0.2f;
        float height          = 0.2f;
        float nozzle_diameter = 0.4f;
        WHEN("flow=Flow(width, height, nozzle_diameter) is called") {
            auto flow = Flow(width, height, nozzle_diameter);
            THEN("flow attributes are correct for the rounded-rectangle and narrow-elipse extrusion model boundary") {
                CHECK(flow.width() == width);
                CHECK(flow.height() == height);
                CHECK(flow.spacing() == Catch::Approx(width - height * (1.0 - M_PI_4)));
                CHECK(flow.spacing() == Catch::Approx(width * M_PI_4));
                CHECK(flow.nozzle_diameter() == nozzle_diameter);
                CHECK(flow.bridge() == false);
                CHECK(flow.narrow() == false);
                CHECK(flow.mm3_per_mm() == Catch::Approx(width * height * M_PI_4));
                CHECK(flow.diameter() == Catch::Approx(width));
                CHECK(flow.density() == Catch::Approx(1.0 / M_PI_4));
                CHECK(flow.overlap_factor() == Catch::Approx(1.0));
                CHECK(flow.vertical_overlap_factor() == Catch::Approx(1.0));
                CHECK(flow.flow_ratio() == Catch::Approx(1.0));
            }
        }
    }

    /// Check for width < height using the narrow-elipse flow model.
    GIVEN("Input width=0.1 is less than height=0.2 for nozzle_diameter=0.4") {
        float width           = 0.1f;
        float height          = 0.2f;
        float nozzle_diameter = 0.4f;
        WHEN("flow=Flow(width, height, nozzle_diameter) is called") {
            auto flow = Flow(width, height, nozzle_diameter);
            THEN("flow attributes are correct for the narrow-elipse extrusion model") {
                CHECK(flow.width() == width);
                CHECK(flow.height() == height);
                CHECK(flow.spacing() == Catch::Approx(width * M_PI_4));
                CHECK(flow.nozzle_diameter() == nozzle_diameter);
                CHECK(flow.bridge() == false);
                CHECK(flow.narrow() == true);
                CHECK(flow.mm3_per_mm() == Catch::Approx(width * height * M_PI_4));
                CHECK(flow.diameter() == Catch::Approx(std::sqrt(width * height)));
                CHECK(flow.density() == Catch::Approx(1.0 / M_PI_4));
                CHECK(flow.overlap_factor() == Catch::Approx(1.0));
                CHECK(flow.vertical_overlap_factor() == Catch::Approx(1.0));
                CHECK(flow.flow_ratio() == Catch::Approx(1.0));
            }
        }
    }
}

SCENARIO("Flow: Flow math for bridge-lines", "[Flow]") {
    GIVEN("Input diameter=0.1 for nozzle_diameter=0.4") {
        float nozzle_diameter = 0.4f;
        float diameter        = 0.1f;
        WHEN("flow=Flow(diameter, nozzle_diameter) is called") {
            auto flow = Flow(diameter, nozzle_diameter);
            THEN("flow attributes are correct for the circular extrusion model") {
                CHECK(flow.width() == diameter);
                CHECK(flow.height() == diameter);
                CHECK(flow.spacing() == diameter);
                CHECK(flow.nozzle_diameter() == nozzle_diameter);
                CHECK(flow.bridge() == true);
                CHECK(flow.narrow() == false);
                CHECK(flow.mm3_per_mm() == Catch::Approx(diameter * diameter * M_PI_4));
                CHECK(flow.diameter() == diameter);
                CHECK(flow.density() == Catch::Approx(1.0));
                CHECK(flow.overlap_factor() == Catch::Approx(0.0));
                CHECK(flow.vertical_overlap_factor() == Catch::Approx(0.0));
                CHECK(flow.flow_ratio() == Catch::Approx(M_PI_4));
            }
        }
    }
}

SCENARIO("Flow: argument validity checking and handling", "[Flow]") {
    GIVEN("Input nozzle_diameter=0.4") {
        float nozzle_diameter = 0.4f;
        WHEN("diameter=0.1 < nozzle_diameter") {
            float diameter = 0.1f;
            THEN("Initializing bridge-line flows does not throw") {
                CHECK_NOTHROW(Flow(diameter, nozzle_diameter));
            }
        }
        WHEN("diameter=0.4 == nozzle_diameter") {
            float diameter = 0.4f;
            THEN("Initializing bridge-line flows does not throw") {
                CHECK_NOTHROW(Flow(diameter, nozzle_diameter));
            }
        }
        WHEN("diameter=0.5 > nozzle_diameter") {
            float diameter = 0.5f;
            THEN("Initializing bridge-line flows throws FlowErrorHeightTooLarge") {
                CHECK_THROWS_AS(Flow(diameter, nozzle_diameter), FlowErrorHeightTooLarge);
            }
        }
        WHEN("diameter=-0.1 is negative") {
            float diameter = -0.1f;
            THEN("Initializing bridge-line flows throws FlowErrorNegativeHeight") {
                CHECK_THROWS_AS(Flow(diameter, nozzle_diameter), FlowErrorNegativeHeight);
            }
        }
        WHEN("height=0.2 < nozzle_diameter") {
            float height = 0.2f;
            THEN("Initializing normal-line flows does not throw") {
                CHECK_NOTHROW(Flow(0.1, height, nozzle_diameter));
                CHECK_NOTHROW(Flow(0.4, height, nozzle_diameter));
                CHECK_NOTHROW(Flow(0.5, height, nozzle_diameter));
            }
        }
        WHEN("height=0.4 == nozzle_diameter") {
            float height = 0.4f;
            THEN("Initializing normal-line flows does not throw") {
                CHECK_NOTHROW(Flow(0.1, height, nozzle_diameter));
                CHECK_NOTHROW(Flow(0.4, height, nozzle_diameter));
                CHECK_NOTHROW(Flow(0.5, height, nozzle_diameter));
            }
        }
        WHEN("height=0.5 > nozzle_diameter") {
            float height = 0.5f;
            THEN("Initializing normal-line flows throws FlowErrorHeightTooLarge") {
                CHECK_THROWS_AS(Flow(0.1, height, nozzle_diameter), FlowErrorHeightTooLarge);
                CHECK_THROWS_AS(Flow(0.4, height, nozzle_diameter), FlowErrorHeightTooLarge);
                CHECK_THROWS_AS(Flow(0.5, height, nozzle_diameter), FlowErrorHeightTooLarge);
            }
        }
        WHEN("height=-0.1 is negative") {
            float height = -0.1f;
            THEN("Initializing normal-line flows throws FlowErrorNegativeHeight") {
                CHECK_THROWS_AS(Flow(0.1, height, nozzle_diameter), FlowErrorNegativeHeight);
                CHECK_THROWS_AS(Flow(0.4, height, nozzle_diameter), FlowErrorNegativeHeight);
                CHECK_THROWS_AS(Flow(0.5, height, nozzle_diameter), FlowErrorNegativeHeight);
            }
        }
    }
    GIVEN("Input nozzle_diameter=0 for unknown") {
        float nozzle_diameter = 0.0f;
        WHEN("diameter=0.5 > nozzle_diameter") {
            float diameter = 0.5f;
            THEN("Initializing bridge-line flows does not throw") {
                CHECK_NOTHROW(Flow(diameter, nozzle_diameter));
            }
        }
        WHEN("height=0.5 > nozzle_diameter") {
            float height = 0.5f;
            THEN("Initializing normal-line flows does not throw") {
                CHECK_NOTHROW(Flow(0.1, height, nozzle_diameter));
                CHECK_NOTHROW(Flow(0.4, height, nozzle_diameter));
                CHECK_NOTHROW(Flow(0.5, height, nozzle_diameter));
            }
        }
    }
}
