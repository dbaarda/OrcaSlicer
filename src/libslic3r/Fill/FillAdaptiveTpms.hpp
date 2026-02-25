#ifndef slic3r_FillAdaptiveTpms_hpp_
#define slic3r_FillAdaptiveTpms_hpp_

#include "FillBase.hpp"

namespace Slic3r {

class PrintObject;

namespace FillAdaptiveTpms {

class TpmsDensityField;
// To keep the definition opaque, we have to define a custom deleter.
void TpmsDensityFieldDeleter(TpmsDensityField* p);
using TpmsDensityFieldPtr = std::unique_ptr<TpmsDensityField, TpmsDensityFieldDeleter>;

TpmsDensityFieldPtr build_tpms_density_field(const PrintObject& print_object, const std::function<void()>& throw_on_cancel_callback);

class Filler : public Slic3r::Fill
{
public:
    ~Filler() override = default;
    bool is_self_crossing() override { return false; }

    TpmsDensityField* tpms_density_field{nullptr};

protected:
    Fill* clone() const override { return new Filler(*this); }

    void _fill_surface_single(const FillParams&              params,
                              unsigned int                   thickness_layers,
                              const std::pair<float, Point>& direction,
                              ExPolygon                      expolygon,
                              Polylines&                     polylines_out) override;

    // Let the G-code export reoder the infill lines.
    bool no_sort() const override { return false; }
};

} // namespace FillAdaptiveTpms
} // namespace Slic3r

#endif // slic3r_FillAdaptiveTpms_hpp_
