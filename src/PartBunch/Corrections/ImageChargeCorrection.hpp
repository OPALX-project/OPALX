#ifndef OPALX_IMAGE_CHARGE_CORRECTION_HPP
#define OPALX_IMAGE_CHARGE_CORRECTION_HPP

#include <cstddef>
#include <optional>

#include "PartBunch/Corrections/BoundaryCorrection.hpp"

/**
 * @brief Explicit image-charge Dirichlet correction policy.
 *
 * The policy describes the second solve pass for the image distribution. Execution stays in the
 * binned solver: scatter only image particles, solve with the standard backend request, accumulate
 * with a negative magnetic-field sign, and allow Dirichlet-plane diagnostics.
 */
template <typename T, unsigned Dim>
class ImageChargeCorrection final : public BoundaryCorrection<T, Dim> {
public:
    /**
     * @brief Enable or disable the correction and store the physical z plane.
     */
    void configure(bool enabled, double zPlane) {
        enabled_ = enabled;
        zPlane_  = zPlane;
    }

    const char* name() const override { return "image-charge"; }
    BoundaryCorrectionKind kind() const override { return BoundaryCorrectionKind::ImageCharge; }
    bool isEnabled() const override { return enabled_; }
    double planeZ() const override { return zPlane_; }

    std::optional<BoundaryCorrectionPass<T, Dim>> makePass(
            const BinContext<T, Dim>&, const Vector_t<double, Dim>&,
            std::size_t) const override {
        if (!enabled_) {
            return std::nullopt;
        }

        BoundaryCorrectionPass<T, Dim> pass;
        pass.scatterMode = ImageScatterMode::ImageOnly;
        pass.bFieldSign = -1.0;
        pass.flipAxis = -1;
        pass.dumpDirichletPlane = true;
        pass.forceSkipFieldDump = true;
        return pass;
    }

private:
    bool enabled_ = false;
    double zPlane_ = 0.0;
};

#endif  // OPALX_IMAGE_CHARGE_CORRECTION_HPP
