#ifndef OPALX_IMAGE_CHARGE_SCATTER_CONTROLLER_H
#define OPALX_IMAGE_CHARGE_SCATTER_CONTROLLER_H

/**
 * @brief Temporary compatibility holder for image-charge enablement and mirror geometry.
 *
 * PicScatterGather owns deposition and temporary particle transformations. This value holder
 * remains until correction configuration moves completely into CorrectionPlan.
 */
template <typename T, unsigned Dim>
class ImageChargeScatterController {
    static_assert(Dim == 3, "ImageChargeScatterController currently supports Dim == 3 only.");

public:
    ImageChargeScatterController() = default;

    ImageChargeScatterController(bool enabled, double zPlane)
        : enabled_m(enabled), zPlane_m(zPlane) {}

    void configure(bool enabled, double zPlane) {
        enabled_m = enabled;
        zPlane_m  = zPlane;
    }

    [[nodiscard]] bool isEnabled() const { return enabled_m; }
    [[nodiscard]] double getZPlane() const { return zPlane_m; }

private:
    bool enabled_m  = false;
    double zPlane_m = 0.0;
};

#endif  // OPALX_IMAGE_CHARGE_SCATTER_CONTROLLER_H
