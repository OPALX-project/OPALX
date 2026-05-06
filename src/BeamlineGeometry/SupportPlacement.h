#ifndef OPALX_SupportPlacement_HH
#define OPALX_SupportPlacement_HH

#include "Algorithms/CoordinateSystemTrafo.h"

/**
 * @class SupportPlacement
 * @brief Optional offset from the canonical body frame to the support frame.
 *
 * This is the bridge to the document's distinction between the canonical local
 * chart and the support region \f$\Omega_i\f$ on which field, aperture, or
 * material data are defined. The first redesign stage keeps this as a rigid
 * body-to-support transform.
 */
class SupportPlacement {
public:
    SupportPlacement() = default;
    explicit SupportPlacement(const CoordinateSystemTrafo& bodyToSupport)
        : bodyToSupport_m(bodyToSupport) {}

    const CoordinateSystemTrafo& getBodyToSupport() const { return bodyToSupport_m; }

private:
    CoordinateSystemTrafo bodyToSupport_m;
};

#endif
