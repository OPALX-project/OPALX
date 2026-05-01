#ifndef OPALX_Misalignment_HH
#define OPALX_Misalignment_HH

#include "Algorithms/CoordinateSystemTrafo.h"

/**
 * @class Misalignment
 * @brief Local nominal-to-actual correction transform.
 *
 * The transform maps coordinates from the nominal element body frame into the
 * actual body frame after local survey or misalignment corrections.
 */
class Misalignment {
public:
    Misalignment() = default;
    explicit Misalignment(const CoordinateSystemTrafo& nominalToActual)
        : nominalToActual_m(nominalToActual) {}

    const CoordinateSystemTrafo& getNominalToActual() const { return nominalToActual_m; }

private:
    CoordinateSystemTrafo nominalToActual_m;
};

#endif
