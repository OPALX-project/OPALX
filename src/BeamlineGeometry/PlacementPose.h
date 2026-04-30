#ifndef OPALX_PlacementPose_HH
#define OPALX_PlacementPose_HH

#include "Algorithms/CoordinateSystemTrafo.h"

/**
 * @class PlacementPose
 * @brief Nominal rigid placement transform.
 *
 * In the language of the placement note, this stores the nominal rigid
 * placement transform \f$T_i\f$ that embeds the canonical local chart of an
 * element into its parent frame before any local correction is applied.
 */
class PlacementPose {
public:
    PlacementPose() = default;
    explicit PlacementPose(const CoordinateSystemTrafo& parentToNominal)
        : parentToNominal_m(parentToNominal) {}

    const CoordinateSystemTrafo& getParentToNominal() const { return parentToNominal_m; }

private:
    CoordinateSystemTrafo parentToNominal_m;
};

#endif
