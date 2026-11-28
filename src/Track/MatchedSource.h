#ifndef OPAL_MATCHED_SOURCE_H
#define OPAL_MATCHED_SOURCE_H

#include "Structure/EmissionSource.h"
#include "Utilities/OpalException.h"

/// Validate a MOMENTS source before allocation and again before sampling.
/// Counts, charge and polarization remain owned by the existing beam/source path.
inline void validateMatchedSource(const EmissionSource& source) {
    if (source.getT0() != 0 || source.getEmissionModel() != "NONE" || source.getKineticEnergy() != 0
        || source.getZeroFaceR0Z() || source.getShiftedGreensFunction())
        throw OpalException(
                "MOMENTS",
                "Requires source T0=0, EMISSIONMODEL=NONE, EKIN=0 and no emission/image plane.");
    const auto r = source.getR0(), p = source.getP0();
    for (unsigned i = 0; i < 3; ++i)
        if (r[i] != 0 || p[i] != 0)
            throw OpalException("MOMENTS", "Source R0/P0 offsets must be zero.");
}
#endif
