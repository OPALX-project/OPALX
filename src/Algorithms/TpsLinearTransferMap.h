// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#ifndef OPAL_TPS_LINEAR_TRANSFER_MAP_H
#define OPAL_TPS_LINEAR_TRANSFER_MAP_H

#include <string>
#include "Structure/LinearTransferMap.h"
#include "tps/Map.hpp"

namespace opalx::maps {
    /**
     * @brief Represent a segment's differential matrix as a first-order TPS map.
     *
     * Local variables are perturbations in (x,x',y,y',zeta,delta), with metres in
     * entries 0,2,4 and dimensionless slopes/delta in 1,3,5. Chart scales are one;
     * the polynomial is F_i(xi)=sum_j M_ij xi_j, with no affine offset. This preserves
     * the existing differential contract, not a claim of exact closed-orbit closure.
     *
     * Boundary frames store lab-to-local axes by row and compensated reference
     * position/time rounded to double. Reference slopes are retained even when the
     * return momentum is not parallel to the fixed section normal. The original
     * LinearTransferMapReference remains authoritative for compensated low parts.
     * Adjacent charts must match exactly for tps::compose. A nonempty occurrenceId
     * identifies one unique traversal within the composition, not an element owner.
     *
     * Requires an active Kokkos runtime and serialized host TPS algebra. No tracking,
     * field evaluation, finite differencing or reference mutation is performed.
     * @throws OpalException Invalid matrix, reference metadata or empty occurrence ID.
     */
    tps::Map<double> toTpsLinearTransferMap(
            const LinearTransferMapDraft&, const std::string& occurrenceId);
}  // namespace opalx::maps
#endif
