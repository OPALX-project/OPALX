/**
 * @file BinningPlan.cpp
 * @brief Explicit instantiation of production Cartesian PIC bin iteration.
 */

#include "SpaceCharge/Pic3D/BinningPlan.h"

namespace opalx::spacecharge {

    template class BinningPlan<double, 3>;

}  // namespace opalx::spacecharge
