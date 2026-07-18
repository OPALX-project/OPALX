/**
 * @file CorrectionPlan.cpp
 * @brief Explicit instantiation of production Cartesian PIC correction planning.
 */

#include "SpaceCharge/Pic3D/CorrectionPlan.h"

namespace opalx::spacecharge {

    template class CorrectionPlan<double, 3>;

}  // namespace opalx::spacecharge
