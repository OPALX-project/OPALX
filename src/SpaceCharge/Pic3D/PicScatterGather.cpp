/**
 * @file PicScatterGather.cpp
 * @brief Explicit instantiation of production Cartesian PIC scatter/gather.
 */

#include "SpaceCharge/Pic3D/PicScatterGather.h"

namespace opalx::spacecharge {

    template class PicScatterGather<double, 3>;

}  // namespace opalx::spacecharge
