/**
 * @file PicWorkspace.cpp
 * @brief Explicit instantiation of the production Cartesian PIC workspace.
 */

#include "SpaceCharge/Pic/PicWorkspace.h"

namespace opalx::spacecharge {

    template class PicWorkspace<double, 3>;

}  // namespace opalx::spacecharge
