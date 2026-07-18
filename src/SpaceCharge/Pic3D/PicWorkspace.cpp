/**
 * @file PicWorkspace.cpp
 * @brief Explicit instantiation of the production Cartesian PIC workspace.
 */

#include "SpaceCharge/Pic3D/PicWorkspace.h"

namespace opalx::spacecharge {

    template class PicWorkspace<double, 3>;

}  // namespace opalx::spacecharge
