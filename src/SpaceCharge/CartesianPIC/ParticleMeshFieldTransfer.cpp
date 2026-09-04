/**
 * @file ParticleMeshFieldTransfer.cpp
 * @brief Explicit instantiation of production Cartesian PIC scatter/gather.
 */

#include "SpaceCharge/CartesianPIC/ParticleMeshFieldTransfer.h"

namespace opalx::spacecharge {

    template class ParticleMeshFieldTransfer<double, 3>;

}  // namespace opalx::spacecharge
