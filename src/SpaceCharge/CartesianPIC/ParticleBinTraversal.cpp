/**
 * @file ParticleBinTraversal.cpp
 * @brief Explicit instantiation of production Cartesian PIC bin iteration.
 */

#include "SpaceCharge/CartesianPIC/ParticleBinTraversal.h"

namespace opalx::spacecharge {

    template class ParticleBinTraversal<double, 3>;

}  // namespace opalx::spacecharge
