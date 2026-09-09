/**
 * @file CartesianPIC3DFieldStorage.cpp
 * @brief Explicit instantiation of production CartesianPIC3D field storage.
 */

#include "SpaceCharge/CartesianPIC3D/CartesianPIC3DFieldStorage.h"

namespace opalx::spacecharge {

    template class CartesianPIC3DFieldStorage<double, 3>;

}  // namespace opalx::spacecharge
