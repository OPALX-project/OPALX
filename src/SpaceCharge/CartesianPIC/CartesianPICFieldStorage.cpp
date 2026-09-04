/**
 * @file CartesianPICFieldStorage.cpp
 * @brief Explicit instantiation of production Cartesian PIC field storage.
 */

#include "SpaceCharge/CartesianPIC/CartesianPICFieldStorage.h"

namespace opalx::spacecharge {

    template class CartesianPICFieldStorage<double, 3>;

}  // namespace opalx::spacecharge
