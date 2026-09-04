/**
 * @file RelativisticFieldComposer.cpp
 * @brief Explicit instantiation of the production Cartesian PIC field composer.
 */

#include "SpaceCharge/CartesianPIC/RelativisticFieldComposer.h"

namespace opalx::spacecharge {

    template class RelativisticFieldComposer<double, 3>;

}  // namespace opalx::spacecharge
