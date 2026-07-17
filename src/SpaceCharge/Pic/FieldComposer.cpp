/**
 * @file FieldComposer.cpp
 * @brief Explicit instantiation of the production Cartesian PIC field composer.
 */

#include "SpaceCharge/Pic/FieldComposer.h"

namespace opalx::spacecharge {

    template class FieldComposer<double, 3>;

}  // namespace opalx::spacecharge
