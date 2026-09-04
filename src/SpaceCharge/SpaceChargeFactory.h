/** @file SpaceChargeFactory.h @brief Declares run-lifetime solver construction. */

#ifndef OPALX_SPACE_CHARGE_FACTORY_H
#define OPALX_SPACE_CHARGE_FACTORY_H

#include "PartBunch/PartBunchFwd.h"
#include "SpaceCharge/SpaceChargeConfig.h"
#include "SpaceCharge/SpaceChargeSolver.h"

#include <memory>

class DataSink;

namespace opalx::spacecharge {

    [[nodiscard]] std::unique_ptr<SpaceChargeSolver> makeSpaceChargeSolver(
            SpaceChargeConfig config, PartBunch_t& bunch, DataSink* dataSink);

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_FACTORY_H
