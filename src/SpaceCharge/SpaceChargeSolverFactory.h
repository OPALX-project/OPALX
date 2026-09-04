/**
 * @file SpaceChargeSolverFactory.h
 * @brief Constructs the space-charge solver selected by immutable setup configuration.
 */

#ifndef OPALX_SPACE_CHARGE_SOLVER_FACTORY_H
#define OPALX_SPACE_CHARGE_SOLVER_FACTORY_H

#include "PartBunch/PartBunchFwd.h"
#include "SpaceCharge/SpaceChargeConfig.h"
#include "SpaceCharge/SpaceChargeSolver.h"

#include <memory>

class DataSink;

namespace opalx::spacecharge {

    /** @brief Construction-time selection point for concrete space-charge algorithms. */
    class SpaceChargeSolverFactory {
    public:
        [[nodiscard]] static std::unique_ptr<SpaceChargeSolver> create(
                SpaceChargeConfig config, PartBunch_t& bunch, DataSink* dataSink);
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_SOLVER_FACTORY_H
