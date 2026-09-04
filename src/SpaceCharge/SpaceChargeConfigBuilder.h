/** @file SpaceChargeConfigBuilder.h @brief Converts parser state into runtime configuration. */

#ifndef OPALX_SPACE_CHARGE_CONFIG_BUILDER_H
#define OPALX_SPACE_CHARGE_CONFIG_BUILDER_H

#include "SpaceCharge/SpaceChargeConfig.h"

#include <vector>

class EmissionSource;
class FieldSolverCmd;

namespace opalx::spacecharge {

    [[nodiscard]] SpaceChargeConfig buildSpaceChargeConfig(
            const FieldSolverCmd& fieldSolver,
            const std::vector<std::vector<EmissionSource*>>& emissionSources);

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_CONFIG_BUILDER_H
