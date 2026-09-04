/**
 * @file SpaceChargeConfigBuilder.h
 * @brief Converts setup-time parser objects into immutable space-charge configuration.
 */

#ifndef OPALX_SPACE_CHARGE_CONFIG_BUILDER_H
#define OPALX_SPACE_CHARGE_CONFIG_BUILDER_H

#include "SpaceCharge/SpaceChargeConfig.h"

#include <vector>

class EmissionSource;
class FieldSolverCmd;

namespace opalx::spacecharge {

    /** @brief One-time conversion boundary from parser objects to space-charge configuration. */
    class SpaceChargeConfigBuilder {
    public:
        /**
         * @brief Snapshot the field-solver, binning, option, and emission-source configuration.
         *
         * No parser pointer is retained in the returned value.
         */
        [[nodiscard]] static SpaceChargeConfig build(
                const FieldSolverCmd& fieldSolver,
                const std::vector<std::vector<EmissionSource*>>& emissionSources);
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_CONFIG_BUILDER_H
