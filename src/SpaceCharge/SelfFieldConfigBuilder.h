/**
 * @file SelfFieldConfigBuilder.h
 * @brief Converts setup-time parser objects into immutable self-field configuration.
 */

#ifndef OPALX_SPACE_CHARGE_SELF_FIELD_CONFIG_BUILDER_H
#define OPALX_SPACE_CHARGE_SELF_FIELD_CONFIG_BUILDER_H

#include "SpaceCharge/SelfFieldConfig.h"

#include <vector>

class EmissionSource;
class FieldSolverCmd;

namespace opalx::spacecharge {

    /** @brief One-time setup adapter between parser objects and the self-field subsystem. */
    class SelfFieldConfigBuilder {
    public:
        /**
         * @brief Snapshot the field-solver, binning, option, and emission-source configuration.
         *
         * No parser pointer is retained in the returned value.
         */
        [[nodiscard]] static SelfFieldConfig build(
                const FieldSolverCmd& fieldSolver,
                const std::vector<std::vector<EmissionSource*>>& emissionSources);
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_SELF_FIELD_CONFIG_BUILDER_H
