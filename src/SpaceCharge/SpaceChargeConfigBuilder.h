/** @file SpaceChargeConfigBuilder.h @brief Converts the RUN description into the space-charge
 * configuration. */

#ifndef OPALX_SPACE_CHARGE_CONFIG_BUILDER_H
#define OPALX_SPACE_CHARGE_CONFIG_BUILDER_H

#include "SpaceCharge/SpaceChargeConfig.h"

#include <optional>
#include <vector>

#include "Track/ObjectValues.h"

namespace opalx::spacecharge {

    /**
     * @brief Make the validated space-charge configuration of a RUN from its descriptions.
     * @param fieldSolver The FIELDSOLVER of the RUN.
     * @param binning The BINNING that FIELDSOLVER BINS names, if any.
     * @param emissionSources The emission sources of all beams (Dirichlet-plane settings).
     * @param distributions The distributions these sources name; a FLATTOP or OPALFLATTOP one
     *        changes the decomposition used when the particle layout is rebuilt.
     * @param repartFreq Steps between load-balancing checks (OPTION REPARTFREQ); zero or less
     *        disables them.
     * @param loadBalancingThreshold OPTION LOADBALANCINGTHRESHOLD.
     * @param loadBalancingPairWeight OPTION LOADBALANCINGPAIRWEIGHT.
     */
    [[nodiscard]] SpaceChargeConfig buildSpaceChargeConfig(
            const ObjectValues& fieldSolver, const std::optional<ObjectValues>& binning,
            const std::vector<ObjectValues>& emissionSources,
            const std::vector<ObjectValues>& distributions, int repartFreq,
            double loadBalancingThreshold, double loadBalancingPairWeight);

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_CONFIG_BUILDER_H
