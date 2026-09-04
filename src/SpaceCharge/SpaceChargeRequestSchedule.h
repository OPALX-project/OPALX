/**
 * @file SpaceChargeRequestSchedule.h
 * @brief Resolves immutable run configuration into per-step physics requests.
 */

#ifndef OPALX_SPACE_CHARGE_REQUEST_SCHEDULE_H
#define OPALX_SPACE_CHARGE_REQUEST_SCHEDULE_H

#include "SpaceCharge/SpaceChargeConfig.h"
#include "SpaceCharge/SpaceChargeSolveContext.h"

#include <cstddef>

namespace opalx::spacecharge {

    /** @brief Tracker-owned snapshot of step-dependent space-charge request rules. */
    class SpaceChargeRequestSchedule {
    public:
        SpaceChargeRequestSchedule() = default;
        explicit SpaceChargeRequestSchedule(const SpaceChargeConfig& config);

        [[nodiscard]] SpaceChargeRequest requestForStep(std::size_t step) const;
        [[nodiscard]] SpaceChargeCorrectionRequest configuredCorrection() const {
            return correction_m;
        }

    private:
        SpaceChargeCorrectionRequest correction_m;
        std::size_t correctionMaximumSteps_m = 0;
        bool useBinning_m                    = false;
        bool writePotential_m                = false;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_REQUEST_SCHEDULE_H
