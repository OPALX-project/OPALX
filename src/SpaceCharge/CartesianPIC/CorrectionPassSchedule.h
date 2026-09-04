/**
 * @file CorrectionPassSchedule.h
 * @brief Expands validated correction requests into ordered Cartesian PIC pass tags.
 */

#ifndef OPALX_SPACE_CHARGE_CARTESIAN_PIC_CORRECTION_PASS_SCHEDULE_H
#define OPALX_SPACE_CHARGE_CARTESIAN_PIC_CORRECTION_PASS_SCHEDULE_H

#include "SpaceCharge/CartesianPIC/CartesianPICPass.h"
#include "SpaceCharge/SpaceChargeConfig.h"
#include "SpaceCharge/SpaceChargeSolveContext.h"

#include <array>
#include <cstddef>

namespace opalx::spacecharge {

    struct CorrectionPassSequence final {
        static constexpr std::size_t maximumPassCount = 2;

        [[nodiscard]] const CartesianPICPass* begin() const { return passes.data(); }
        [[nodiscard]] const CartesianPICPass* end() const { return passes.data() + passCount; }

        std::array<CartesianPICPass, maximumPassCount> passes{};
        std::size_t passCount = 0;
        SpaceChargeCorrectionRequest configuredCorrection;
        std::size_t maximumSteps = 0;
        bool correctionExpired   = false;
        bool shiftedIgnored      = false;
    };

    /**
     * @brief Expands a request already validated by SpaceChargeSolver into ordered pass tags.
     *
     * This class retains only immutable correction configuration and whether traversal is binned.
     * It performs no second request-validation layer.
     */
    class CorrectionPassSchedule final {
    public:
        CorrectionPassSchedule(CorrectionConfig config, bool binned);

        [[nodiscard]] CorrectionPassSequence passesForStep(
                const SpaceChargeRequest& requested, std::size_t step) const;

    private:
        CorrectionConfig config_m;
        bool binned_m = false;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_CARTESIAN_PIC_CORRECTION_PASS_SCHEDULE_H
