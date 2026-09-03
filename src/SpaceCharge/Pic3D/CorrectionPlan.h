/**
 * @file CorrectionPlan.h
 * @brief Expands validated correction requests into ordered Cartesian PIC pass tags.
 */

#ifndef OPALX_SPACE_CHARGE_PIC3D_CORRECTION_PLAN_H
#define OPALX_SPACE_CHARGE_PIC3D_CORRECTION_PLAN_H

#include "SpaceCharge/Pic3D/SolvePass.h"
#include "SpaceCharge/SelfFieldConfig.h"
#include "SpaceCharge/SolveContext.h"

#include <array>
#include <cstddef>

namespace opalx::spacecharge {

    struct PreparedCorrection final {
        static constexpr std::size_t maximumPassCount = 2;

        [[nodiscard]] const SolvePass* begin() const { return passes.data(); }
        [[nodiscard]] const SolvePass* end() const { return passes.data() + passCount; }

        std::array<SolvePass, maximumPassCount> passes{};
        std::size_t passCount = 0;
        CorrectionRequest configuredCorrection;
        std::size_t maximumSteps = 0;
        bool correctionExpired   = false;
        bool shiftedIgnored      = false;
    };

    /**
     * @brief Expands a request already validated by SelfFieldSystem into ordered pass tags.
     *
     * This class retains only immutable correction configuration and whether traversal is binned.
     * It performs no second request-validation layer.
     */
    class CorrectionPlan final {
    public:
        CorrectionPlan(CorrectionConfig config, bool binned);

        [[nodiscard]] PreparedCorrection prepare(
                const RequestedPhysics& requested, std::size_t step) const;

    private:
        CorrectionConfig config_m;
        bool binned_m = false;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_PIC3D_CORRECTION_PLAN_H
