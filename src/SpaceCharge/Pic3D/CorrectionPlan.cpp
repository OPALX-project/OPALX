/**
 * @file CorrectionPlan.cpp
 * @brief Explicit instantiation of production Cartesian PIC correction planning.
 */

#include "SpaceCharge/Pic3D/CorrectionPlan.h"

#include <utility>

namespace opalx::spacecharge {

    CorrectionPlan::CorrectionPlan(CorrectionConfig config, bool binned)
        : config_m(std::move(config)), binned_m(binned) {}

    PreparedCorrection CorrectionPlan::prepare(
            const RequestedPhysics& requested, std::size_t step) const {
        PreparedCorrection prepared;
        prepared.configuredCorrection = {config_m.kind(), config_m.planeZ()};
        prepared.maximumSteps         = config_m.maximumSteps();
        prepared.correctionExpired    = config_m.enabled() && config_m.maximumSteps() != 0
                                     && step >= config_m.maximumSteps();

        const CorrectionKind active = requested.correction.kind;
        if (!binned_m) {
            prepared.passes[0]      = active == CorrectionKind::ImageCharge
                                              ? SolvePass::DirectPrimaryWithImage
                                              : SolvePass::DirectPrimary;
            prepared.passCount      = 1;
            prepared.shiftedIgnored = active == CorrectionKind::ShiftedGreen;
            return prepared;
        }

        prepared.passes[prepared.passCount++] = SolvePass::BinnedPrimary;
        if (active == CorrectionKind::ImageCharge) {
            prepared.passes[prepared.passCount++] = SolvePass::BinnedImage;
        } else if (active == CorrectionKind::ShiftedGreen) {
            prepared.passes[prepared.passCount++] = SolvePass::BinnedShiftedGreen;
        }
        return prepared;
    }

}  // namespace opalx::spacecharge
