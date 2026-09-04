/**
 * @file CorrectionPassSchedule.cpp
 * @brief Implements the Cartesian PIC correction-pass schedule.
 */

#include "SpaceCharge/CartesianPIC/CorrectionPassSchedule.h"

#include <utility>

namespace opalx::spacecharge {

    CorrectionPassSchedule::CorrectionPassSchedule(CorrectionConfig config, bool binned)
        : config_m(std::move(config)), binned_m(binned) {}

    CorrectionPassSequence CorrectionPassSchedule::passesForStep(
            const SpaceChargeRequest& requested, std::size_t step) const {
        CorrectionPassSequence sequence;
        sequence.configuredCorrection = {config_m.kind(), config_m.planeZ()};
        sequence.maximumSteps         = config_m.maximumSteps();
        sequence.correctionExpired    = config_m.enabled() && config_m.maximumSteps() != 0
                                     && step >= config_m.maximumSteps();

        const SpaceChargeCorrectionType active = requested.correction.kind;
        if (!binned_m) {
            sequence.passes[0]      = active == SpaceChargeCorrectionType::ImageCharge
                                              ? CartesianPICPass::DirectPrimaryWithImage
                                              : CartesianPICPass::DirectPrimary;
            sequence.passCount      = 1;
            sequence.shiftedIgnored = active == SpaceChargeCorrectionType::ShiftedGreen;
            return sequence;
        }

        sequence.passes[sequence.passCount++] = CartesianPICPass::BinnedPrimary;
        if (active == SpaceChargeCorrectionType::ImageCharge) {
            sequence.passes[sequence.passCount++] = CartesianPICPass::BinnedImage;
        } else if (active == SpaceChargeCorrectionType::ShiftedGreen) {
            sequence.passes[sequence.passCount++] = CartesianPICPass::BinnedShiftedGreen;
        }
        return sequence;
    }

}  // namespace opalx::spacecharge
