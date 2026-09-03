/**
 * @file SelfFieldRequestPolicy.cpp
 * @brief Implements per-step request resolution from immutable configuration.
 */

#include "SpaceCharge/SelfFieldRequestPolicy.h"

#include <type_traits>

namespace opalx::spacecharge {

    SelfFieldRequestPolicy::SelfFieldRequestPolicy(const SelfFieldConfig& config) {
        std::visit(
                [this](const auto& algorithmConfig) {
                    using Config = std::decay_t<decltype(algorithmConfig)>;
                    if constexpr (std::is_same_v<Config, Pic3DConfig>) {
                        const CorrectionConfig& correction = algorithmConfig.correction();
                        correction_m             = {correction.kind(), correction.planeZ()};
                        correctionMaximumSteps_m = correction.maximumSteps();
                        useBinning_m             = algorithmConfig.binning().has_value();
                        writePotential_m         = correction.kind() == CorrectionKind::ImageCharge
                                           && correction.planeDumpFrequency() != 0;
                    }
                },
                config.algorithmConfig());
    }

    RequestedPhysics SelfFieldRequestPolicy::forStep(std::size_t step) const {
        RequestedPhysics requested;
        requested.useBinning = useBinning_m;
        const bool correctionActive =
                correction_m.kind != CorrectionKind::None
                && (correctionMaximumSteps_m == 0 || step < correctionMaximumSteps_m);
        if (correctionActive) {
            requested.correction     = correction_m;
            requested.writePotential = writePotential_m;
        }
        return requested;
    }

}  // namespace opalx::spacecharge
