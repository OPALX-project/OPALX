/**
 * @file SelfFieldConfig.cpp
 * @brief Implements immutable self-field configuration snapshots.
 */

#include "SpaceCharge/SelfFieldConfig.h"

#include "Utilities/OpalException.h"

#include <algorithm>
#include <type_traits>

namespace opalx::spacecharge {

    CorrectionConfig::CorrectionConfig(CorrectionConfigValues values)
        : values_m(std::move(values)) {
        if (values_m.kind != CorrectionKind::ImageCharge && values_m.planeDumpFrequency != 0) {
            throw OpalException(
                    "CorrectionConfig::CorrectionConfig",
                    "ZEROFACE plane dumping is only available for the image-charge correction.");
        }
    }

    BinningConfig::BinningConfig(BinningConfigValues values) : values_m(std::move(values)) {
        if (values_m.maximumBins == 0) {
            throw OpalException(
                    "BinningConfig::BinningConfig", "Binning requires at least one bin.");
        }
        if (!values_m.dumpFile.empty() && values_m.dumpFrequency == 0) {
            throw OpalException(
                    "BinningConfig::BinningConfig",
                    "DUMPBINSFREQ must be positive when bin dumping is enabled.");
        }
    }

    Pic3DConfig::Pic3DConfig(Pic3DConfigValues values) : values_m(std::move(values)) {
        const bool hasEmptyDimension = std::any_of(
                values_m.meshSize.begin(), values_m.meshSize.end(), [](std::size_t size) {
                    return size == 0;
                });
        if (hasEmptyDimension) {
            throw OpalException("Pic3DConfig::Pic3DConfig", "Mesh dimensions must be positive.");
        }
        if (values_m.boundingBoxIncreasePercent < 0.0) {
            throw OpalException(
                    "Pic3DConfig::Pic3DConfig", "The bounding-box increase must not be negative.");
        }
    }

    SelfFieldConfig::SelfFieldConfig(SelfFieldAlgorithmConfig algorithmConfig)
        : algorithmConfig_m(std::move(algorithmConfig)) {}

    SelfFieldAlgorithmKind SelfFieldConfig::algorithmKind() const {
        return std::visit(
                [](const auto& config) {
                    using Config = std::decay_t<decltype(config)>;
                    if constexpr (std::is_same_v<Config, Pic3DConfig>) {
                        return SelfFieldAlgorithmKind::Pic3D;
                    }
                },
                algorithmConfig_m);
    }

}  // namespace opalx::spacecharge
