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
        if (!(values_m.loadBalancingThreshold >= 0.0)) {
            throw OpalException(
                    "Pic3DConfig::Pic3DConfig",
                    "The load-balancing threshold must not be negative.");
        }
        if (values_m.backend == PoissonBackendKind::P3M && !(values_m.p3mCutoff > 0.0)) {
            throw OpalException(
                    "Pic3DConfig::Pic3DConfig", "The P3M cutoff radius must be positive.");
        }
        if (values_m.backend != PoissonBackendKind::P3M && values_m.p3mCutoff != 0.0) {
            throw OpalException(
                    "Pic3DConfig::Pic3DConfig",
                    "A P3M cutoff may only be configured for the P3M backend.");
        }
        if (values_m.backend == PoissonBackendKind::P3M) {
            const bool allOpen = std::all_of(
                    values_m.boundaryConditions.begin(), values_m.boundaryConditions.end(),
                    [](BoundaryConditionKind kind) {
                        return kind == BoundaryConditionKind::Open;
                    });
            const bool allPeriodic = std::all_of(
                    values_m.boundaryConditions.begin(), values_m.boundaryConditions.end(),
                    [](BoundaryConditionKind kind) {
                        return kind == BoundaryConditionKind::Periodic;
                    });
            if (!allOpen && !allPeriodic) {
                throw OpalException(
                        "Pic3DConfig::Pic3DConfig",
                        "P3M requires uniform OPEN or PERIODIC boundary conditions.");
            }
            if (values_m.binning.has_value()) {
                throw OpalException("Pic3DConfig::Pic3DConfig", "P3M does not support binning.");
            }
            if (values_m.correction.enabled()) {
                throw OpalException(
                        "Pic3DConfig::Pic3DConfig",
                        "P3M does not support source-plane corrections.");
            }
        }
    }

    Pic2d5Config::Pic2d5Config(Pic2d5ConfigValues values) : values_m(std::move(values)) {
        const bool hasEmptyDimension = std::any_of(
                values_m.meshSize.begin(), values_m.meshSize.end(), [](std::size_t size) {
                    return size == 0;
                });
        if (hasEmptyDimension) {
            throw OpalException("Pic2d5Config::Pic2d5Config", "Mesh dimensions must be positive.");
        }
        if (!(values_m.pipeSizeX > 0.0) || !(values_m.pipeSizeY > 0.0)) {
            throw OpalException(
                    "Pic2d5Config::Pic2d5Config",
                    "The transverse pipe dimensions must be positive.");
        }
        if (!(values_m.beamRadius > 0.0)) {
            throw OpalException("Pic2d5Config::Pic2d5Config", "The beam radius must be positive.");
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
                    } else {
                        return SelfFieldAlgorithmKind::Pic2d5;
                    }
                },
                algorithmConfig_m);
    }

    ParticleLayoutConfig SelfFieldConfig::particleLayoutConfig() const {
        return std::visit(
                [](const auto& config) {
                    using Config = std::decay_t<decltype(config)>;
                    ParticleLayoutConfig result;
                    if constexpr (std::is_same_v<Config, Pic3DConfig>) {
                        result.periodic = std::all_of(
                                config.boundaryConditions().begin(),
                                config.boundaryConditions().end(), [](BoundaryConditionKind kind) {
                                    return kind == BoundaryConditionKind::Periodic;
                                });
                        if (config.backend() == PoissonBackendKind::P3M) {
                            result.kind          = ParticleLayoutKind::SpatialOverlap;
                            result.overlapCutoff = config.p3mCutoff();
                        }
                    }
                    return result;
                },
                algorithmConfig_m);
    }

}  // namespace opalx::spacecharge
