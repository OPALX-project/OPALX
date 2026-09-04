/**
 * @file SpaceChargeConfig.cpp
 * @brief Implements immutable space-charge configuration snapshots.
 */

#include "SpaceCharge/SpaceChargeConfig.h"

#include "Utilities/OpalException.h"

#include <algorithm>
#include <type_traits>

namespace opalx::spacecharge {

    CorrectionConfig::CorrectionConfig(CorrectionConfig::Parameters parameters)
        : parameters_m(std::move(parameters)) {
        if (parameters_m.kind != SpaceChargeCorrectionType::ImageCharge
            && parameters_m.planeDumpFrequency != 0) {
            throw OpalException(
                    "CorrectionConfig::CorrectionConfig",
                    "ZEROFACE plane dumping is only available for the image-charge correction.");
        }
    }

    BinningConfig::BinningConfig(BinningConfig::Parameters parameters)
        : parameters_m(std::move(parameters)) {
        if (parameters_m.maximumBins == 0) {
            throw OpalException(
                    "BinningConfig::BinningConfig", "Binning requires at least one bin.");
        }
        if (!parameters_m.dumpFile.empty() && parameters_m.dumpFrequency == 0) {
            throw OpalException(
                    "BinningConfig::BinningConfig",
                    "DUMPBINSFREQ must be positive when bin dumping is enabled.");
        }
    }

    CartesianPICConfig::CartesianPICConfig(CartesianPICConfig::Parameters parameters)
        : parameters_m(std::move(parameters)) {
        const bool hasEmptyDimension = std::any_of(
                parameters_m.meshSize.begin(), parameters_m.meshSize.end(), [](std::size_t size) {
                    return size == 0;
                });
        if (hasEmptyDimension) {
            throw OpalException(
                    "CartesianPICConfig::CartesianPICConfig", "Mesh dimensions must be positive.");
        }
        if (parameters_m.boundingBoxIncreasePercent < 0.0) {
            throw OpalException(
                    "CartesianPICConfig::CartesianPICConfig",
                    "The bounding-box increase must not be negative.");
        }
        if (!(parameters_m.loadBalancingThreshold >= 0.0)) {
            throw OpalException(
                    "CartesianPICConfig::CartesianPICConfig",
                    "The load-balancing threshold must not be negative.");
        }
        if (parameters_m.backend == PoissonSolverType::P3M && !(parameters_m.p3mCutoff > 0.0)) {
            throw OpalException(
                    "CartesianPICConfig::CartesianPICConfig",
                    "The P3M cutoff radius must be positive.");
        }
        if (parameters_m.backend != PoissonSolverType::P3M && parameters_m.p3mCutoff != 0.0) {
            throw OpalException(
                    "CartesianPICConfig::CartesianPICConfig",
                    "A P3M cutoff may only be configured for the P3M backend.");
        }
        if (parameters_m.backend == PoissonSolverType::P3M) {
            const bool allOpen = std::all_of(
                    parameters_m.boundaryConditions.begin(), parameters_m.boundaryConditions.end(),
                    [](FieldBoundaryCondition kind) {
                        return kind == FieldBoundaryCondition::Open;
                    });
            const bool allPeriodic = std::all_of(
                    parameters_m.boundaryConditions.begin(), parameters_m.boundaryConditions.end(),
                    [](FieldBoundaryCondition kind) {
                        return kind == FieldBoundaryCondition::Periodic;
                    });
            if (!allOpen && !allPeriodic) {
                throw OpalException(
                        "CartesianPICConfig::CartesianPICConfig",
                        "P3M requires uniform OPEN or PERIODIC boundary conditions.");
            }
            if (parameters_m.binning.has_value()) {
                throw OpalException(
                        "CartesianPICConfig::CartesianPICConfig", "P3M does not support binning.");
            }
            if (parameters_m.correction.enabled()) {
                throw OpalException(
                        "CartesianPICConfig::CartesianPICConfig",
                        "P3M does not support source-plane corrections.");
            }
        }
        if (parameters_m.correction.kind() == SpaceChargeCorrectionType::ShiftedGreen
            && parameters_m.backend != PoissonSolverType::Open) {
            throw OpalException(
                    "CartesianPICConfig::CartesianPICConfig",
                    "The shifted-Green correction requires the OPEN backend.");
        }
    }

    FFT2D5Config::FFT2D5Config(FFT2D5Config::Parameters parameters)
        : parameters_m(std::move(parameters)) {
        const bool hasEmptyDimension = std::any_of(
                parameters_m.meshSize.begin(), parameters_m.meshSize.end(), [](std::size_t size) {
                    return size == 0;
                });
        if (hasEmptyDimension) {
            throw OpalException("FFT2D5Config::FFT2D5Config", "Mesh dimensions must be positive.");
        }
        if (!(parameters_m.pipeSizeX > 0.0) || !(parameters_m.pipeSizeY > 0.0)) {
            throw OpalException(
                    "FFT2D5Config::FFT2D5Config",
                    "The transverse pipe dimensions must be positive.");
        }
        if (!(parameters_m.beamRadius > 0.0)) {
            throw OpalException("FFT2D5Config::FFT2D5Config", "The beam radius must be positive.");
        }
    }

    SpaceChargeConfig::SpaceChargeConfig(
            SpaceChargeAlgorithmConfig algorithmConfig,
            CartesianDomainConfig3D cartesianDomainConfig)
        : algorithmConfig_m(std::move(algorithmConfig)),
          cartesianDomainConfig_m(std::move(cartesianDomainConfig)) {
        const bool invalidMesh = std::any_of(
                cartesianDomainConfig_m.meshSize.begin(), cartesianDomainConfig_m.meshSize.end(),
                [](std::size_t extent) {
                    return extent == 0;
                });
        if (invalidMesh || cartesianDomainConfig_m.boundingBoxIncreasePercent < 0.0) {
            throw OpalException(
                    "SpaceChargeConfig::SpaceChargeConfig",
                    "The particle-storage mesh and bounding box must be valid.");
        }

        std::visit(
                [this](const auto& config) {
                    using Config = std::decay_t<decltype(config)>;
                    if (config.meshSize() != cartesianDomainConfig_m.meshSize) {
                        throw OpalException(
                                "SpaceChargeConfig::SpaceChargeConfig",
                                "The solver and particle-storage mesh sizes differ.");
                    }
                    if constexpr (std::is_same_v<Config, CartesianPICConfig>) {
                        const bool p3m         = config.backend() == PoissonSolverType::P3M;
                        const bool allPeriodic = std::all_of(
                                config.boundaryConditions().begin(),
                                config.boundaryConditions().end(), [](FieldBoundaryCondition kind) {
                                    return kind == FieldBoundaryCondition::Periodic;
                                });
                        const bool wrongLayout = p3m
                                                 != (cartesianDomainConfig_m.layoutType
                                                     == ParticleLayoutType::SpatialOverlap);
                        const bool wrongCutoff =
                                p3m ? cartesianDomainConfig_m.overlapCutoff != config.p3mCutoff()
                                    : cartesianDomainConfig_m.overlapCutoff != 0.0;
                        if (config.parallelDimensions() != cartesianDomainConfig_m.decomposition
                            || config.boundingBoxIncreasePercent()
                                       != cartesianDomainConfig_m.boundingBoxIncreasePercent
                            || allPeriodic != cartesianDomainConfig_m.periodicParticleBoundary
                            || wrongLayout || wrongCutoff) {
                            throw OpalException(
                                    "SpaceChargeConfig::SpaceChargeConfig",
                                    "The Cartesian solver and particle-storage setup differ.");
                        }
                    } else {
                        if (cartesianDomainConfig_m.layoutType != ParticleLayoutType::Spatial
                            || cartesianDomainConfig_m.overlapCutoff != 0.0
                            || cartesianDomainConfig_m.periodicParticleBoundary) {
                            throw OpalException(
                                    "SpaceChargeConfig::SpaceChargeConfig",
                                    "FFT2D5 requires the standard non-wrapping particle layout.");
                        }
                    }
                },
                algorithmConfig_m);
    }

    SpaceChargeAlgorithmType SpaceChargeConfig::algorithmType() const {
        return std::visit(
                [](const auto& config) {
                    using Config = std::decay_t<decltype(config)>;
                    if constexpr (std::is_same_v<Config, CartesianPICConfig>) {
                        return SpaceChargeAlgorithmType::CartesianPIC;
                    } else {
                        return SpaceChargeAlgorithmType::FFT2D5;
                    }
                },
                algorithmConfig_m);
    }

}  // namespace opalx::spacecharge
