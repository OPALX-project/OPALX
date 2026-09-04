/**
 * @file SpaceChargeConfig.cpp
 * @brief Validates space-charge configuration and derives particle-domain setup.
 */

#include "SpaceCharge/SpaceChargeConfig.h"

#include "Utilities/OpalException.h"

#include <algorithm>
#include <type_traits>

namespace opalx::spacecharge {
    namespace {

        void validateGrid(const CartesianGridConfig& grid) {
            if (std::any_of(grid.meshSize.begin(), grid.meshSize.end(), [](std::size_t size) {
                    return size == 0;
                })) {
                throw OpalException(
                        "validateSpaceChargeConfig", "Mesh dimensions must be positive.");
            }
            if (grid.boundingBoxIncreasePercent < 0.0) {
                throw OpalException(
                        "validateSpaceChargeConfig",
                        "The bounding-box increase must not be negative.");
            }
        }

        void validateCorrection(const CorrectionConfig& correction) {
            switch (correction.kind) {
                case SpaceChargeCorrectionType::None:
                case SpaceChargeCorrectionType::ImageCharge:
                case SpaceChargeCorrectionType::ShiftedGreen:
                    break;
                default:
                    throw OpalException(
                            "validateSpaceChargeConfig", "Unknown source-plane correction.");
            }
            if (correction.kind != SpaceChargeCorrectionType::ImageCharge
                && correction.planeDumpFrequency != 0) {
                throw OpalException(
                        "validateSpaceChargeConfig",
                        "Source-plane potential output requires image-charge correction.");
            }
        }

        void validateBinning(const BinningConfig& binning) {
            if (binning.maximumBins == 0) {
                throw OpalException(
                        "validateSpaceChargeConfig", "Binning requires at least one bin.");
            }
            if (!binning.dumpFile.empty() && binning.dumpFrequency == 0) {
                throw OpalException(
                        "validateSpaceChargeConfig",
                        "DUMPBINSFREQ must be positive when bin dumping is enabled.");
            }
            if (binning.parameter != BinningVariable::VelocityZ
                && binning.parameter != BinningVariable::GammaZ) {
                throw OpalException(
                        "validateSpaceChargeConfig", "Binning supports only VELOCITYZ and GAMMAZ.");
            }
        }

        void validateCartesian(const CartesianPICConfig& config) {
            validateGrid(config.grid);
            validateCorrection(config.correction);
            if (config.binning.has_value()) {
                validateBinning(*config.binning);
            }
            if (!(config.loadBalancingThreshold >= 0.0)) {
                throw OpalException(
                        "validateSpaceChargeConfig",
                        "The load-balancing threshold must not be negative.");
            }
            switch (config.backend) {
                case PoissonSolverType::None:
                case PoissonSolverType::PeriodicFFT:
                case PoissonSolverType::Open:
                case PoissonSolverType::ConjugateGradient:
                case PoissonSolverType::P3M:
                    break;
                default:
                    throw OpalException("validateSpaceChargeConfig", "Unknown Poisson backend.");
            }
            switch (config.greenFunction) {
                case GreenFunctionType::Standard:
                case GreenFunctionType::Integrated:
                    break;
                default:
                    throw OpalException(
                            "validateSpaceChargeConfig", "Unknown open Green function.");
            }
            for (const FieldBoundaryCondition boundary : config.boundaryConditions) {
                if (boundary != FieldBoundaryCondition::Open
                    && boundary != FieldBoundaryCondition::Dirichlet
                    && boundary != FieldBoundaryCondition::Periodic) {
                    throw OpalException(
                            "validateSpaceChargeConfig", "Unknown field boundary condition.");
                }
            }
            if (config.backend == PoissonSolverType::P3M && !(config.p3mCutoff > 0.0)) {
                throw OpalException(
                        "validateSpaceChargeConfig", "The P3M cutoff radius must be positive.");
            }
            if (config.backend != PoissonSolverType::P3M && config.p3mCutoff != 0.0) {
                throw OpalException(
                        "validateSpaceChargeConfig",
                        "A P3M cutoff may only be configured for the P3M backend.");
            }
            const bool allOpen = std::all_of(
                    config.boundaryConditions.begin(), config.boundaryConditions.end(),
                    [](FieldBoundaryCondition boundary) {
                        return boundary == FieldBoundaryCondition::Open;
                    });
            const bool allPeriodic = std::all_of(
                    config.boundaryConditions.begin(), config.boundaryConditions.end(),
                    [](FieldBoundaryCondition boundary) {
                        return boundary == FieldBoundaryCondition::Periodic;
                    });
            if (config.backend == PoissonSolverType::P3M) {
                if (!allOpen && !allPeriodic) {
                    throw OpalException(
                            "validateSpaceChargeConfig",
                            "P3M requires uniform OPEN or PERIODIC boundary conditions.");
                }
                if (config.binning.has_value() || config.correction.enabled()) {
                    throw OpalException(
                            "validateSpaceChargeConfig",
                            "P3M does not support binning or source-plane corrections.");
                }
            }
            if (config.correction.kind == SpaceChargeCorrectionType::ShiftedGreen
                && (config.backend != PoissonSolverType::Open || !config.binning.has_value())) {
                throw OpalException(
                        "validateSpaceChargeConfig",
                        "Shifted-Green correction requires OPEN with particle binning.");
            }
            if (config.correction.planeDumpFrequency != 0 && config.binning.has_value()) {
                throw OpalException(
                        "validateSpaceChargeConfig",
                        "Source-plane potential diagnostics are not supported with binning.");
            }
        }

        void validateFFT2D5(const FFT2D5Config& config) {
            validateGrid(config.grid);
            if (!(config.pipeSizeX > 0.0) || !(config.pipeSizeY > 0.0)
                || !(config.beamRadius > 0.0)) {
                throw OpalException(
                        "validateSpaceChargeConfig",
                        "FFT2D5 pipe dimensions and beam radius must be positive.");
            }
            switch (config.longitudinalFieldMode) {
                case FFT2D5LongitudinalFieldMode::Open:
                case FFT2D5LongitudinalFieldMode::Cylindrical:
                case FFT2D5LongitudinalFieldMode::Plates:
                case FFT2D5LongitudinalFieldMode::None:
                    break;
                default:
                    throw OpalException(
                            "validateSpaceChargeConfig", "Unknown FFT2D5 longitudinal mode.");
            }
        }

    }  // namespace

    void validateSpaceChargeConfig(const SpaceChargeConfig& config) {
        std::visit(
                [](const auto& selected) {
                    using Config = std::decay_t<decltype(selected)>;
                    if constexpr (std::is_same_v<Config, CartesianPICConfig>) {
                        validateCartesian(selected);
                    } else {
                        validateFFT2D5(selected);
                    }
                },
                config);
    }

    SpaceChargeAlgorithmType algorithmType(const SpaceChargeConfig& config) {
        return std::holds_alternative<CartesianPICConfig>(config)
                       ? SpaceChargeAlgorithmType::CartesianPIC
                       : SpaceChargeAlgorithmType::FFT2D5;
    }

    CartesianDomainConfig3D makeCartesianDomainConfig(const SpaceChargeConfig& config) {
        CartesianDomainConfig3D domain;
        std::visit(
                [&domain](const auto& selected) {
                    domain.meshSize                   = selected.grid.meshSize;
                    domain.decomposition              = selected.grid.decomposition;
                    domain.boundingBoxIncreasePercent = selected.grid.boundingBoxIncreasePercent;
                    using Config                      = std::decay_t<decltype(selected)>;
                    if constexpr (std::is_same_v<Config, CartesianPICConfig>) {
                        domain.periodicParticleBoundary = std::all_of(
                                selected.boundaryConditions.begin(),
                                selected.boundaryConditions.end(),
                                [](FieldBoundaryCondition boundary) {
                                    return boundary == FieldBoundaryCondition::Periodic;
                                });
                        if (selected.backend == PoissonSolverType::P3M) {
                            domain.layoutType    = ParticleLayoutType::SpatialOverlap;
                            domain.overlapCutoff = selected.p3mCutoff;
                        }
                    }
                },
                config);
        return domain;
    }

}  // namespace opalx::spacecharge
