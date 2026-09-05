/**
 * @file SpaceChargeConfig.cpp
 * @brief Validates space-charge configuration and derives particle-domain setup.
 */

#include "SpaceCharge/SpaceChargeConfig.h"

#include "Utilities/OpalException.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>

namespace opalx::spacecharge {
    namespace {

        void validateGrid(const CartesianGridConfig& grid) {
            if (std::any_of(grid.meshSize.begin(), grid.meshSize.end(), [](std::size_t size) {
                    return size == 0
                           || size > static_cast<std::size_t>(std::numeric_limits<int>::max());
                })) {
                throw OpalException(
                        "validateSpaceChargeConfig",
                        "Mesh dimensions must fit in a positive IPPL Index.");
            }
            if (!std::isfinite(grid.boundingBoxIncreasePercent)
                || grid.boundingBoxIncreasePercent < 0.0) {
                throw OpalException(
                        "validateSpaceChargeConfig",
                        "The bounding-box increase must not be negative.");
            }
        }

        void validateCorrection(const CorrectionConfig& correction) {
            if (!std::isfinite(correction.planeZ)) {
                throw OpalException(
                        "validateSpaceChargeConfig", "Source-plane position must be finite.");
            }
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
            if (!std::isfinite(binning.desiredWidth) || !(binning.desiredWidth > 0.0)
                || !std::isfinite(binning.alpha) || !std::isfinite(binning.beta)) {
                throw OpalException(
                        "validateSpaceChargeConfig",
                        "Binning needs a finite positive width and finite cost parameters.");
            }
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
            if (!std::isfinite(config.loadBalancingThreshold)
                || config.loadBalancingThreshold < 0.0) {
                throw OpalException(
                        "validateSpaceChargeConfig",
                        "The load-balancing threshold must not be negative.");
            }
            validatePoissonSolverConfig(makePoissonSolverConfig(config));
            if (config.backend == PoissonSolverType::P3M
                && (config.binning.has_value() || config.correction.enabled())) {
                throw OpalException(
                        "validateSpaceChargeConfig",
                        "P3M does not support binning or source-plane corrections.");
            }
            if (config.correction.enabled() && config.backend != PoissonSolverType::Open) {
                throw OpalException(
                        "validateSpaceChargeConfig",
                        "Source-plane corrections require the OPEN Poisson backend.");
            }
            if (config.correction.planeDumpFrequency != 0 && config.binning.has_value()) {
                throw OpalException(
                        "validateSpaceChargeConfig",
                        "Source-plane potential diagnostics are not supported with binning.");
            }
        }

        void validateFFT2D5(const FFT2D5Config& config) {
            validateGrid(config.grid);
            if (std::any_of(
                        config.grid.decomposition.begin(), config.grid.decomposition.end(),
                        [](bool parallel) {
                            return parallel;
                        })) {
                throw OpalException(
                        "validateSpaceChargeConfig",
                        "FFT2D5 uses serial slices; set PARFFTX/PARFFTY/PARFFTZ to FALSE.");
            }
            if (!std::isfinite(config.pipeSizeX) || !std::isfinite(config.pipeSizeY)
                || !std::isfinite(config.beamRadius) || !(config.pipeSizeX > 0.0)
                || !(config.pipeSizeY > 0.0) || !(config.beamRadius > 0.0)) {
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

    PoissonSolverConfig makePoissonSolverConfig(const CartesianPICConfig& config) {
        return {config.backend, config.greenFunction, config.p3mCutoff, config.boundaryConditions};
    }

    void validatePoissonSolverConfig(const PoissonSolverConfig& config) {
        const auto& boundaries = config.boundaryConditions;
        const bool allOpen = std::all_of(boundaries.begin(), boundaries.end(), [](auto boundary) {
            return boundary == FieldBoundaryCondition::Open;
        });
        const bool allPeriodic =
                std::all_of(boundaries.begin(), boundaries.end(), [](auto boundary) {
                    return boundary == FieldBoundaryCondition::Periodic;
                });
        if (config.type == PoissonSolverType::ConjugateGradient) {
            throw OpalException(
                    "validateSpaceChargeConfig",
                    "The CG backend is recognized but not implemented.");
        }
        if (!allOpen && !allPeriodic) {
            throw OpalException(
                    "validateSpaceChargeConfig",
                    "Implemented Poisson backends require uniform OPEN or PERIODIC "
                    "domain boundaries. DIRICHLET mesh faces and mixed boundaries are "
                    "not implemented. Source-plane corrections are configured separately.");
        }
        switch (config.type) {
            case PoissonSolverType::None:
            case PoissonSolverType::P3M:
                break;
            case PoissonSolverType::Open:
                if (!allOpen) {
                    throw OpalException(
                            "validateSpaceChargeConfig",
                            "TYPE=OPEN requires OPEN domain boundaries.");
                }
                break;
            case PoissonSolverType::PeriodicFFT:
                if (!allPeriodic) {
                    throw OpalException(
                            "validateSpaceChargeConfig",
                            "TYPE=FFT requires PERIODIC domain boundaries.");
                }
                break;
            default:
                throw OpalException("validateSpaceChargeConfig", "Unknown Poisson backend.");
        }
        if (config.greenFunction != GreenFunctionType::Standard
            && config.greenFunction != GreenFunctionType::Integrated) {
            throw OpalException("validateSpaceChargeConfig", "Unknown open Green function.");
        }
        if (!std::isfinite(config.p3mCutoff)
            || (config.type == PoissonSolverType::P3M ? config.p3mCutoff <= 0.0
                                                      : config.p3mCutoff != 0.0)) {
            throw OpalException(
                    "validateSpaceChargeConfig",
                    "P3M requires a finite positive cutoff; other backends require no cutoff.");
        }
    }

    void validateSpaceChargeConfig(const SpaceChargeConfig& config) {
        std::visit(
                [](const auto& selected) {
                    using Config = std::decay_t<decltype(selected)>;
                    if constexpr (std::is_same_v<Config, CartesianPICConfig>) {
                        validateCartesian(selected);
                    } else if constexpr (std::is_same_v<Config, FFT2D5Config>) {
                        validateFFT2D5(selected);
                    } else {
                        static_assert(
                                std::is_same_v<Config, void>, "Add validation for this algorithm");
                    }
                },
                config);
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
