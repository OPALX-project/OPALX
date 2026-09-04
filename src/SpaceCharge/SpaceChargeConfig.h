/**
 * @file SpaceChargeConfig.h
 * @brief Immutable configuration snapshots for space-charge algorithms.
 */

#ifndef OPALX_SPACE_CHARGE_CONFIG_H
#define OPALX_SPACE_CHARGE_CONFIG_H

#include "PartBunch/CartesianDomainConfig.h"
#include "SpaceCharge/SpaceChargeCapabilities.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <variant>

namespace opalx::spacecharge {

    /** @brief Poisson backend selected for the Cartesian 3D PIC algorithm. */
    enum class PoissonSolverType : std::uint8_t { None, PeriodicFFT, Open, ConjugateGradient, P3M };

    /** @brief Boundary condition snapshot independent of IPPL boundary types. */
    enum class FieldBoundaryCondition : std::uint8_t { Open, Dirichlet, Periodic };

    /** @brief Free-space Green function used by the open Poisson backend. */
    enum class GreenFunctionType : std::uint8_t { Standard, Integrated };

    /** @brief Bunch attribute used to construct energy or longitudinal bins. */
    enum class BinningVariable : std::uint8_t { VelocityZ, PositionZ, MomentumZ, GammaZ };

    /** @brief Longitudinal-field approximation used by the 2.5D algorithm. */
    enum class FFT2D5LongitudinalFieldMode : std::uint8_t { Open, Cylindrical, Plates, None };

    /**
     * @brief Immutable-by-API source-plane correction configuration.
     *
     * The selected source and its fixed plane belong here. Whether the correction is active on a
     * particular step remains per-call state in SpaceChargeSolveContext.
     */
    class CorrectionConfig {
    public:
        /** @brief Mutable construction parameters resolved during run setup. */
        struct Parameters {
            SpaceChargeCorrectionType kind = SpaceChargeCorrectionType::None;
            double planeZ                  = 0.0;
            std::size_t planeDumpFrequency = 0;
            std::size_t maximumSteps       = 0;
        };

        CorrectionConfig() = default;
        explicit CorrectionConfig(Parameters parameters);

        [[nodiscard]] SpaceChargeCorrectionType kind() const { return parameters_m.kind; }
        [[nodiscard]] bool enabled() const {
            return parameters_m.kind != SpaceChargeCorrectionType::None;
        }
        [[nodiscard]] double planeZ() const { return parameters_m.planeZ; }
        [[nodiscard]] std::size_t planeDumpFrequency() const {
            return parameters_m.planeDumpFrequency;
        }
        [[nodiscard]] std::size_t maximumSteps() const { return parameters_m.maximumSteps; }

    private:
        Parameters parameters_m;
    };

    /**
     * @brief Immutable-by-API binning configuration used by a concrete algorithm.
     *
     * SpaceChargeConfigBuilder fills Parameters while reading parser objects. Runtime algorithms
     * retain only the resulting BinningConfig and never retain parser pointers.
     */
    class BinningConfig {
    public:
        /** @brief Mutable parser-independent construction parameters. */
        struct Parameters {
            std::string name;
            std::size_t maximumBins   = 128;
            double desiredWidth       = 0.1;
            double alpha              = 1.0;
            double beta               = 1.5;
            BinningVariable parameter = BinningVariable::VelocityZ;
            bool adaptive             = true;
            std::string dumpFile;
            std::size_t dumpFrequency       = 1;
            std::size_t tablePrintFrequency = 10;
        };

        explicit BinningConfig(Parameters parameters);

        [[nodiscard]] const std::string& name() const { return parameters_m.name; }
        [[nodiscard]] std::size_t maximumBins() const { return parameters_m.maximumBins; }
        [[nodiscard]] double desiredWidth() const { return parameters_m.desiredWidth; }
        [[nodiscard]] double alpha() const { return parameters_m.alpha; }
        [[nodiscard]] double beta() const { return parameters_m.beta; }
        [[nodiscard]] BinningVariable parameter() const { return parameters_m.parameter; }
        [[nodiscard]] bool adaptive() const { return parameters_m.adaptive; }
        [[nodiscard]] const std::string& dumpFile() const { return parameters_m.dumpFile; }
        [[nodiscard]] std::size_t dumpFrequency() const { return parameters_m.dumpFrequency; }
        [[nodiscard]] std::size_t tablePrintFrequency() const {
            return parameters_m.tablePrintFrequency;
        }

    private:
        Parameters parameters_m;
    };

    /**
     * @brief Immutable-by-API configuration of the complete Cartesian 3D PIC algorithm.
     *
     * This snapshot contains only construction-time values. Per-step correction and beam state
     * belong in SpaceChargeSolveContext.
     */
    class CartesianPICConfig {
    public:
        /** @brief Mutable construction parameters for Cartesian PIC. */
        struct Parameters {
            PoissonSolverType backend = PoissonSolverType::None;
            std::array<std::size_t, 3> meshSize{8, 8, 8};
            std::array<bool, 3> parallelDimensions{true, true, true};
            std::optional<std::array<bool, 3>> layoutRebuildParallelDimensions;
            std::array<FieldBoundaryCondition, 3> boundaryConditions{
                    FieldBoundaryCondition::Open, FieldBoundaryCondition::Open,
                    FieldBoundaryCondition::Open};
            GreenFunctionType greenFunction   = GreenFunctionType::Integrated;
            double p3mCutoff                  = 0.0;
            double boundingBoxIncreasePercent = 2.0;
            std::optional<BinningConfig> binning;
            std::size_t repartitionFrequency = 0;
            double loadBalancingThreshold    = 0.05;
            CorrectionConfig correction;
        };

        explicit CartesianPICConfig(Parameters parameters);

        [[nodiscard]] PoissonSolverType backend() const { return parameters_m.backend; }
        [[nodiscard]] const std::array<std::size_t, 3>& meshSize() const {
            return parameters_m.meshSize;
        }
        [[nodiscard]] const std::array<bool, 3>& parallelDimensions() const {
            return parameters_m.parallelDimensions;
        }
        /**
         * @brief Decomposition used when a correction changes the global mesh extent.
         *
         * Emitted FlatTop beams preserve their legacy longitudinal-only decomposition during
         * image-domain resize without allowing the sampler to mutate live field storage.
         */
        [[nodiscard]] const std::array<bool, 3>& layoutRebuildParallelDimensions() const {
            return parameters_m.layoutRebuildParallelDimensions.has_value()
                           ? *parameters_m.layoutRebuildParallelDimensions
                           : parameters_m.parallelDimensions;
        }
        [[nodiscard]] const std::array<FieldBoundaryCondition, 3>& boundaryConditions() const {
            return parameters_m.boundaryConditions;
        }
        [[nodiscard]] GreenFunctionType greenFunction() const { return parameters_m.greenFunction; }
        [[nodiscard]] double p3mCutoff() const { return parameters_m.p3mCutoff; }
        [[nodiscard]] double boundingBoxIncreasePercent() const {
            return parameters_m.boundingBoxIncreasePercent;
        }
        [[nodiscard]] const std::optional<BinningConfig>& binning() const {
            return parameters_m.binning;
        }
        [[nodiscard]] std::size_t repartitionFrequency() const {
            return parameters_m.repartitionFrequency;
        }
        [[nodiscard]] double loadBalancingThreshold() const {
            return parameters_m.loadBalancingThreshold;
        }
        [[nodiscard]] const CorrectionConfig& correction() const { return parameters_m.correction; }

    private:
        Parameters parameters_m;
    };

    /** @brief Immutable configuration of the reference-path 2.5D PIC algorithm. */
    class FFT2D5Config {
    public:
        /** @brief Mutable construction parameters for the independent FFT2D5 algorithm. */
        struct Parameters {
            std::array<std::size_t, 3> meshSize{8, 8, 8};
            FFT2D5LongitudinalFieldMode longitudinalFieldMode = FFT2D5LongitudinalFieldMode::Open;
            double pipeSizeX                                  = 1.0;
            double pipeSizeY                                  = 1.0;
            double beamRadius                                 = 1.0;
            bool closedRing                                   = false;
            bool scatterLongitudinally                        = true;
            std::string referencePathFile;
        };

        explicit FFT2D5Config(Parameters parameters);

        [[nodiscard]] const std::array<std::size_t, 3>& meshSize() const {
            return parameters_m.meshSize;
        }
        [[nodiscard]] FFT2D5LongitudinalFieldMode longitudinalFieldMode() const {
            return parameters_m.longitudinalFieldMode;
        }
        [[nodiscard]] double pipeSizeX() const { return parameters_m.pipeSizeX; }
        [[nodiscard]] double pipeSizeY() const { return parameters_m.pipeSizeY; }
        [[nodiscard]] double beamRadius() const { return parameters_m.beamRadius; }
        [[nodiscard]] bool closedRing() const { return parameters_m.closedRing; }
        [[nodiscard]] bool scatterLongitudinally() const {
            return parameters_m.scatterLongitudinally;
        }
        [[nodiscard]] const std::string& referencePathFile() const {
            return parameters_m.referencePathFile;
        }

    private:
        Parameters parameters_m;
    };

    using SpaceChargeAlgorithmConfig = std::variant<CartesianPICConfig, FFT2D5Config>;

    /**
     * @brief Run-lifetime configuration owned by SpaceChargeSolver.
     *
     * Add a new alternative to SpaceChargeAlgorithmConfig when a future independent algorithm is
     * integrated. Callers and SpaceChargeSolveContext do not need to change for that addition.
     */
    class SpaceChargeConfig {
    public:
        SpaceChargeConfig(
                SpaceChargeAlgorithmConfig algorithmConfig,
                CartesianDomainConfig3D cartesianDomainConfig);

        [[nodiscard]] SpaceChargeAlgorithmType algorithmType() const;
        [[nodiscard]] const CartesianDomainConfig3D& cartesianDomainConfig() const {
            return cartesianDomainConfig_m;
        }
        [[nodiscard]] const SpaceChargeAlgorithmConfig& algorithmConfig() const {
            return algorithmConfig_m;
        }

        template <typename Config>
        [[nodiscard]] const Config& get() const {
            return std::get<Config>(algorithmConfig_m);
        }

    private:
        SpaceChargeAlgorithmConfig algorithmConfig_m;
        CartesianDomainConfig3D cartesianDomainConfig_m;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_CONFIG_H
