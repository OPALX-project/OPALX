/**
 * @file SelfFieldConfig.h
 * @brief Immutable configuration snapshots for self-field algorithms.
 */

#ifndef OPALX_SPACE_CHARGE_SELF_FIELD_CONFIG_H
#define OPALX_SPACE_CHARGE_SELF_FIELD_CONFIG_H

#include "SpaceCharge/SolverCapabilities.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <variant>

namespace opalx::spacecharge {

    /** @brief Poisson backend selected for the Cartesian 3D PIC algorithm. */
    enum class PoissonBackendKind : std::uint8_t { None, FftPeriodic, Open, ConjugateGradient };

    /** @brief Boundary condition snapshot independent of IPPL boundary types. */
    enum class BoundaryConditionKind : std::uint8_t { Open, Dirichlet, Periodic };

    /** @brief Free-space Green function used by the open Poisson backend. */
    enum class GreenFunctionKind : std::uint8_t { Standard, Integrated };

    /** @brief Bunch attribute used to construct energy or longitudinal bins. */
    enum class BinningParameterKind : std::uint8_t { VelocityZ, PositionZ, MomentumZ, GammaZ };

    /** @brief Static source-plane correction values resolved during run setup. */
    struct CorrectionConfigValues {
        CorrectionKind kind            = CorrectionKind::None;
        double planeZ                  = 0.0;
        std::size_t planeDumpFrequency = 0;
        std::size_t maximumSteps       = 0;
    };

    /**
     * @brief Immutable-by-API source-plane correction configuration.
     *
     * The selected source and its fixed plane belong here. Whether the correction is active on a
     * particular step remains per-call state in SolveContext.
     */
    class CorrectionConfig {
    public:
        explicit CorrectionConfig(CorrectionConfigValues values = {});

        [[nodiscard]] CorrectionKind kind() const { return values_m.kind; }
        [[nodiscard]] bool enabled() const { return values_m.kind != CorrectionKind::None; }
        [[nodiscard]] double planeZ() const { return values_m.planeZ; }
        [[nodiscard]] std::size_t planeDumpFrequency() const { return values_m.planeDumpFrequency; }
        [[nodiscard]] std::size_t maximumSteps() const { return values_m.maximumSteps; }

    private:
        CorrectionConfigValues values_m;
    };

    /**
     * @brief Mutable construction values used to create an immutable BinningConfig snapshot.
     *
     * SelfFieldConfigBuilder may fill this record while reading parser objects. Runtime solver
     * objects retain only the resulting BinningConfig and never retain parser pointers.
     */
    struct BinningConfigValues {
        std::string name;
        std::size_t maximumBins        = 128;
        double desiredWidth            = 0.1;
        double alpha                   = 1.0;
        double beta                    = 1.5;
        BinningParameterKind parameter = BinningParameterKind::VelocityZ;
        bool adaptive                  = true;
        std::string dumpFile;
        std::size_t dumpFrequency       = 1;
        std::size_t tablePrintFrequency = 10;
    };

    /** @brief Immutable-by-API binning configuration used by a concrete algorithm. */
    class BinningConfig {
    public:
        explicit BinningConfig(BinningConfigValues values);

        [[nodiscard]] const std::string& name() const { return values_m.name; }
        [[nodiscard]] std::size_t maximumBins() const { return values_m.maximumBins; }
        [[nodiscard]] double desiredWidth() const { return values_m.desiredWidth; }
        [[nodiscard]] double alpha() const { return values_m.alpha; }
        [[nodiscard]] double beta() const { return values_m.beta; }
        [[nodiscard]] BinningParameterKind parameter() const { return values_m.parameter; }
        [[nodiscard]] bool adaptive() const { return values_m.adaptive; }
        [[nodiscard]] const std::string& dumpFile() const { return values_m.dumpFile; }
        [[nodiscard]] std::size_t dumpFrequency() const { return values_m.dumpFrequency; }
        [[nodiscard]] std::size_t tablePrintFrequency() const {
            return values_m.tablePrintFrequency;
        }

    private:
        BinningConfigValues values_m;
    };

    /**
     * @brief Mutable construction values for a Cartesian 3D PIC configuration snapshot.
     */
    struct Pic3DConfigValues {
        PoissonBackendKind backend = PoissonBackendKind::None;
        std::array<std::size_t, 3> meshSize{8, 8, 8};
        std::array<bool, 3> parallelDimensions{true, true, true};
        std::optional<std::array<bool, 3>> layoutRebuildParallelDimensions;
        std::array<BoundaryConditionKind, 3> boundaryConditions{
                BoundaryConditionKind::Open, BoundaryConditionKind::Open,
                BoundaryConditionKind::Open};
        GreenFunctionKind greenFunction   = GreenFunctionKind::Integrated;
        double boundingBoxIncreasePercent = 2.0;
        std::optional<BinningConfig> binning;
        std::size_t repartitionFrequency = 0;
        double loadBalancingThreshold    = 0.05;
        CorrectionConfig correction;
    };

    /**
     * @brief Immutable-by-API configuration of the complete Cartesian 3D PIC algorithm.
     *
     * This snapshot contains only construction-time values. Per-step correction and beam state
     * belong in SolveContext.
     */
    class Pic3DConfig {
    public:
        explicit Pic3DConfig(Pic3DConfigValues values);

        [[nodiscard]] PoissonBackendKind backend() const { return values_m.backend; }
        [[nodiscard]] const std::array<std::size_t, 3>& meshSize() const {
            return values_m.meshSize;
        }
        [[nodiscard]] const std::array<bool, 3>& parallelDimensions() const {
            return values_m.parallelDimensions;
        }
        /**
         * @brief Decomposition used when a correction changes the global mesh extent.
         *
         * Emitted FlatTop beams preserve their legacy longitudinal-only decomposition during
         * image-domain resize without allowing the sampler to mutate the live workspace.
         */
        [[nodiscard]] const std::array<bool, 3>& layoutRebuildParallelDimensions() const {
            return values_m.layoutRebuildParallelDimensions.has_value()
                           ? *values_m.layoutRebuildParallelDimensions
                           : values_m.parallelDimensions;
        }
        [[nodiscard]] const std::array<BoundaryConditionKind, 3>& boundaryConditions() const {
            return values_m.boundaryConditions;
        }
        [[nodiscard]] GreenFunctionKind greenFunction() const { return values_m.greenFunction; }
        [[nodiscard]] double boundingBoxIncreasePercent() const {
            return values_m.boundingBoxIncreasePercent;
        }
        [[nodiscard]] const std::optional<BinningConfig>& binning() const {
            return values_m.binning;
        }
        [[nodiscard]] std::size_t repartitionFrequency() const {
            return values_m.repartitionFrequency;
        }
        [[nodiscard]] double loadBalancingThreshold() const {
            return values_m.loadBalancingThreshold;
        }
        [[nodiscard]] const CorrectionConfig& correction() const { return values_m.correction; }

    private:
        Pic3DConfigValues values_m;
    };

    using SelfFieldAlgorithmConfig = std::variant<Pic3DConfig>;

    /**
     * @brief Run-lifetime configuration owned by SelfFieldSystem.
     *
     * Add a new alternative to SelfFieldAlgorithmConfig when a future independent algorithm is
     * integrated. Callers and SolveContext do not need to change for that addition.
     */
    class SelfFieldConfig {
    public:
        explicit SelfFieldConfig(SelfFieldAlgorithmConfig algorithmConfig);

        [[nodiscard]] SelfFieldAlgorithmKind algorithmKind() const;
        [[nodiscard]] const SelfFieldAlgorithmConfig& algorithmConfig() const {
            return algorithmConfig_m;
        }

        template <typename Config>
        [[nodiscard]] const Config& get() const {
            return std::get<Config>(algorithmConfig_m);
        }

    private:
        SelfFieldAlgorithmConfig algorithmConfig_m;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_SELF_FIELD_CONFIG_H
