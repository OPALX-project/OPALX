/**
 * @file SpaceChargeConfig.h
 * @brief Parser-independent configuration values for space-charge algorithms.
 */

#ifndef OPALX_SPACE_CHARGE_CONFIG_H
#define OPALX_SPACE_CHARGE_CONFIG_H

#include "PartBunch/CartesianDomainConfig.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <variant>

namespace opalx::spacecharge {

    enum class SpaceChargeCorrectionType : std::uint8_t { None, ImageCharge, ShiftedGreen };
    enum class PoissonSolverType : std::uint8_t { None, PeriodicFFT, Open, ConjugateGradient, P3M };
    /** @brief Boundary of the Poisson domain, independent of source-plane corrections. */
    enum class FieldBoundaryCondition : std::uint8_t { Open, Dirichlet, Periodic };
    enum class GreenFunctionType : std::uint8_t { Standard, Integrated };
    enum class BinningVariable : std::uint8_t { VelocityZ, PositionZ, MomentumZ, GammaZ };
    enum class FFT2D5LongitudinalFieldMode : std::uint8_t { Open, Cylindrical, Plates, None };

    /** @brief Backend-independent settings for a 3D Poisson solve. */
    struct PoissonSolverConfig {
        PoissonSolverType type          = PoissonSolverType::None;
        GreenFunctionType greenFunction = GreenFunctionType::Integrated;
        double p3mCutoff                = 0.0;  ///< P3M cutoff radius in metres; otherwise zero.
        std::array<FieldBoundaryCondition, 3> boundaryConditions{
                FieldBoundaryCondition::Open, FieldBoundaryCondition::Open,
                FieldBoundaryCondition::Open};
    };

    /** @brief Cartesian mesh extents, MPI decomposition, and particle-bound margin. */
    struct CartesianGridConfig {
        std::array<std::size_t, 3> meshSize{8, 8, 8};         ///< Grid points per axis.
        std::array<bool, 3> decomposition{true, true, true};  ///< Distributed MPI axes.
        double boundingBoxIncreasePercent = 2.0;  ///< Particle-span percentage added per side.
    };

    /** @brief Source-plane field correction; planeZ is in metres in the Cartesian solve frame. */
    struct CorrectionConfig {
        SpaceChargeCorrectionType kind = SpaceChargeCorrectionType::None;
        double planeZ                  = 0.0;
        std::size_t planeDumpFrequency = 0;  ///< Steps between plane dumps; zero disables them.
        std::size_t maximumSteps       = 0;  ///< Number of corrected steps; zero never expires.

        [[nodiscard]] bool enabled() const { return kind != SpaceChargeCorrectionType::None; }
    };

    /** @brief Particle binning and its diagnostics. */
    struct BinningConfig {
        std::string name;                 ///< BINNING definition name used in diagnostics.
        std::size_t maximumBins   = 128;  ///< Initial uniform histogram size.
        double desiredWidth       = 0.1;  ///< Target normalized width in the merge cost.
        double alpha              = 1.0;  ///< Weight penalizing wide bins.
        double beta               = 1.5;  ///< Weight penalizing deviation from desiredWidth.
        BinningVariable parameter = BinningVariable::VelocityZ;
        bool adaptive             = true;  ///< Merge the initial uniform histogram.
        std::string dumpFile;
        std::size_t dumpFrequency       = 1;   ///< Steps between JSON snapshots.
        std::size_t tablePrintFrequency = 10;  ///< Steps between console tables.
    };

    /** @brief Complete runtime configuration for Cartesian PIC. */
    struct CartesianPICConfig {
        CartesianGridConfig grid;
        PoissonSolverType backend = PoissonSolverType::None;
        /** @brief MPI decomposition used after mesh resize; absent uses grid.decomposition. */
        std::optional<std::array<bool, 3>> layoutRebuildDecomposition;
        std::array<FieldBoundaryCondition, 3> boundaryConditions{
                FieldBoundaryCondition::Open, FieldBoundaryCondition::Open,
                FieldBoundaryCondition::Open};
        GreenFunctionType greenFunction = GreenFunctionType::Integrated;
        double p3mCutoff                = 0.0;  ///< P3M cutoff radius in metres.
        std::optional<BinningConfig> binning;
        std::size_t repartitionFrequency = 0;  ///< Steps between ORB checks; zero disables them.
        /** @brief Per-rank count deviation, divided by global count, that triggers ORB. */
        double loadBalancingThreshold = 0.05;
        CorrectionConfig correction;

        [[nodiscard]] const std::array<bool, 3>& layoutDecomposition() const {
            return layoutRebuildDecomposition.has_value() ? *layoutRebuildDecomposition
                                                          : grid.decomposition;
        }
    };

    struct FFT2D5Config {
        CartesianGridConfig grid{.decomposition = {false, false, false}};
        FFT2D5LongitudinalFieldMode longitudinalFieldMode = FFT2D5LongitudinalFieldMode::Open;
        double pipeSizeX                                  = 1.0;
        double pipeSizeY                                  = 1.0;
        double beamRadius                                 = 1.0;
        bool closedRing                                   = false;
        bool scatterLongitudinally                        = true;
        std::string referencePathFile;
    };

    using SpaceChargeConfig = std::variant<CartesianPICConfig, FFT2D5Config>;

    /** @brief Reject unsupported or inconsistent Poisson settings. */
    void validatePoissonSolverConfig(const PoissonSolverConfig& config);
    /** @brief Extract Poisson settings from Cartesian PIC configuration. */
    [[nodiscard]] PoissonSolverConfig makePoissonSolverConfig(const CartesianPICConfig& config);

    /** @brief Reject unsupported or inconsistent space-charge settings. */
    void validateSpaceChargeConfig(const SpaceChargeConfig& config);
    /** @brief Derive the initial PartBunch Cartesian domain settings. */
    [[nodiscard]] CartesianDomainConfig3D makeCartesianDomainConfig(
            const SpaceChargeConfig& config);

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_CONFIG_H
