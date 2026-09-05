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

    struct PoissonSolverConfig {
        PoissonSolverType type          = PoissonSolverType::None;
        GreenFunctionType greenFunction = GreenFunctionType::Integrated;
        double p3mCutoff                = 0.0;
        std::array<FieldBoundaryCondition, 3> boundaryConditions{
                FieldBoundaryCondition::Open, FieldBoundaryCondition::Open,
                FieldBoundaryCondition::Open};
    };

    struct CartesianGridConfig {
        std::array<std::size_t, 3> meshSize{8, 8, 8};
        std::array<bool, 3> decomposition{true, true, true};
        double boundingBoxIncreasePercent = 2.0;
    };

    /** @brief Source-plane field correction; planeZ is in metres in the Cartesian solve frame. */
    struct CorrectionConfig {
        SpaceChargeCorrectionType kind = SpaceChargeCorrectionType::None;
        double planeZ                  = 0.0;
        std::size_t planeDumpFrequency = 0;
        std::size_t maximumSteps       = 0;

        [[nodiscard]] bool enabled() const { return kind != SpaceChargeCorrectionType::None; }
    };

    struct BinningConfig {
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

    struct CartesianPICConfig {
        CartesianGridConfig grid;
        PoissonSolverType backend = PoissonSolverType::None;
        std::optional<std::array<bool, 3>> layoutRebuildDecomposition;
        std::array<FieldBoundaryCondition, 3> boundaryConditions{
                FieldBoundaryCondition::Open, FieldBoundaryCondition::Open,
                FieldBoundaryCondition::Open};
        GreenFunctionType greenFunction = GreenFunctionType::Integrated;
        double p3mCutoff                = 0.0;
        std::optional<BinningConfig> binning;
        std::size_t repartitionFrequency = 0;
        double loadBalancingThreshold    = 0.05;
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

    void validatePoissonSolverConfig(const PoissonSolverConfig& config);
    [[nodiscard]] PoissonSolverConfig makePoissonSolverConfig(const CartesianPICConfig& config);

    void validateSpaceChargeConfig(const SpaceChargeConfig& config);
    [[nodiscard]] CartesianDomainConfig3D makeCartesianDomainConfig(
            const SpaceChargeConfig& config);

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_CONFIG_H
