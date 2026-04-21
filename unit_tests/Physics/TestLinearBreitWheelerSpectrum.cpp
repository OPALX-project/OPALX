/**
 * ile TestLinearBreitWheelerSpectrum.cpp
 * rief CAIN-backed regression tests for linear Breit-Wheeler spectra.
 *
 * The Breit-Wheeler validation path is staged in two layers:
 *
 * - a fixed-geometry, unpolarized benchmark against stored CAIN histograms,
 * - a finite incoming-photon-beam benchmark that folds the same kernel over a
 *   Gaussian photon-beam divergence before comparing to CAIN.
 *
 * The tests use normalized histogram area, mean, and L1-distance tolerances
 * instead of pointwise equality so they remain stable under expected Monte
 * Carlo noise while still detecting shape regressions in the generated energy,
 * angle, and joint `E`-`theta` spectra.
 */

#include "LinearBreitWheelerBenchmarkCommon.h"
#include "Utilities/Options.h"
#include "gtest/gtest.h"

#include <cstddef>
#include <filesystem>
#include <stdexcept>
#include <vector>

namespace {
std::filesystem::path referenceElectronEnergyPath() {
    return std::filesystem::path(OPALX_TEST_SOURCE_DIR) / "data"
        / "cain_linear_breit_wheeler_head_on_electron_energy.csv";
}

std::filesystem::path referencePositronEnergyPath() {
    return std::filesystem::path(OPALX_TEST_SOURCE_DIR) / "data"
        / "cain_linear_breit_wheeler_head_on_positron_energy.csv";
}

std::filesystem::path referenceElectronThetaPath() {
    return std::filesystem::path(OPALX_TEST_SOURCE_DIR) / "data"
        / "cain_linear_breit_wheeler_head_on_electron_theta.csv";
}

std::filesystem::path referencePositronThetaPath() {
    return std::filesystem::path(OPALX_TEST_SOURCE_DIR) / "data"
        / "cain_linear_breit_wheeler_head_on_positron_theta.csv";
}

std::filesystem::path referenceElectronJointPath() {
    return std::filesystem::path(OPALX_TEST_SOURCE_DIR) / "data"
        / "cain_linear_breit_wheeler_head_on_electron_joint.csv";
}

std::filesystem::path referencePositronJointPath() {
    return std::filesystem::path(OPALX_TEST_SOURCE_DIR) / "data"
        / "cain_linear_breit_wheeler_head_on_positron_joint.csv";
}

std::filesystem::path referenceFinitePhotonBeamElectronEnergyPath() {
    return std::filesystem::path(OPALX_TEST_SOURCE_DIR) / "data"
        / "cain_linear_breit_wheeler_finite_photon_beam_electron_energy.csv";
}

std::filesystem::path referenceFinitePhotonBeamPositronEnergyPath() {
    return std::filesystem::path(OPALX_TEST_SOURCE_DIR) / "data"
        / "cain_linear_breit_wheeler_finite_photon_beam_positron_energy.csv";
}

std::filesystem::path referenceFinitePhotonBeamElectronThetaPath() {
    return std::filesystem::path(OPALX_TEST_SOURCE_DIR) / "data"
        / "cain_linear_breit_wheeler_finite_photon_beam_electron_theta.csv";
}

std::filesystem::path referenceFinitePhotonBeamPositronThetaPath() {
    return std::filesystem::path(OPALX_TEST_SOURCE_DIR) / "data"
        / "cain_linear_breit_wheeler_finite_photon_beam_positron_theta.csv";
}

std::filesystem::path referenceFinitePhotonBeamElectronJointPath() {
    return std::filesystem::path(OPALX_TEST_SOURCE_DIR) / "data"
        / "cain_linear_breit_wheeler_finite_photon_beam_electron_joint.csv";
}

std::filesystem::path referenceFinitePhotonBeamPositronJointPath() {
    return std::filesystem::path(OPALX_TEST_SOURCE_DIR) / "data"
        / "cain_linear_breit_wheeler_finite_photon_beam_positron_joint.csv";
}

std::filesystem::path referenceFinitePhotonBeamEnergySpreadElectronEnergyPath() {
    return std::filesystem::path(OPALX_TEST_SOURCE_DIR) / "data"
        / "cain_linear_breit_wheeler_finite_photon_beam_energy_spread_electron_energy.csv";
}

std::filesystem::path referenceFinitePhotonBeamEnergySpreadPositronEnergyPath() {
    return std::filesystem::path(OPALX_TEST_SOURCE_DIR) / "data"
        / "cain_linear_breit_wheeler_finite_photon_beam_energy_spread_positron_energy.csv";
}

std::filesystem::path referenceFinitePhotonBeamEnergySpreadElectronThetaPath() {
    return std::filesystem::path(OPALX_TEST_SOURCE_DIR) / "data"
        / "cain_linear_breit_wheeler_finite_photon_beam_energy_spread_electron_theta.csv";
}

std::filesystem::path referenceFinitePhotonBeamEnergySpreadPositronThetaPath() {
    return std::filesystem::path(OPALX_TEST_SOURCE_DIR) / "data"
        / "cain_linear_breit_wheeler_finite_photon_beam_energy_spread_positron_theta.csv";
}

std::filesystem::path referenceFinitePhotonBeamEnergySpreadElectronJointPath() {
    return std::filesystem::path(OPALX_TEST_SOURCE_DIR) / "data"
        / "cain_linear_breit_wheeler_finite_photon_beam_energy_spread_electron_joint.csv";
}

std::filesystem::path referenceFinitePhotonBeamEnergySpreadPositronJointPath() {
    return std::filesystem::path(OPALX_TEST_SOURCE_DIR) / "data"
        / "cain_linear_breit_wheeler_finite_photon_beam_energy_spread_positron_joint.csv";
}

std::filesystem::path referenceOverlapElectronEnergyPath() {
    return std::filesystem::path(OPALX_TEST_SOURCE_DIR) / "data"
        / "cain_linear_breit_wheeler_overlap_electron_energy.csv";
}

std::filesystem::path referenceOverlapPositronEnergyPath() {
    return std::filesystem::path(OPALX_TEST_SOURCE_DIR) / "data"
        / "cain_linear_breit_wheeler_overlap_positron_energy.csv";
}

std::filesystem::path referenceOverlapElectronThetaPath() {
    return std::filesystem::path(OPALX_TEST_SOURCE_DIR) / "data"
        / "cain_linear_breit_wheeler_overlap_electron_theta.csv";
}

std::filesystem::path referenceOverlapPositronThetaPath() {
    return std::filesystem::path(OPALX_TEST_SOURCE_DIR) / "data"
        / "cain_linear_breit_wheeler_overlap_positron_theta.csv";
}

constexpr std::size_t kFastSampleCount = 10000;
constexpr std::size_t kFastBinCount = 20;

/**
 * @brief Build the unit-test histogram configuration.
 *
 * The CAIN references are stored at 80 bins. The unit test intentionally uses
 * 20 bins and rebins the references by preserving integrated probability,
 * @f$\rho'_i \Delta x'_i = \sum_j \rho_j \Delta x_j@f$, so the test remains a
 * shape regression while staying below the 30 s unit-test budget.
 */
LinearBreitWheelerBenchmark::HistogramConfig fastHistogramConfig() {
    LinearBreitWheelerBenchmark::HistogramConfig config;
    config.bins = kFastBinCount;
    return config;
}

/**
 * @brief Build the finite-photon-beam unit-test histogram configuration.
 */
LinearBreitWheelerBenchmark::FinitePhotonBeamConfig fastFinitePhotonBeamConfig() {
    LinearBreitWheelerBenchmark::FinitePhotonBeamConfig config;
    config.bins = kFastBinCount;
    return config;
}

/**
 * @brief Build the fixed-geometry joint unit-test histogram configuration.
 */
LinearBreitWheelerBenchmark::JointHistogramConfig fastJointHistogramConfig() {
    LinearBreitWheelerBenchmark::JointHistogramConfig config;
    config.energyBins = kFastBinCount;
    config.thetaBins = kFastBinCount;
    return config;
}

/**
 * @brief Build the finite-photon-beam joint unit-test histogram configuration.
 */
LinearBreitWheelerBenchmark::FinitePhotonBeamJointConfig fastFinitePhotonBeamJointConfig() {
    LinearBreitWheelerBenchmark::FinitePhotonBeamJointConfig config;
    config.energyBins = kFastBinCount;
    config.thetaBins = kFastBinCount;
    return config;
}

/**
 * @brief Coarsen a one-dimensional reference histogram while preserving area.
 */
LinearBreitWheelerBenchmark::Histogram coarsenHistogram(
    const LinearBreitWheelerBenchmark::Histogram& source,
    std::size_t targetBins) {
    if (targetBins == 0 || source.centers.size() % targetBins != 0) {
        throw std::runtime_error("LinearBreitWheelerSpectrum test: incompatible 1D rebinning.");
    }

    const std::size_t groupSize = source.centers.size() / targetBins;
    LinearBreitWheelerBenchmark::Histogram result;
    result.centers.resize(targetBins);
    result.density.assign(targetBins, 0.0);
    result.counts.assign(targetBins, 0.0);
    result.binWidth = source.binWidth * static_cast<double>(groupSize);
    result.totalWeight = source.totalWeight;

    for (std::size_t i = 0; i < targetBins; ++i) {
        const std::size_t begin = i * groupSize;
        result.centers[i] = 0.5 * (source.centers[begin] + source.centers[begin + groupSize - 1]);
        for (std::size_t j = 0; j < groupSize; ++j) {
            const std::size_t sourceBin = begin + j;
            result.density[i] += source.density[sourceBin] * source.binWidth / result.binWidth;
            result.counts[i] += source.counts[sourceBin];
        }
    }
    return result;
}

/**
 * @brief Coarsen a joint @f$(E,\theta)@f$ reference histogram while preserving area.
 */
LinearBreitWheelerBenchmark::JointHistogram coarsenJointHistogram(
    const LinearBreitWheelerBenchmark::JointHistogram& source,
    std::size_t targetEnergyBins,
    std::size_t targetThetaBins) {
    if (targetEnergyBins == 0 || targetThetaBins == 0
        || source.energyCentersGeV.size() % targetEnergyBins != 0
        || source.thetaCentersRad.size() % targetThetaBins != 0) {
        throw std::runtime_error("LinearBreitWheelerSpectrum test: incompatible joint rebinning.");
    }

    const std::size_t energyGroupSize = source.energyCentersGeV.size() / targetEnergyBins;
    const std::size_t thetaGroupSize = source.thetaCentersRad.size() / targetThetaBins;
    LinearBreitWheelerBenchmark::JointHistogram result;
    result.energyCentersGeV.resize(targetEnergyBins);
    result.thetaCentersRad.resize(targetThetaBins);
    result.densityPerGeVRad.assign(targetEnergyBins,
                                   std::vector<double>(targetThetaBins, 0.0));
    result.counts.assign(targetEnergyBins,
                         std::vector<double>(targetThetaBins, 0.0));
    result.energyBinWidthGeV = source.energyBinWidthGeV * static_cast<double>(energyGroupSize);
    result.thetaBinWidthRad = source.thetaBinWidthRad * static_cast<double>(thetaGroupSize);
    result.totalWeight = source.totalWeight;

    for (std::size_t i = 0; i < targetEnergyBins; ++i) {
        const std::size_t energyBegin = i * energyGroupSize;
        result.energyCentersGeV[i] = 0.5 * (source.energyCentersGeV[energyBegin]
                                            + source.energyCentersGeV[energyBegin + energyGroupSize - 1]);
    }
    for (std::size_t j = 0; j < targetThetaBins; ++j) {
        const std::size_t thetaBegin = j * thetaGroupSize;
        result.thetaCentersRad[j] = 0.5 * (source.thetaCentersRad[thetaBegin]
                                           + source.thetaCentersRad[thetaBegin + thetaGroupSize - 1]);
    }

    const double sourceCellArea = source.energyBinWidthGeV * source.thetaBinWidthRad;
    const double resultCellArea = result.energyBinWidthGeV * result.thetaBinWidthRad;
    for (std::size_t i = 0; i < targetEnergyBins; ++i) {
        for (std::size_t j = 0; j < targetThetaBins; ++j) {
            for (std::size_t di = 0; di < energyGroupSize; ++di) {
                for (std::size_t dj = 0; dj < thetaGroupSize; ++dj) {
                    const std::size_t sourceEnergyBin = i * energyGroupSize + di;
                    const std::size_t sourceThetaBin = j * thetaGroupSize + dj;
                    result.densityPerGeVRad[i][j] += source.densityPerGeVRad[sourceEnergyBin][sourceThetaBin]
                        * sourceCellArea / resultCellArea;
                    result.counts[i][j] += source.counts[sourceEnergyBin][sourceThetaBin];
                }
            }
        }
    }
    return result;
}

LinearBreitWheelerBenchmark::Histogram readLowResolutionHistogram(
    const std::filesystem::path& inputPath) {
    return coarsenHistogram(LinearBreitWheelerBenchmark::readHistogramCSV(inputPath),
                            kFastBinCount);
}

LinearBreitWheelerBenchmark::JointHistogram readLowResolutionJointHistogram(
    const std::filesystem::path& inputPath) {
    return coarsenJointHistogram(LinearBreitWheelerBenchmark::readJointHistogramCSV(inputPath),
                                 kFastBinCount,
                                 kFastBinCount);
}
}  // namespace

TEST(TestLinearBreitWheelerSpectrum, ElectronEnergySpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed = 13579;

    auto config = fastHistogramConfig();
    config.minValue = 0.0;
    config.maxValue = 0.5;
    const auto opalx = LinearBreitWheelerBenchmark::sampleHistogram(
        config,
        LinearBreitWheelerBenchmark::FinalState::Electron,
        LinearBreitWheelerBenchmark::Observable::Energy,
        kFastSampleCount);
    const auto cain = readLowResolutionHistogram(referenceElectronEnergyPath());

    ASSERT_EQ(opalx.centers.size(), cain.centers.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(opalx), 1.0, 1.5e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(cain), 1.0, 5.0e-4);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramMean(opalx),
                LinearBreitWheelerBenchmark::histogramMean(cain),
                LinearBreitWheelerBenchmark::histogramMean(cain) * 5.0e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::histogramL1Distance(opalx, cain), 0.12);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum, PositronEnergySpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed = 13579;

    auto config = fastHistogramConfig();
    config.minValue = 0.0;
    config.maxValue = 0.5;
    const auto opalx = LinearBreitWheelerBenchmark::sampleHistogram(
        config,
        LinearBreitWheelerBenchmark::FinalState::Positron,
        LinearBreitWheelerBenchmark::Observable::Energy,
        kFastSampleCount);
    const auto cain = readLowResolutionHistogram(referencePositronEnergyPath());

    ASSERT_EQ(opalx.centers.size(), cain.centers.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(opalx), 1.0, 1.5e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(cain), 1.0, 5.0e-4);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramMean(opalx),
                LinearBreitWheelerBenchmark::histogramMean(cain),
                LinearBreitWheelerBenchmark::histogramMean(cain) * 5.0e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::histogramL1Distance(opalx, cain), 0.12);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum, ElectronAngularSpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed = 13579;

    auto config = fastHistogramConfig();
    config.minValue = 0.0;
    config.maxValue = 0.0045;
    const auto opalx = LinearBreitWheelerBenchmark::sampleHistogram(
        config,
        LinearBreitWheelerBenchmark::FinalState::Electron,
        LinearBreitWheelerBenchmark::Observable::Theta,
        kFastSampleCount);
    const auto cain = readLowResolutionHistogram(referenceElectronThetaPath());

    ASSERT_EQ(opalx.centers.size(), cain.centers.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(opalx), 1.0, 2.0e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(cain), 1.0, 5.0e-4);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramMean(opalx),
                LinearBreitWheelerBenchmark::histogramMean(cain),
                LinearBreitWheelerBenchmark::histogramMean(cain) * 7.0e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::histogramL1Distance(opalx, cain), 0.15);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum, PositronAngularSpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed = 13579;

    auto config = fastHistogramConfig();
    config.minValue = 0.0;
    config.maxValue = 0.0045;
    const auto opalx = LinearBreitWheelerBenchmark::sampleHistogram(
        config,
        LinearBreitWheelerBenchmark::FinalState::Positron,
        LinearBreitWheelerBenchmark::Observable::Theta,
        kFastSampleCount);
    const auto cain = readLowResolutionHistogram(referencePositronThetaPath());

    ASSERT_EQ(opalx.centers.size(), cain.centers.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(opalx), 1.0, 2.0e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(cain), 1.0, 5.0e-4);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramMean(opalx),
                LinearBreitWheelerBenchmark::histogramMean(cain),
                LinearBreitWheelerBenchmark::histogramMean(cain) * 7.0e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::histogramL1Distance(opalx, cain), 0.15);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum, ElectronJointSpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed = 13579;

    auto config = fastJointHistogramConfig();
    const auto opalx = LinearBreitWheelerBenchmark::sampleJointHistogram(
        config,
        LinearBreitWheelerBenchmark::FinalState::Electron,
        kFastSampleCount);
    const auto cain = readLowResolutionJointHistogram(referenceElectronJointPath());

    ASSERT_EQ(opalx.energyCentersGeV.size(), cain.energyCentersGeV.size());
    ASSERT_EQ(opalx.thetaCentersRad.size(), cain.thetaCentersRad.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramArea(opalx), 1.0, 2.0e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramArea(cain), 1.0, 2.0e-3);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(opalx),
                LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(cain),
                LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(cain) * 5.0e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(opalx),
                LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(cain),
                LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(cain) * 7.0e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::jointHistogramL1Distance(opalx, cain), 0.16);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum, PositronJointSpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed = 13579;

    auto config = fastJointHistogramConfig();
    const auto opalx = LinearBreitWheelerBenchmark::sampleJointHistogram(
        config,
        LinearBreitWheelerBenchmark::FinalState::Positron,
        kFastSampleCount);
    const auto cain = readLowResolutionJointHistogram(referencePositronJointPath());

    ASSERT_EQ(opalx.energyCentersGeV.size(), cain.energyCentersGeV.size());
    ASSERT_EQ(opalx.thetaCentersRad.size(), cain.thetaCentersRad.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramArea(opalx), 1.0, 2.0e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramArea(cain), 1.0, 2.0e-3);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(opalx),
                LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(cain),
                LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(cain) * 5.0e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(opalx),
                LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(cain),
                LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(cain) * 7.0e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::jointHistogramL1Distance(opalx, cain), 0.16);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum, FinitePhotonBeamElectronEnergySpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed = 13579;

    auto config = fastFinitePhotonBeamConfig();
    config.minValue = 0.0;
    config.maxValue = 0.5;
    config.sigmaThetaXRad = 1.0e-3;
    config.sigmaThetaYRad = 1.0e-3;
    const auto opalx = LinearBreitWheelerBenchmark::sampleFinitePhotonBeamHistogram(
        config,
        LinearBreitWheelerBenchmark::FinalState::Electron,
        LinearBreitWheelerBenchmark::Observable::Energy,
        kFastSampleCount);
    const auto cain = readLowResolutionHistogram(referenceFinitePhotonBeamElectronEnergyPath());

    ASSERT_EQ(opalx.centers.size(), cain.centers.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(opalx), 1.0, 2.0e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(cain), 1.0, 5.0e-4);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramMean(opalx),
                LinearBreitWheelerBenchmark::histogramMean(cain),
                LinearBreitWheelerBenchmark::histogramMean(cain) * 8.0e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::histogramL1Distance(opalx, cain), 0.18);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum, FinitePhotonBeamPositronEnergySpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed = 13579;

    auto config = fastFinitePhotonBeamConfig();
    config.minValue = 0.0;
    config.maxValue = 0.5;
    config.sigmaThetaXRad = 1.0e-3;
    config.sigmaThetaYRad = 1.0e-3;
    const auto opalx = LinearBreitWheelerBenchmark::sampleFinitePhotonBeamHistogram(
        config,
        LinearBreitWheelerBenchmark::FinalState::Positron,
        LinearBreitWheelerBenchmark::Observable::Energy,
        kFastSampleCount);
    const auto cain = readLowResolutionHistogram(referenceFinitePhotonBeamPositronEnergyPath());

    ASSERT_EQ(opalx.centers.size(), cain.centers.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(opalx), 1.0, 2.0e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(cain), 1.0, 5.0e-4);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramMean(opalx),
                LinearBreitWheelerBenchmark::histogramMean(cain),
                LinearBreitWheelerBenchmark::histogramMean(cain) * 8.0e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::histogramL1Distance(opalx, cain), 0.18);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum, FinitePhotonBeamElectronAngularSpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed = 13579;

    auto config = fastFinitePhotonBeamConfig();
    config.minValue = 0.0;
    config.maxValue = 0.0060;
    config.sigmaThetaXRad = 1.0e-3;
    config.sigmaThetaYRad = 1.0e-3;
    const auto opalx = LinearBreitWheelerBenchmark::sampleFinitePhotonBeamHistogram(
        config,
        LinearBreitWheelerBenchmark::FinalState::Electron,
        LinearBreitWheelerBenchmark::Observable::Theta,
        kFastSampleCount);
    const auto cain = readLowResolutionHistogram(referenceFinitePhotonBeamElectronThetaPath());

    ASSERT_EQ(opalx.centers.size(), cain.centers.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(opalx), 1.0, 2.5e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(cain), 1.0, 4.0e-3);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramMean(opalx),
                LinearBreitWheelerBenchmark::histogramMean(cain),
                LinearBreitWheelerBenchmark::histogramMean(cain) * 8.0e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::histogramL1Distance(opalx, cain), 0.20);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum, FinitePhotonBeamPositronAngularSpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed = 13579;

    auto config = fastFinitePhotonBeamConfig();
    config.minValue = 0.0;
    config.maxValue = 0.0060;
    config.sigmaThetaXRad = 1.0e-3;
    config.sigmaThetaYRad = 1.0e-3;
    const auto opalx = LinearBreitWheelerBenchmark::sampleFinitePhotonBeamHistogram(
        config,
        LinearBreitWheelerBenchmark::FinalState::Positron,
        LinearBreitWheelerBenchmark::Observable::Theta,
        kFastSampleCount);
    const auto cain = readLowResolutionHistogram(referenceFinitePhotonBeamPositronThetaPath());

    ASSERT_EQ(opalx.centers.size(), cain.centers.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(opalx), 1.0, 2.5e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(cain), 1.0, 4.0e-3);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramMean(opalx),
                LinearBreitWheelerBenchmark::histogramMean(cain),
                LinearBreitWheelerBenchmark::histogramMean(cain) * 8.0e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::histogramL1Distance(opalx, cain), 0.20);

    Options::seed = previousSeed;
}


TEST(TestLinearBreitWheelerSpectrum, FinitePhotonBeamElectronJointSpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed = 13579;

    auto config = fastFinitePhotonBeamJointConfig();
    config.energyMinGeV = 0.0;
    config.energyMaxGeV = 0.5;
    config.thetaMinRad = 0.0;
    config.thetaMaxRad = 0.0060;
    config.sigmaThetaXRad = 1.0e-3;
    config.sigmaThetaYRad = 1.0e-3;
    const auto opalx = LinearBreitWheelerBenchmark::sampleFinitePhotonBeamJointHistogram(
        config,
        LinearBreitWheelerBenchmark::FinalState::Electron,
        kFastSampleCount);
    const auto cain = readLowResolutionJointHistogram(referenceFinitePhotonBeamElectronJointPath());

    ASSERT_EQ(opalx.energyCentersGeV.size(), cain.energyCentersGeV.size());
    ASSERT_EQ(opalx.thetaCentersRad.size(), cain.thetaCentersRad.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramArea(opalx), 1.0, 2.5e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramArea(cain), 1.0, 4.0e-3);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(opalx),
                LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(cain),
                LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(cain) * 8.0e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(opalx),
                LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(cain),
                LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(cain) * 8.5e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::jointHistogramL1Distance(opalx, cain), 0.46);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum, FinitePhotonBeamPositronJointSpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed = 13579;

    auto config = fastFinitePhotonBeamJointConfig();
    config.energyMinGeV = 0.0;
    config.energyMaxGeV = 0.5;
    config.thetaMinRad = 0.0;
    config.thetaMaxRad = 0.0060;
    config.sigmaThetaXRad = 1.0e-3;
    config.sigmaThetaYRad = 1.0e-3;
    const auto opalx = LinearBreitWheelerBenchmark::sampleFinitePhotonBeamJointHistogram(
        config,
        LinearBreitWheelerBenchmark::FinalState::Positron,
        kFastSampleCount);
    const auto cain = readLowResolutionJointHistogram(referenceFinitePhotonBeamPositronJointPath());

    ASSERT_EQ(opalx.energyCentersGeV.size(), cain.energyCentersGeV.size());
    ASSERT_EQ(opalx.thetaCentersRad.size(), cain.thetaCentersRad.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramArea(opalx), 1.0, 2.5e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramArea(cain), 1.0, 4.0e-3);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(opalx),
                LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(cain),
                LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(cain) * 8.0e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(opalx),
                LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(cain),
                LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(cain) * 8.5e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::jointHistogramL1Distance(opalx, cain), 0.46);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum, FinitePhotonBeamEnergySpreadElectronEnergySpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed = 13579;

    auto config = fastFinitePhotonBeamConfig();
    config.minValue = 0.0;
    config.maxValue = 0.5;
    config.sigmaThetaXRad = 1.0e-3;
    config.sigmaThetaYRad = 1.0e-3;
    config.relativeEnergySpread = 1.0e-3;
    const auto opalx = LinearBreitWheelerBenchmark::sampleFinitePhotonBeamHistogram(
        config,
        LinearBreitWheelerBenchmark::FinalState::Electron,
        LinearBreitWheelerBenchmark::Observable::Energy,
        kFastSampleCount);
    const auto cain = readLowResolutionHistogram(referenceFinitePhotonBeamEnergySpreadElectronEnergyPath());

    ASSERT_EQ(opalx.centers.size(), cain.centers.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(opalx), 1.0, 2.0e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(cain), 1.0, 5.0e-4);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramMean(opalx),
                LinearBreitWheelerBenchmark::histogramMean(cain),
                LinearBreitWheelerBenchmark::histogramMean(cain) * 8.0e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::histogramL1Distance(opalx, cain), 0.18);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum, FinitePhotonBeamEnergySpreadPositronEnergySpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed = 13579;

    auto config = fastFinitePhotonBeamConfig();
    config.minValue = 0.0;
    config.maxValue = 0.5;
    config.sigmaThetaXRad = 1.0e-3;
    config.sigmaThetaYRad = 1.0e-3;
    config.relativeEnergySpread = 1.0e-3;
    const auto opalx = LinearBreitWheelerBenchmark::sampleFinitePhotonBeamHistogram(
        config,
        LinearBreitWheelerBenchmark::FinalState::Positron,
        LinearBreitWheelerBenchmark::Observable::Energy,
        kFastSampleCount);
    const auto cain = readLowResolutionHistogram(referenceFinitePhotonBeamEnergySpreadPositronEnergyPath());

    ASSERT_EQ(opalx.centers.size(), cain.centers.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(opalx), 1.0, 2.0e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(cain), 1.0, 5.0e-4);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramMean(opalx),
                LinearBreitWheelerBenchmark::histogramMean(cain),
                LinearBreitWheelerBenchmark::histogramMean(cain) * 8.0e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::histogramL1Distance(opalx, cain), 0.18);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum, FinitePhotonBeamEnergySpreadElectronAngularSpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed = 13579;

    auto config = fastFinitePhotonBeamConfig();
    config.minValue = 0.0;
    config.maxValue = 0.0060;
    config.sigmaThetaXRad = 1.0e-3;
    config.sigmaThetaYRad = 1.0e-3;
    config.relativeEnergySpread = 1.0e-3;
    const auto opalx = LinearBreitWheelerBenchmark::sampleFinitePhotonBeamHistogram(
        config,
        LinearBreitWheelerBenchmark::FinalState::Electron,
        LinearBreitWheelerBenchmark::Observable::Theta,
        kFastSampleCount);
    const auto cain = readLowResolutionHistogram(referenceFinitePhotonBeamEnergySpreadElectronThetaPath());

    ASSERT_EQ(opalx.centers.size(), cain.centers.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(opalx), 1.0, 2.5e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(cain), 1.0, 4.0e-3);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramMean(opalx),
                LinearBreitWheelerBenchmark::histogramMean(cain),
                LinearBreitWheelerBenchmark::histogramMean(cain) * 8.0e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::histogramL1Distance(opalx, cain), 0.20);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum, FinitePhotonBeamEnergySpreadPositronAngularSpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed = 13579;

    auto config = fastFinitePhotonBeamConfig();
    config.minValue = 0.0;
    config.maxValue = 0.0060;
    config.sigmaThetaXRad = 1.0e-3;
    config.sigmaThetaYRad = 1.0e-3;
    config.relativeEnergySpread = 1.0e-3;
    const auto opalx = LinearBreitWheelerBenchmark::sampleFinitePhotonBeamHistogram(
        config,
        LinearBreitWheelerBenchmark::FinalState::Positron,
        LinearBreitWheelerBenchmark::Observable::Theta,
        kFastSampleCount);
    const auto cain = readLowResolutionHistogram(referenceFinitePhotonBeamEnergySpreadPositronThetaPath());

    ASSERT_EQ(opalx.centers.size(), cain.centers.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(opalx), 1.0, 2.5e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(cain), 1.0, 4.0e-3);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramMean(opalx),
                LinearBreitWheelerBenchmark::histogramMean(cain),
                LinearBreitWheelerBenchmark::histogramMean(cain) * 8.0e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::histogramL1Distance(opalx, cain), 0.20);

    Options::seed = previousSeed;
}


TEST(TestLinearBreitWheelerSpectrum, FinitePhotonBeamEnergySpreadElectronJointSpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed = 13579;

    auto config = fastFinitePhotonBeamJointConfig();
    config.energyMinGeV = 0.0;
    config.energyMaxGeV = 0.5;
    config.thetaMinRad = 0.0;
    config.thetaMaxRad = 0.0060;
    config.sigmaThetaXRad = 1.0e-3;
    config.sigmaThetaYRad = 1.0e-3;
    config.relativeEnergySpread = 1.0e-3;
    const auto opalx = LinearBreitWheelerBenchmark::sampleFinitePhotonBeamJointHistogram(
        config,
        LinearBreitWheelerBenchmark::FinalState::Electron,
        kFastSampleCount);
    const auto cain = readLowResolutionJointHistogram(referenceFinitePhotonBeamEnergySpreadElectronJointPath());

    ASSERT_EQ(opalx.energyCentersGeV.size(), cain.energyCentersGeV.size());
    ASSERT_EQ(opalx.thetaCentersRad.size(), cain.thetaCentersRad.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramArea(opalx), 1.0, 2.5e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramArea(cain), 1.0, 4.0e-3);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(opalx),
                LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(cain),
                LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(cain) * 8.5e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(opalx),
                LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(cain),
                LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(cain) * 8.5e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::jointHistogramL1Distance(opalx, cain), 0.46);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum, FinitePhotonBeamEnergySpreadPositronJointSpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed = 13579;

    auto config = fastFinitePhotonBeamJointConfig();
    config.energyMinGeV = 0.0;
    config.energyMaxGeV = 0.5;
    config.thetaMinRad = 0.0;
    config.thetaMaxRad = 0.0060;
    config.sigmaThetaXRad = 1.0e-3;
    config.sigmaThetaYRad = 1.0e-3;
    config.relativeEnergySpread = 1.0e-3;
    const auto opalx = LinearBreitWheelerBenchmark::sampleFinitePhotonBeamJointHistogram(
        config,
        LinearBreitWheelerBenchmark::FinalState::Positron,
        kFastSampleCount);
    const auto cain = readLowResolutionJointHistogram(referenceFinitePhotonBeamEnergySpreadPositronJointPath());

    ASSERT_EQ(opalx.energyCentersGeV.size(), cain.energyCentersGeV.size());
    ASSERT_EQ(opalx.thetaCentersRad.size(), cain.thetaCentersRad.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramArea(opalx), 1.0, 2.5e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramArea(cain), 1.0, 4.0e-3);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(opalx),
                LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(cain),
                LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(cain) * 8.5e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(opalx),
                LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(cain),
                LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(cain) * 8.5e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::jointHistogramL1Distance(opalx, cain), 0.46);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum, OverlapWeightedElectronEnergySpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed = 13579;

    auto config = fastFinitePhotonBeamConfig();
    config.minValue = 0.0;
    config.maxValue = 0.5;
    config.sigmaThetaXRad = 1.0e-3;
    config.sigmaThetaYRad = 1.0e-3;
    config.sigmaX_m = 1.0e-6;
    config.sigmaY_m = 1.0e-6;
    config.sigmaS_m = 5.0e-12 * Physics::c;
    config.laserRayleighX_m = 1.0e-6;
    config.laserRayleighY_m = 1.0e-6;
    config.laserSigmaT_m = 5.0e-12 * Physics::c;
    config.overlapWeighting = true;
    const auto opalx = LinearBreitWheelerBenchmark::sampleFinitePhotonBeamHistogram(
        config,
        LinearBreitWheelerBenchmark::FinalState::Electron,
        LinearBreitWheelerBenchmark::Observable::Energy,
        kFastSampleCount);
    const auto cain = readLowResolutionHistogram(referenceOverlapElectronEnergyPath());

    ASSERT_EQ(opalx.centers.size(), cain.centers.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(opalx), 1.0, 2.0e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramMean(opalx),
                LinearBreitWheelerBenchmark::histogramMean(cain),
                LinearBreitWheelerBenchmark::histogramMean(cain) * 1.2e-1);
    EXPECT_LT(LinearBreitWheelerBenchmark::histogramL1Distance(opalx, cain), 0.28);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum, OverlapWeightedPositronEnergySpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed = 13579;

    auto config = fastFinitePhotonBeamConfig();
    config.minValue = 0.0;
    config.maxValue = 0.5;
    config.sigmaThetaXRad = 1.0e-3;
    config.sigmaThetaYRad = 1.0e-3;
    config.sigmaX_m = 1.0e-6;
    config.sigmaY_m = 1.0e-6;
    config.sigmaS_m = 5.0e-12 * Physics::c;
    config.laserRayleighX_m = 1.0e-6;
    config.laserRayleighY_m = 1.0e-6;
    config.laserSigmaT_m = 5.0e-12 * Physics::c;
    config.overlapWeighting = true;
    const auto opalx = LinearBreitWheelerBenchmark::sampleFinitePhotonBeamHistogram(
        config,
        LinearBreitWheelerBenchmark::FinalState::Positron,
        LinearBreitWheelerBenchmark::Observable::Energy,
        kFastSampleCount);
    const auto cain = readLowResolutionHistogram(referenceOverlapPositronEnergyPath());

    ASSERT_EQ(opalx.centers.size(), cain.centers.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(opalx), 1.0, 2.0e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramMean(opalx),
                LinearBreitWheelerBenchmark::histogramMean(cain),
                LinearBreitWheelerBenchmark::histogramMean(cain) * 1.2e-1);
    EXPECT_LT(LinearBreitWheelerBenchmark::histogramL1Distance(opalx, cain), 0.28);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum, OverlapWeightedElectronAngularSpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed = 13579;

    auto config = fastFinitePhotonBeamConfig();
    config.minValue = 0.0;
    config.maxValue = 0.0060;
    config.sigmaThetaXRad = 1.0e-3;
    config.sigmaThetaYRad = 1.0e-3;
    config.sigmaX_m = 1.0e-6;
    config.sigmaY_m = 1.0e-6;
    config.sigmaS_m = 5.0e-12 * Physics::c;
    config.laserRayleighX_m = 1.0e-6;
    config.laserRayleighY_m = 1.0e-6;
    config.laserSigmaT_m = 5.0e-12 * Physics::c;
    config.overlapWeighting = true;
    const auto opalx = LinearBreitWheelerBenchmark::sampleFinitePhotonBeamHistogram(
        config,
        LinearBreitWheelerBenchmark::FinalState::Electron,
        LinearBreitWheelerBenchmark::Observable::Theta,
        kFastSampleCount);
    const auto cain = readLowResolutionHistogram(referenceOverlapElectronThetaPath());

    ASSERT_EQ(opalx.centers.size(), cain.centers.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(opalx), 1.0, 2.5e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramMean(opalx),
                LinearBreitWheelerBenchmark::histogramMean(cain),
                LinearBreitWheelerBenchmark::histogramMean(cain) * 1.2e-1);
    EXPECT_LT(LinearBreitWheelerBenchmark::histogramL1Distance(opalx, cain), 0.32);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum, OverlapWeightedPositronAngularSpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed = 13579;

    auto config = fastFinitePhotonBeamConfig();
    config.minValue = 0.0;
    config.maxValue = 0.0060;
    config.sigmaThetaXRad = 1.0e-3;
    config.sigmaThetaYRad = 1.0e-3;
    config.sigmaX_m = 1.0e-6;
    config.sigmaY_m = 1.0e-6;
    config.sigmaS_m = 5.0e-12 * Physics::c;
    config.laserRayleighX_m = 1.0e-6;
    config.laserRayleighY_m = 1.0e-6;
    config.laserSigmaT_m = 5.0e-12 * Physics::c;
    config.overlapWeighting = true;
    const auto opalx = LinearBreitWheelerBenchmark::sampleFinitePhotonBeamHistogram(
        config,
        LinearBreitWheelerBenchmark::FinalState::Positron,
        LinearBreitWheelerBenchmark::Observable::Theta,
        kFastSampleCount);
    const auto cain = readLowResolutionHistogram(referenceOverlapPositronThetaPath());

    ASSERT_EQ(opalx.centers.size(), cain.centers.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(opalx), 1.0, 2.5e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramMean(opalx),
                LinearBreitWheelerBenchmark::histogramMean(cain),
                LinearBreitWheelerBenchmark::histogramMean(cain) * 1.2e-1);
    EXPECT_LT(LinearBreitWheelerBenchmark::histogramL1Distance(opalx, cain), 0.32);

    Options::seed = previousSeed;
}
