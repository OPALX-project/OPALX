/**
 * \file TestLinearBreitWheelerSpectrum.cpp
 * \brief CAIN-backed regression tests for linear Breit-Wheeler spectra.
 *
 * The first Breit-Wheeler validation layer mirrors the earlier linear-Compton
 * strategy: keep the physics host-only and compare sampled OPALX spectra
 * against stored CAIN histograms for a fixed-geometry, unpolarized benchmark.
 *
 * The tests use normalized histogram area, mean, and L1-distance tolerances
 * instead of pointwise equality so they remain stable under expected Monte
 * Carlo noise while still detecting shape regressions in the generated energy,
 * angle, and joint `E`-`theta` spectra.
 */

#include "LinearBreitWheelerBenchmarkCommon.h"
#include "Utilities/Options.h"
#include "gtest/gtest.h"

#include <filesystem>

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
}  // namespace

TEST(TestLinearBreitWheelerSpectrum, ElectronEnergySpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed = 13579;

    LinearBreitWheelerBenchmark::HistogramConfig config;
    config.minValue = 0.0;
    config.maxValue = 0.5;
    const auto opalx = LinearBreitWheelerBenchmark::sampleHistogram(
        config,
        LinearBreitWheelerBenchmark::FinalState::Electron,
        LinearBreitWheelerBenchmark::Observable::Energy,
        250000);
    const auto cain = LinearBreitWheelerBenchmark::readHistogramCSV(referenceElectronEnergyPath());

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

    LinearBreitWheelerBenchmark::HistogramConfig config;
    config.minValue = 0.0;
    config.maxValue = 0.5;
    const auto opalx = LinearBreitWheelerBenchmark::sampleHistogram(
        config,
        LinearBreitWheelerBenchmark::FinalState::Positron,
        LinearBreitWheelerBenchmark::Observable::Energy,
        250000);
    const auto cain = LinearBreitWheelerBenchmark::readHistogramCSV(referencePositronEnergyPath());

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

    LinearBreitWheelerBenchmark::HistogramConfig config;
    config.minValue = 0.0;
    config.maxValue = 0.0045;
    const auto opalx = LinearBreitWheelerBenchmark::sampleHistogram(
        config,
        LinearBreitWheelerBenchmark::FinalState::Electron,
        LinearBreitWheelerBenchmark::Observable::Theta,
        250000);
    const auto cain = LinearBreitWheelerBenchmark::readHistogramCSV(referenceElectronThetaPath());

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

    LinearBreitWheelerBenchmark::HistogramConfig config;
    config.minValue = 0.0;
    config.maxValue = 0.0045;
    const auto opalx = LinearBreitWheelerBenchmark::sampleHistogram(
        config,
        LinearBreitWheelerBenchmark::FinalState::Positron,
        LinearBreitWheelerBenchmark::Observable::Theta,
        250000);
    const auto cain = LinearBreitWheelerBenchmark::readHistogramCSV(referencePositronThetaPath());

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

    LinearBreitWheelerBenchmark::JointHistogramConfig config;
    const auto opalx = LinearBreitWheelerBenchmark::sampleJointHistogram(
        config,
        LinearBreitWheelerBenchmark::FinalState::Electron,
        250000);
    const auto cain = LinearBreitWheelerBenchmark::readJointHistogramCSV(referenceElectronJointPath());

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

    LinearBreitWheelerBenchmark::JointHistogramConfig config;
    const auto opalx = LinearBreitWheelerBenchmark::sampleJointHistogram(
        config,
        LinearBreitWheelerBenchmark::FinalState::Positron,
        250000);
    const auto cain = LinearBreitWheelerBenchmark::readJointHistogramCSV(referencePositronJointPath());

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
