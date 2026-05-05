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
}  // namespace

TEST(TestLinearBreitWheelerSpectrum, ElectronEnergySpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed          = 13579;

    LinearBreitWheelerBenchmark::HistogramConfig config;
    config.minValue  = 0.0;
    config.maxValue  = 0.5;
    const auto opalx = LinearBreitWheelerBenchmark::sampleHistogram(
            config, LinearBreitWheelerBenchmark::FinalState::Electron,
            LinearBreitWheelerBenchmark::Observable::Energy, 250000);
    const auto cain = LinearBreitWheelerBenchmark::readHistogramCSV(referenceElectronEnergyPath());

    ASSERT_EQ(opalx.centers.size(), cain.centers.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(opalx), 1.0, 1.5e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(cain), 1.0, 5.0e-4);
    EXPECT_NEAR(
            LinearBreitWheelerBenchmark::histogramMean(opalx),
            LinearBreitWheelerBenchmark::histogramMean(cain),
            LinearBreitWheelerBenchmark::histogramMean(cain) * 5.0e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::histogramL1Distance(opalx, cain), 0.12);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum, PositronEnergySpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed          = 13579;

    LinearBreitWheelerBenchmark::HistogramConfig config;
    config.minValue  = 0.0;
    config.maxValue  = 0.5;
    const auto opalx = LinearBreitWheelerBenchmark::sampleHistogram(
            config, LinearBreitWheelerBenchmark::FinalState::Positron,
            LinearBreitWheelerBenchmark::Observable::Energy, 250000);
    const auto cain = LinearBreitWheelerBenchmark::readHistogramCSV(referencePositronEnergyPath());

    ASSERT_EQ(opalx.centers.size(), cain.centers.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(opalx), 1.0, 1.5e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(cain), 1.0, 5.0e-4);
    EXPECT_NEAR(
            LinearBreitWheelerBenchmark::histogramMean(opalx),
            LinearBreitWheelerBenchmark::histogramMean(cain),
            LinearBreitWheelerBenchmark::histogramMean(cain) * 5.0e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::histogramL1Distance(opalx, cain), 0.12);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum, ElectronAngularSpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed          = 13579;

    LinearBreitWheelerBenchmark::HistogramConfig config;
    config.minValue  = 0.0;
    config.maxValue  = 0.0045;
    const auto opalx = LinearBreitWheelerBenchmark::sampleHistogram(
            config, LinearBreitWheelerBenchmark::FinalState::Electron,
            LinearBreitWheelerBenchmark::Observable::Theta, 250000);
    const auto cain = LinearBreitWheelerBenchmark::readHistogramCSV(referenceElectronThetaPath());

    ASSERT_EQ(opalx.centers.size(), cain.centers.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(opalx), 1.0, 2.0e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(cain), 1.0, 5.0e-4);
    EXPECT_NEAR(
            LinearBreitWheelerBenchmark::histogramMean(opalx),
            LinearBreitWheelerBenchmark::histogramMean(cain),
            LinearBreitWheelerBenchmark::histogramMean(cain) * 7.0e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::histogramL1Distance(opalx, cain), 0.15);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum, PositronAngularSpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed          = 13579;

    LinearBreitWheelerBenchmark::HistogramConfig config;
    config.minValue  = 0.0;
    config.maxValue  = 0.0045;
    const auto opalx = LinearBreitWheelerBenchmark::sampleHistogram(
            config, LinearBreitWheelerBenchmark::FinalState::Positron,
            LinearBreitWheelerBenchmark::Observable::Theta, 250000);
    const auto cain = LinearBreitWheelerBenchmark::readHistogramCSV(referencePositronThetaPath());

    ASSERT_EQ(opalx.centers.size(), cain.centers.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(opalx), 1.0, 2.0e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(cain), 1.0, 5.0e-4);
    EXPECT_NEAR(
            LinearBreitWheelerBenchmark::histogramMean(opalx),
            LinearBreitWheelerBenchmark::histogramMean(cain),
            LinearBreitWheelerBenchmark::histogramMean(cain) * 7.0e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::histogramL1Distance(opalx, cain), 0.15);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum, ElectronJointSpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed          = 13579;

    LinearBreitWheelerBenchmark::JointHistogramConfig config;
    const auto opalx = LinearBreitWheelerBenchmark::sampleJointHistogram(
            config, LinearBreitWheelerBenchmark::FinalState::Electron, 250000);
    const auto cain =
            LinearBreitWheelerBenchmark::readJointHistogramCSV(referenceElectronJointPath());

    ASSERT_EQ(opalx.energyCentersGeV.size(), cain.energyCentersGeV.size());
    ASSERT_EQ(opalx.thetaCentersRad.size(), cain.thetaCentersRad.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramArea(opalx), 1.0, 2.0e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramArea(cain), 1.0, 2.0e-3);
    EXPECT_NEAR(
            LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(opalx),
            LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(cain),
            LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(cain) * 5.0e-2);
    EXPECT_NEAR(
            LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(opalx),
            LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(cain),
            LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(cain) * 7.0e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::jointHistogramL1Distance(opalx, cain), 0.16);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum, PositronJointSpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed          = 13579;

    LinearBreitWheelerBenchmark::JointHistogramConfig config;
    const auto opalx = LinearBreitWheelerBenchmark::sampleJointHistogram(
            config, LinearBreitWheelerBenchmark::FinalState::Positron, 250000);
    const auto cain =
            LinearBreitWheelerBenchmark::readJointHistogramCSV(referencePositronJointPath());

    ASSERT_EQ(opalx.energyCentersGeV.size(), cain.energyCentersGeV.size());
    ASSERT_EQ(opalx.thetaCentersRad.size(), cain.thetaCentersRad.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramArea(opalx), 1.0, 2.0e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramArea(cain), 1.0, 2.0e-3);
    EXPECT_NEAR(
            LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(opalx),
            LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(cain),
            LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(cain) * 5.0e-2);
    EXPECT_NEAR(
            LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(opalx),
            LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(cain),
            LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(cain) * 7.0e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::jointHistogramL1Distance(opalx, cain), 0.16);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum, FinitePhotonBeamElectronEnergySpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed          = 13579;

    LinearBreitWheelerBenchmark::FinitePhotonBeamConfig config;
    config.minValue       = 0.0;
    config.maxValue       = 0.5;
    config.sigmaThetaXRad = 1.0e-3;
    config.sigmaThetaYRad = 1.0e-3;
    const auto opalx      = LinearBreitWheelerBenchmark::sampleFinitePhotonBeamHistogram(
            config, LinearBreitWheelerBenchmark::FinalState::Electron,
            LinearBreitWheelerBenchmark::Observable::Energy, 250000);
    const auto cain = LinearBreitWheelerBenchmark::readHistogramCSV(
            referenceFinitePhotonBeamElectronEnergyPath());

    ASSERT_EQ(opalx.centers.size(), cain.centers.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(opalx), 1.0, 2.0e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(cain), 1.0, 5.0e-4);
    EXPECT_NEAR(
            LinearBreitWheelerBenchmark::histogramMean(opalx),
            LinearBreitWheelerBenchmark::histogramMean(cain),
            LinearBreitWheelerBenchmark::histogramMean(cain) * 8.0e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::histogramL1Distance(opalx, cain), 0.18);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum, FinitePhotonBeamPositronEnergySpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed          = 13579;

    LinearBreitWheelerBenchmark::FinitePhotonBeamConfig config;
    config.minValue       = 0.0;
    config.maxValue       = 0.5;
    config.sigmaThetaXRad = 1.0e-3;
    config.sigmaThetaYRad = 1.0e-3;
    const auto opalx      = LinearBreitWheelerBenchmark::sampleFinitePhotonBeamHistogram(
            config, LinearBreitWheelerBenchmark::FinalState::Positron,
            LinearBreitWheelerBenchmark::Observable::Energy, 250000);
    const auto cain = LinearBreitWheelerBenchmark::readHistogramCSV(
            referenceFinitePhotonBeamPositronEnergyPath());

    ASSERT_EQ(opalx.centers.size(), cain.centers.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(opalx), 1.0, 2.0e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(cain), 1.0, 5.0e-4);
    EXPECT_NEAR(
            LinearBreitWheelerBenchmark::histogramMean(opalx),
            LinearBreitWheelerBenchmark::histogramMean(cain),
            LinearBreitWheelerBenchmark::histogramMean(cain) * 8.0e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::histogramL1Distance(opalx, cain), 0.18);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum, FinitePhotonBeamElectronAngularSpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed          = 13579;

    LinearBreitWheelerBenchmark::FinitePhotonBeamConfig config;
    config.minValue       = 0.0;
    config.maxValue       = 0.0060;
    config.sigmaThetaXRad = 1.0e-3;
    config.sigmaThetaYRad = 1.0e-3;
    const auto opalx      = LinearBreitWheelerBenchmark::sampleFinitePhotonBeamHistogram(
            config, LinearBreitWheelerBenchmark::FinalState::Electron,
            LinearBreitWheelerBenchmark::Observable::Theta, 250000);
    const auto cain = LinearBreitWheelerBenchmark::readHistogramCSV(
            referenceFinitePhotonBeamElectronThetaPath());

    ASSERT_EQ(opalx.centers.size(), cain.centers.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(opalx), 1.0, 2.5e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(cain), 1.0, 4.0e-3);
    EXPECT_NEAR(
            LinearBreitWheelerBenchmark::histogramMean(opalx),
            LinearBreitWheelerBenchmark::histogramMean(cain),
            LinearBreitWheelerBenchmark::histogramMean(cain) * 8.0e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::histogramL1Distance(opalx, cain), 0.20);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum, FinitePhotonBeamPositronAngularSpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed          = 13579;

    LinearBreitWheelerBenchmark::FinitePhotonBeamConfig config;
    config.minValue       = 0.0;
    config.maxValue       = 0.0060;
    config.sigmaThetaXRad = 1.0e-3;
    config.sigmaThetaYRad = 1.0e-3;
    const auto opalx      = LinearBreitWheelerBenchmark::sampleFinitePhotonBeamHistogram(
            config, LinearBreitWheelerBenchmark::FinalState::Positron,
            LinearBreitWheelerBenchmark::Observable::Theta, 250000);
    const auto cain = LinearBreitWheelerBenchmark::readHistogramCSV(
            referenceFinitePhotonBeamPositronThetaPath());

    ASSERT_EQ(opalx.centers.size(), cain.centers.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(opalx), 1.0, 2.5e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(cain), 1.0, 4.0e-3);
    EXPECT_NEAR(
            LinearBreitWheelerBenchmark::histogramMean(opalx),
            LinearBreitWheelerBenchmark::histogramMean(cain),
            LinearBreitWheelerBenchmark::histogramMean(cain) * 8.0e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::histogramL1Distance(opalx, cain), 0.20);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum, FinitePhotonBeamElectronJointSpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed          = 13579;

    LinearBreitWheelerBenchmark::FinitePhotonBeamJointConfig config;
    config.energyMinGeV   = 0.0;
    config.energyMaxGeV   = 0.5;
    config.thetaMinRad    = 0.0;
    config.thetaMaxRad    = 0.0060;
    config.sigmaThetaXRad = 1.0e-3;
    config.sigmaThetaYRad = 1.0e-3;
    const auto opalx      = LinearBreitWheelerBenchmark::sampleFinitePhotonBeamJointHistogram(
            config, LinearBreitWheelerBenchmark::FinalState::Electron, 250000);
    const auto cain = LinearBreitWheelerBenchmark::readJointHistogramCSV(
            referenceFinitePhotonBeamElectronJointPath());

    ASSERT_EQ(opalx.energyCentersGeV.size(), cain.energyCentersGeV.size());
    ASSERT_EQ(opalx.thetaCentersRad.size(), cain.thetaCentersRad.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramArea(opalx), 1.0, 2.5e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramArea(cain), 1.0, 4.0e-3);
    EXPECT_NEAR(
            LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(opalx),
            LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(cain),
            LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(cain) * 8.0e-2);
    EXPECT_NEAR(
            LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(opalx),
            LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(cain),
            LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(cain) * 8.5e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::jointHistogramL1Distance(opalx, cain), 0.46);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum, FinitePhotonBeamPositronJointSpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed          = 13579;

    LinearBreitWheelerBenchmark::FinitePhotonBeamJointConfig config;
    config.energyMinGeV   = 0.0;
    config.energyMaxGeV   = 0.5;
    config.thetaMinRad    = 0.0;
    config.thetaMaxRad    = 0.0060;
    config.sigmaThetaXRad = 1.0e-3;
    config.sigmaThetaYRad = 1.0e-3;
    const auto opalx      = LinearBreitWheelerBenchmark::sampleFinitePhotonBeamJointHistogram(
            config, LinearBreitWheelerBenchmark::FinalState::Positron, 250000);
    const auto cain = LinearBreitWheelerBenchmark::readJointHistogramCSV(
            referenceFinitePhotonBeamPositronJointPath());

    ASSERT_EQ(opalx.energyCentersGeV.size(), cain.energyCentersGeV.size());
    ASSERT_EQ(opalx.thetaCentersRad.size(), cain.thetaCentersRad.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramArea(opalx), 1.0, 2.5e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramArea(cain), 1.0, 4.0e-3);
    EXPECT_NEAR(
            LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(opalx),
            LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(cain),
            LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(cain) * 8.0e-2);
    EXPECT_NEAR(
            LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(opalx),
            LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(cain),
            LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(cain) * 8.5e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::jointHistogramL1Distance(opalx, cain), 0.46);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum,
     FinitePhotonBeamEnergySpreadElectronEnergySpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed          = 13579;

    LinearBreitWheelerBenchmark::FinitePhotonBeamConfig config;
    config.minValue             = 0.0;
    config.maxValue             = 0.5;
    config.sigmaThetaXRad       = 1.0e-3;
    config.sigmaThetaYRad       = 1.0e-3;
    config.relativeEnergySpread = 1.0e-3;
    const auto opalx            = LinearBreitWheelerBenchmark::sampleFinitePhotonBeamHistogram(
            config, LinearBreitWheelerBenchmark::FinalState::Electron,
            LinearBreitWheelerBenchmark::Observable::Energy, 250000);
    const auto cain = LinearBreitWheelerBenchmark::readHistogramCSV(
            referenceFinitePhotonBeamEnergySpreadElectronEnergyPath());

    ASSERT_EQ(opalx.centers.size(), cain.centers.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(opalx), 1.0, 2.0e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(cain), 1.0, 5.0e-4);
    EXPECT_NEAR(
            LinearBreitWheelerBenchmark::histogramMean(opalx),
            LinearBreitWheelerBenchmark::histogramMean(cain),
            LinearBreitWheelerBenchmark::histogramMean(cain) * 8.0e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::histogramL1Distance(opalx, cain), 0.18);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum,
     FinitePhotonBeamEnergySpreadPositronEnergySpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed          = 13579;

    LinearBreitWheelerBenchmark::FinitePhotonBeamConfig config;
    config.minValue             = 0.0;
    config.maxValue             = 0.5;
    config.sigmaThetaXRad       = 1.0e-3;
    config.sigmaThetaYRad       = 1.0e-3;
    config.relativeEnergySpread = 1.0e-3;
    const auto opalx            = LinearBreitWheelerBenchmark::sampleFinitePhotonBeamHistogram(
            config, LinearBreitWheelerBenchmark::FinalState::Positron,
            LinearBreitWheelerBenchmark::Observable::Energy, 250000);
    const auto cain = LinearBreitWheelerBenchmark::readHistogramCSV(
            referenceFinitePhotonBeamEnergySpreadPositronEnergyPath());

    ASSERT_EQ(opalx.centers.size(), cain.centers.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(opalx), 1.0, 2.0e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(cain), 1.0, 5.0e-4);
    EXPECT_NEAR(
            LinearBreitWheelerBenchmark::histogramMean(opalx),
            LinearBreitWheelerBenchmark::histogramMean(cain),
            LinearBreitWheelerBenchmark::histogramMean(cain) * 8.0e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::histogramL1Distance(opalx, cain), 0.18);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum,
     FinitePhotonBeamEnergySpreadElectronAngularSpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed          = 13579;

    LinearBreitWheelerBenchmark::FinitePhotonBeamConfig config;
    config.minValue             = 0.0;
    config.maxValue             = 0.0060;
    config.sigmaThetaXRad       = 1.0e-3;
    config.sigmaThetaYRad       = 1.0e-3;
    config.relativeEnergySpread = 1.0e-3;
    const auto opalx            = LinearBreitWheelerBenchmark::sampleFinitePhotonBeamHistogram(
            config, LinearBreitWheelerBenchmark::FinalState::Electron,
            LinearBreitWheelerBenchmark::Observable::Theta, 250000);
    const auto cain = LinearBreitWheelerBenchmark::readHistogramCSV(
            referenceFinitePhotonBeamEnergySpreadElectronThetaPath());

    ASSERT_EQ(opalx.centers.size(), cain.centers.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(opalx), 1.0, 2.5e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(cain), 1.0, 4.0e-3);
    EXPECT_NEAR(
            LinearBreitWheelerBenchmark::histogramMean(opalx),
            LinearBreitWheelerBenchmark::histogramMean(cain),
            LinearBreitWheelerBenchmark::histogramMean(cain) * 8.0e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::histogramL1Distance(opalx, cain), 0.20);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum,
     FinitePhotonBeamEnergySpreadPositronAngularSpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed          = 13579;

    LinearBreitWheelerBenchmark::FinitePhotonBeamConfig config;
    config.minValue             = 0.0;
    config.maxValue             = 0.0060;
    config.sigmaThetaXRad       = 1.0e-3;
    config.sigmaThetaYRad       = 1.0e-3;
    config.relativeEnergySpread = 1.0e-3;
    const auto opalx            = LinearBreitWheelerBenchmark::sampleFinitePhotonBeamHistogram(
            config, LinearBreitWheelerBenchmark::FinalState::Positron,
            LinearBreitWheelerBenchmark::Observable::Theta, 250000);
    const auto cain = LinearBreitWheelerBenchmark::readHistogramCSV(
            referenceFinitePhotonBeamEnergySpreadPositronThetaPath());

    ASSERT_EQ(opalx.centers.size(), cain.centers.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(opalx), 1.0, 2.5e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(cain), 1.0, 4.0e-3);
    EXPECT_NEAR(
            LinearBreitWheelerBenchmark::histogramMean(opalx),
            LinearBreitWheelerBenchmark::histogramMean(cain),
            LinearBreitWheelerBenchmark::histogramMean(cain) * 8.0e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::histogramL1Distance(opalx, cain), 0.20);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum,
     FinitePhotonBeamEnergySpreadElectronJointSpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed          = 13579;

    LinearBreitWheelerBenchmark::FinitePhotonBeamJointConfig config;
    config.energyMinGeV         = 0.0;
    config.energyMaxGeV         = 0.5;
    config.thetaMinRad          = 0.0;
    config.thetaMaxRad          = 0.0060;
    config.sigmaThetaXRad       = 1.0e-3;
    config.sigmaThetaYRad       = 1.0e-3;
    config.relativeEnergySpread = 1.0e-3;
    const auto opalx            = LinearBreitWheelerBenchmark::sampleFinitePhotonBeamJointHistogram(
            config, LinearBreitWheelerBenchmark::FinalState::Electron, 250000);
    const auto cain = LinearBreitWheelerBenchmark::readJointHistogramCSV(
            referenceFinitePhotonBeamEnergySpreadElectronJointPath());

    ASSERT_EQ(opalx.energyCentersGeV.size(), cain.energyCentersGeV.size());
    ASSERT_EQ(opalx.thetaCentersRad.size(), cain.thetaCentersRad.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramArea(opalx), 1.0, 2.5e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramArea(cain), 1.0, 4.0e-3);
    EXPECT_NEAR(
            LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(opalx),
            LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(cain),
            LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(cain) * 8.5e-2);
    EXPECT_NEAR(
            LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(opalx),
            LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(cain),
            LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(cain) * 8.5e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::jointHistogramL1Distance(opalx, cain), 0.46);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum,
     FinitePhotonBeamEnergySpreadPositronJointSpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed          = 13579;

    LinearBreitWheelerBenchmark::FinitePhotonBeamJointConfig config;
    config.energyMinGeV         = 0.0;
    config.energyMaxGeV         = 0.5;
    config.thetaMinRad          = 0.0;
    config.thetaMaxRad          = 0.0060;
    config.sigmaThetaXRad       = 1.0e-3;
    config.sigmaThetaYRad       = 1.0e-3;
    config.relativeEnergySpread = 1.0e-3;
    const auto opalx            = LinearBreitWheelerBenchmark::sampleFinitePhotonBeamJointHistogram(
            config, LinearBreitWheelerBenchmark::FinalState::Positron, 250000);
    const auto cain = LinearBreitWheelerBenchmark::readJointHistogramCSV(
            referenceFinitePhotonBeamEnergySpreadPositronJointPath());

    ASSERT_EQ(opalx.energyCentersGeV.size(), cain.energyCentersGeV.size());
    ASSERT_EQ(opalx.thetaCentersRad.size(), cain.thetaCentersRad.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramArea(opalx), 1.0, 2.5e-2);
    EXPECT_NEAR(LinearBreitWheelerBenchmark::jointHistogramArea(cain), 1.0, 4.0e-3);
    EXPECT_NEAR(
            LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(opalx),
            LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(cain),
            LinearBreitWheelerBenchmark::jointHistogramMeanEnergyGeV(cain) * 8.5e-2);
    EXPECT_NEAR(
            LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(opalx),
            LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(cain),
            LinearBreitWheelerBenchmark::jointHistogramMeanThetaRad(cain) * 8.5e-2);
    EXPECT_LT(LinearBreitWheelerBenchmark::jointHistogramL1Distance(opalx, cain), 0.46);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum, OverlapWeightedElectronEnergySpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed          = 13579;

    LinearBreitWheelerBenchmark::FinitePhotonBeamConfig config;
    config.minValue         = 0.0;
    config.maxValue         = 0.5;
    config.sigmaThetaXRad   = 1.0e-3;
    config.sigmaThetaYRad   = 1.0e-3;
    config.sigmaX_m         = 1.0e-6;
    config.sigmaY_m         = 1.0e-6;
    config.sigmaS_m         = 5.0e-12 * Physics::c;
    config.laserRayleighX_m = 1.0e-6;
    config.laserRayleighY_m = 1.0e-6;
    config.laserSigmaT_m    = 5.0e-12 * Physics::c;
    config.overlapWeighting = true;
    const auto opalx        = LinearBreitWheelerBenchmark::sampleFinitePhotonBeamHistogram(
            config, LinearBreitWheelerBenchmark::FinalState::Electron,
            LinearBreitWheelerBenchmark::Observable::Energy, 250000);
    const auto cain =
            LinearBreitWheelerBenchmark::readHistogramCSV(referenceOverlapElectronEnergyPath());

    ASSERT_EQ(opalx.centers.size(), cain.centers.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(opalx), 1.0, 2.0e-2);
    EXPECT_NEAR(
            LinearBreitWheelerBenchmark::histogramMean(opalx),
            LinearBreitWheelerBenchmark::histogramMean(cain),
            LinearBreitWheelerBenchmark::histogramMean(cain) * 1.2e-1);
    EXPECT_LT(LinearBreitWheelerBenchmark::histogramL1Distance(opalx, cain), 0.28);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum, OverlapWeightedPositronEnergySpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed          = 13579;

    LinearBreitWheelerBenchmark::FinitePhotonBeamConfig config;
    config.minValue         = 0.0;
    config.maxValue         = 0.5;
    config.sigmaThetaXRad   = 1.0e-3;
    config.sigmaThetaYRad   = 1.0e-3;
    config.sigmaX_m         = 1.0e-6;
    config.sigmaY_m         = 1.0e-6;
    config.sigmaS_m         = 5.0e-12 * Physics::c;
    config.laserRayleighX_m = 1.0e-6;
    config.laserRayleighY_m = 1.0e-6;
    config.laserSigmaT_m    = 5.0e-12 * Physics::c;
    config.overlapWeighting = true;
    const auto opalx        = LinearBreitWheelerBenchmark::sampleFinitePhotonBeamHistogram(
            config, LinearBreitWheelerBenchmark::FinalState::Positron,
            LinearBreitWheelerBenchmark::Observable::Energy, 250000);
    const auto cain =
            LinearBreitWheelerBenchmark::readHistogramCSV(referenceOverlapPositronEnergyPath());

    ASSERT_EQ(opalx.centers.size(), cain.centers.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(opalx), 1.0, 2.0e-2);
    EXPECT_NEAR(
            LinearBreitWheelerBenchmark::histogramMean(opalx),
            LinearBreitWheelerBenchmark::histogramMean(cain),
            LinearBreitWheelerBenchmark::histogramMean(cain) * 1.2e-1);
    EXPECT_LT(LinearBreitWheelerBenchmark::histogramL1Distance(opalx, cain), 0.28);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum, OverlapWeightedElectronAngularSpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed          = 13579;

    LinearBreitWheelerBenchmark::FinitePhotonBeamConfig config;
    config.minValue         = 0.0;
    config.maxValue         = 0.0060;
    config.sigmaThetaXRad   = 1.0e-3;
    config.sigmaThetaYRad   = 1.0e-3;
    config.sigmaX_m         = 1.0e-6;
    config.sigmaY_m         = 1.0e-6;
    config.sigmaS_m         = 5.0e-12 * Physics::c;
    config.laserRayleighX_m = 1.0e-6;
    config.laserRayleighY_m = 1.0e-6;
    config.laserSigmaT_m    = 5.0e-12 * Physics::c;
    config.overlapWeighting = true;
    const auto opalx        = LinearBreitWheelerBenchmark::sampleFinitePhotonBeamHistogram(
            config, LinearBreitWheelerBenchmark::FinalState::Electron,
            LinearBreitWheelerBenchmark::Observable::Theta, 250000);
    const auto cain =
            LinearBreitWheelerBenchmark::readHistogramCSV(referenceOverlapElectronThetaPath());

    ASSERT_EQ(opalx.centers.size(), cain.centers.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(opalx), 1.0, 2.5e-2);
    EXPECT_NEAR(
            LinearBreitWheelerBenchmark::histogramMean(opalx),
            LinearBreitWheelerBenchmark::histogramMean(cain),
            LinearBreitWheelerBenchmark::histogramMean(cain) * 1.2e-1);
    EXPECT_LT(LinearBreitWheelerBenchmark::histogramL1Distance(opalx, cain), 0.32);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheelerSpectrum, OverlapWeightedPositronAngularSpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed          = 13579;

    LinearBreitWheelerBenchmark::FinitePhotonBeamConfig config;
    config.minValue         = 0.0;
    config.maxValue         = 0.0060;
    config.sigmaThetaXRad   = 1.0e-3;
    config.sigmaThetaYRad   = 1.0e-3;
    config.sigmaX_m         = 1.0e-6;
    config.sigmaY_m         = 1.0e-6;
    config.sigmaS_m         = 5.0e-12 * Physics::c;
    config.laserRayleighX_m = 1.0e-6;
    config.laserRayleighY_m = 1.0e-6;
    config.laserSigmaT_m    = 5.0e-12 * Physics::c;
    config.overlapWeighting = true;
    const auto opalx        = LinearBreitWheelerBenchmark::sampleFinitePhotonBeamHistogram(
            config, LinearBreitWheelerBenchmark::FinalState::Positron,
            LinearBreitWheelerBenchmark::Observable::Theta, 250000);
    const auto cain =
            LinearBreitWheelerBenchmark::readHistogramCSV(referenceOverlapPositronThetaPath());

    ASSERT_EQ(opalx.centers.size(), cain.centers.size());
    EXPECT_NEAR(LinearBreitWheelerBenchmark::histogramArea(opalx), 1.0, 2.5e-2);
    EXPECT_NEAR(
            LinearBreitWheelerBenchmark::histogramMean(opalx),
            LinearBreitWheelerBenchmark::histogramMean(cain),
            LinearBreitWheelerBenchmark::histogramMean(cain) * 1.2e-1);
    EXPECT_LT(LinearBreitWheelerBenchmark::histogramL1Distance(opalx, cain), 0.32);

    Options::seed = previousSeed;
}
