/**
 * @file TestLinearComptonSpectrum.cpp
 * @brief Regression tests for CAIN-validated linear-Compton energy and angle spectra.
 *
 * This file mixes two validation styles:
 *
 * - direct deterministic comparisons against stored CAIN histograms for the
 *   single-electron benchmark, and
 * - subprocess-driven regeneration of finite-beam benchmark CSV files through
 *   `LinearComptonSpectrumBenchmark`, followed by comparison against stored
 *   CAIN references.
 *
 * The histogram checks use area, mean, and L1-distance tolerances instead of
 * pointwise equality. That keeps the tests stable across small Monte Carlo and
 * quadrature differences while still detecting shape regressions in the photon
 * spectra.
 */

#include "LinearComptonSpectrumCommon.h"
#include "Utilities/Options.h"
#include "gtest/gtest.h"

#include <cstdlib>
#include <filesystem>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {
/** @brief Stored CAIN energy-spectrum reference for the weak-field single-electron case. */
std::filesystem::path referenceSpectrumPath() {
    return std::filesystem::path(OPALX_TEST_SOURCE_DIR)
        / "data"
        / "cain_linear_compton_90deg_xi029_histogram.csv";
}

/** @brief Stored CAIN joint E_gamma-theta_gamma reference for the weak-field single-electron case. */
std::filesystem::path referenceJointSpectrumPath() {
    return std::filesystem::path(OPALX_TEST_SOURCE_DIR)
        / "data"
        / "cain_linear_compton_90deg_xi029_joint_histogram.csv";
}

/** @brief Stored CAIN lab-angle reference for the weak-field single-electron case. */
std::filesystem::path referenceAngularSpectrumPath() {
    return std::filesystem::path(OPALX_TEST_SOURCE_DIR)
        / "data"
        / "cain_linear_compton_90deg_xi029_theta_histogram.csv";
}

/** @brief Stored CAIN energy-spectrum reference for the finite-emittance benchmark. */
std::filesystem::path referenceFiniteBeamSpectrumPath() {
    return std::filesystem::path(OPALX_TEST_SOURCE_DIR)
        / "data"
        / "cain_linear_compton_90deg_xi029_finite_beam_histogram.csv";
}

/** @brief Stored CAIN lab-angle reference for the finite-emittance benchmark. */
std::filesystem::path referenceFiniteBeamAngularSpectrumPath() {
    return std::filesystem::path(OPALX_TEST_SOURCE_DIR)
        / "data"
        / "cain_linear_compton_90deg_xi029_finite_beam_theta_histogram.csv";
}

/** @brief Stored CAIN joint-spectrum reference for the finite-emittance benchmark. */
std::filesystem::path referenceFiniteBeamJointSpectrumPath() {
    return std::filesystem::path(OPALX_TEST_SOURCE_DIR)
        / "data"
        / "cain_linear_compton_90deg_xi029_finite_beam_joint_histogram.csv";
}

/** @brief Stored CAIN energy-spectrum reference for the finite-beam plus energy-spread benchmark. */
std::filesystem::path referenceFiniteBeamEnergySpreadSpectrumPath() {
    return std::filesystem::path(OPALX_TEST_SOURCE_DIR)
        / "data"
        / "cain_linear_compton_90deg_xi029_finite_beam_energy_spread_histogram.csv";
}

/** @brief Stored CAIN lab-angle reference for the finite-beam plus energy-spread benchmark. */
std::filesystem::path referenceFiniteBeamEnergySpreadAngularSpectrumPath() {
    return std::filesystem::path(OPALX_TEST_SOURCE_DIR)
        / "data"
        / "cain_linear_compton_90deg_xi029_finite_beam_energy_spread_theta_histogram.csv";
}

/** @brief Stored CAIN joint-spectrum reference for the finite-beam plus energy-spread benchmark. */
std::filesystem::path referenceFiniteBeamEnergySpreadJointSpectrumPath() {
    return std::filesystem::path(OPALX_TEST_SOURCE_DIR)
        / "data"
        / "cain_linear_compton_90deg_xi029_finite_beam_energy_spread_joint_histogram.csv";
}

/** @brief Stored CAIN energy-spectrum reference for the overlap-restricted finite-beam benchmark. */
std::filesystem::path referenceFiniteBeamOverlapSpectrumPath() {
    return std::filesystem::path(OPALX_TEST_SOURCE_DIR)
        / "data"
        / "cain_linear_compton_90deg_xi029_finite_beam_overlap_histogram.csv";
}

/** @brief Stored CAIN angular reference for the overlap-restricted finite-beam benchmark. */
std::filesystem::path referenceFiniteBeamOverlapAngularSpectrumPath() {
    return std::filesystem::path(OPALX_TEST_SOURCE_DIR)
        / "data"
        / "cain_linear_compton_90deg_xi029_finite_beam_overlap_theta_histogram.csv";
}

/** @brief Path to the standalone benchmark executable used by the finite-beam regression cases. */
std::filesystem::path benchmarkExecutablePath() {
    return std::filesystem::path(OPALX_LINEAR_COMPTON_BENCHMARK);
}

/**
 * @brief Regenerate one finite-beam OPALX histogram through the standalone benchmark executable.
 *
 * The finite-beam tests intentionally invoke the external benchmark binary so
 * they cover the same code path used to regenerate the checked-in CSV
 * references. The helper keeps the command line fixed except for the output
 * filename, the observable selector, and optional extra benchmark parameters
 * such as energy spread.
 */
std::filesystem::path runFiniteBeamBenchmark(const std::string& filename,
                                           bool angular,
                                           bool joint = false,
                                           const std::string& extraOptions = "") {
    const auto output = std::filesystem::temp_directory_path() / filename;
    std::ostringstream command;
    command << benchmarkExecutablePath().string() << ' '
            << output.string()
            << " --finite-beam --beam-particles 100000 --seed 13579";
    if (!extraOptions.empty()) {
        command << ' ' << extraOptions;
    }
    if (joint) {
        command << " --joint";
    } else if (angular) {
        command << " --angular";
    }

    const int returnCode = std::system(command.str().c_str());
    if (returnCode != 0) {
        throw std::runtime_error("Failed to run finite-beam benchmark command: " + command.str());
    }
    return output;
}
}  // namespace

TEST(TestLinearComptonSpectrum, WeakFieldSpectrumMatchesCainReference) {
    LinearComptonBenchmark::SpectrumConfig config;
    const auto opalxSpectrum = LinearComptonBenchmark::integrateLabSpectrum(config);
    const auto cainSpectrum = LinearComptonBenchmark::readSpectrumCSV(referenceSpectrumPath());

    ASSERT_EQ(opalxSpectrum.centersGeV.size(), cainSpectrum.centersGeV.size());
    EXPECT_NEAR(LinearComptonBenchmark::histogramArea(opalxSpectrum), 1.0, 5.0e-4);
    EXPECT_NEAR(LinearComptonBenchmark::histogramArea(cainSpectrum), 1.0, 5.0e-4);

    const double opalxMean = LinearComptonBenchmark::histogramMeanEnergyGeV(opalxSpectrum);
    const double cainMean = LinearComptonBenchmark::histogramMeanEnergyGeV(cainSpectrum);
    EXPECT_NEAR(opalxMean, cainMean, cainMean * 3.0e-2);

    const double l1Distance = LinearComptonBenchmark::histogramL1Distance(opalxSpectrum,
                                                                          cainSpectrum);
    EXPECT_LT(l1Distance, 0.1);
}

TEST(TestLinearComptonSpectrum, WeakFieldSampledSpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed = 13579;

    LinearComptonBenchmark::SpectrumConfig config;
    const auto opalxSpectrum = LinearComptonBenchmark::sampleLabSpectrum(config, 250000);
    const auto cainSpectrum = LinearComptonBenchmark::readSpectrumCSV(referenceSpectrumPath());

    ASSERT_EQ(opalxSpectrum.centersGeV.size(), cainSpectrum.centersGeV.size());
    EXPECT_NEAR(LinearComptonBenchmark::histogramArea(opalxSpectrum), 1.0, 1.5e-2);

    const double opalxMean = LinearComptonBenchmark::histogramMeanEnergyGeV(opalxSpectrum);
    const double cainMean = LinearComptonBenchmark::histogramMeanEnergyGeV(cainSpectrum);
    EXPECT_NEAR(opalxMean, cainMean, cainMean * 3.5e-2);

    const double l1Distance = LinearComptonBenchmark::histogramL1Distance(opalxSpectrum,
                                                                          cainSpectrum);
    EXPECT_LT(l1Distance, 0.14);

    Options::seed = previousSeed;
}

TEST(TestLinearComptonSpectrum, WeakFieldAngularSpectrumMatchesCainReference) {
    LinearComptonBenchmark::AngleConfig config;
    const auto opalxSpectrum = LinearComptonBenchmark::integrateLabAngularSpectrum(config);
    const auto cainSpectrum = LinearComptonBenchmark::readAngleCSV(referenceAngularSpectrumPath());

    ASSERT_EQ(opalxSpectrum.centersRad.size(), cainSpectrum.centersRad.size());
    EXPECT_NEAR(LinearComptonBenchmark::angleHistogramArea(opalxSpectrum), 1.0, 5.0e-4);
    EXPECT_NEAR(LinearComptonBenchmark::angleHistogramArea(cainSpectrum), 1.0, 5.0e-4);

    const double opalxMean = LinearComptonBenchmark::angleHistogramMeanRad(opalxSpectrum);
    const double cainMean = LinearComptonBenchmark::angleHistogramMeanRad(cainSpectrum);
    EXPECT_NEAR(opalxMean, cainMean, cainMean * 5.0e-2);

    const double l1Distance = LinearComptonBenchmark::angleHistogramL1Distance(opalxSpectrum,
                                                                               cainSpectrum);
    EXPECT_LT(l1Distance, 0.12);
}

TEST(TestLinearComptonSpectrum, WeakFieldSampledAngularSpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed = 13579;

    LinearComptonBenchmark::AngleConfig config;
    const auto opalxSpectrum = LinearComptonBenchmark::sampleLabAngularSpectrum(config, 250000);
    const auto cainSpectrum = LinearComptonBenchmark::readAngleCSV(referenceAngularSpectrumPath());

    ASSERT_EQ(opalxSpectrum.centersRad.size(), cainSpectrum.centersRad.size());
    EXPECT_NEAR(LinearComptonBenchmark::angleHistogramArea(opalxSpectrum), 1.0, 2.0e-2);

    const double opalxMean = LinearComptonBenchmark::angleHistogramMeanRad(opalxSpectrum);
    const double cainMean = LinearComptonBenchmark::angleHistogramMeanRad(cainSpectrum);
    EXPECT_NEAR(opalxMean, cainMean, cainMean * 5.0e-2);

    const double l1Distance = LinearComptonBenchmark::angleHistogramL1Distance(opalxSpectrum,
                                                                               cainSpectrum);
    EXPECT_LT(l1Distance, 0.10);

    Options::seed = previousSeed;
}

TEST(TestLinearComptonSpectrum, WeakFieldJointSpectrumMatchesCainReference) {
    LinearComptonBenchmark::JointConfig config;
    const auto opalxSpectrum = LinearComptonBenchmark::integrateLabJointSpectrum(config);
    const auto cainSpectrum = LinearComptonBenchmark::readJointCSV(referenceJointSpectrumPath());

    ASSERT_EQ(opalxSpectrum.energyCentersGeV.size(), cainSpectrum.energyCentersGeV.size());
    ASSERT_EQ(opalxSpectrum.thetaCentersRad.size(), cainSpectrum.thetaCentersRad.size());
    EXPECT_NEAR(LinearComptonBenchmark::jointHistogramArea(opalxSpectrum), 1.0, 2.0e-3);
    EXPECT_NEAR(LinearComptonBenchmark::jointHistogramArea(cainSpectrum), 1.0, 2.0e-3);

    const double opalxMeanEnergy = LinearComptonBenchmark::jointHistogramMeanEnergyGeV(opalxSpectrum);
    const double cainMeanEnergy = LinearComptonBenchmark::jointHistogramMeanEnergyGeV(cainSpectrum);
    EXPECT_NEAR(opalxMeanEnergy, cainMeanEnergy, cainMeanEnergy * 3.0e-2);

    const double opalxMeanTheta = LinearComptonBenchmark::jointHistogramMeanThetaRad(opalxSpectrum);
    const double cainMeanTheta = LinearComptonBenchmark::jointHistogramMeanThetaRad(cainSpectrum);
    EXPECT_NEAR(opalxMeanTheta, cainMeanTheta, cainMeanTheta * 5.0e-2);

    const double l1Distance = LinearComptonBenchmark::jointHistogramL1Distance(opalxSpectrum,
                                                                               cainSpectrum);
    EXPECT_LT(l1Distance, 0.16);
}

TEST(TestLinearComptonSpectrum, WeakFieldSampledJointSpectrumMatchesCainReference) {
    const int previousSeed = Options::seed;
    Options::seed = 13579;

    LinearComptonBenchmark::JointConfig config;
    const auto opalxSpectrum = LinearComptonBenchmark::sampleLabJointSpectrum(config, 250000);
    const auto cainSpectrum = LinearComptonBenchmark::readJointCSV(referenceJointSpectrumPath());

    ASSERT_EQ(opalxSpectrum.energyCentersGeV.size(), cainSpectrum.energyCentersGeV.size());
    ASSERT_EQ(opalxSpectrum.thetaCentersRad.size(), cainSpectrum.thetaCentersRad.size());
    EXPECT_NEAR(LinearComptonBenchmark::jointHistogramArea(opalxSpectrum), 1.0, 1.5e-2);

    const double opalxMeanEnergy = LinearComptonBenchmark::jointHistogramMeanEnergyGeV(opalxSpectrum);
    const double cainMeanEnergy = LinearComptonBenchmark::jointHistogramMeanEnergyGeV(cainSpectrum);
    EXPECT_NEAR(opalxMeanEnergy, cainMeanEnergy, cainMeanEnergy * 3.5e-2);

    const double opalxMeanTheta = LinearComptonBenchmark::jointHistogramMeanThetaRad(opalxSpectrum);
    const double cainMeanTheta = LinearComptonBenchmark::jointHistogramMeanThetaRad(cainSpectrum);
    EXPECT_NEAR(opalxMeanTheta, cainMeanTheta, cainMeanTheta * 6.0e-2);

    const double l1Distance = LinearComptonBenchmark::jointHistogramL1Distance(opalxSpectrum,
                                                                               cainSpectrum);
    EXPECT_LT(l1Distance, 0.09);

    Options::seed = previousSeed;
}

/**
 * The finite-beam regressions allow a small normalization and shape tolerance
 * because they combine Gaussian beam sampling with stochastic Compton event
 * generation. The mean and L1 checks are the primary regression guards.
 */
TEST(TestLinearComptonSpectrum, FiniteBeamSpectrumMatchesCainReference) {
    const auto output = runFiniteBeamBenchmark("opalx-linear-compton-90deg-xi029-finite-beam-histogram.csv",
                                               false);
    const auto opalxSpectrum = LinearComptonBenchmark::readSpectrumCSV(output);
    const auto cainSpectrum = LinearComptonBenchmark::readSpectrumCSV(referenceFiniteBeamSpectrumPath());

    ASSERT_EQ(opalxSpectrum.centersGeV.size(), cainSpectrum.centersGeV.size());
    EXPECT_NEAR(LinearComptonBenchmark::histogramArea(opalxSpectrum), 1.0, 2.0e-2);
    EXPECT_NEAR(LinearComptonBenchmark::histogramArea(cainSpectrum), 1.0, 5.0e-4);

    const double opalxMean = LinearComptonBenchmark::histogramMeanEnergyGeV(opalxSpectrum);
    const double cainMean = LinearComptonBenchmark::histogramMeanEnergyGeV(cainSpectrum);
    EXPECT_NEAR(opalxMean, cainMean, cainMean * 2.0e-2);

    const double l1Distance = LinearComptonBenchmark::histogramL1Distance(opalxSpectrum,
                                                                          cainSpectrum);
    EXPECT_LT(l1Distance, 0.05);
}

TEST(TestLinearComptonSpectrum, FiniteBeamAngularSpectrumMatchesCainReference) {
    const auto output = runFiniteBeamBenchmark("opalx-linear-compton-90deg-xi029-finite-beam-theta-histogram.csv",
                                               true);
    const auto opalxSpectrum = LinearComptonBenchmark::readAngleCSV(output);
    const auto cainSpectrum = LinearComptonBenchmark::readAngleCSV(referenceFiniteBeamAngularSpectrumPath());

    ASSERT_EQ(opalxSpectrum.centersRad.size(), cainSpectrum.centersRad.size());
    EXPECT_NEAR(LinearComptonBenchmark::angleHistogramArea(opalxSpectrum), 1.0, 2.5e-2);
    EXPECT_NEAR(LinearComptonBenchmark::angleHistogramArea(cainSpectrum), 1.0, 5.0e-4);

    const double opalxMean = LinearComptonBenchmark::angleHistogramMeanRad(opalxSpectrum);
    const double cainMean = LinearComptonBenchmark::angleHistogramMeanRad(cainSpectrum);
    EXPECT_NEAR(opalxMean, cainMean, cainMean * 2.0e-2);

    const double l1Distance = LinearComptonBenchmark::angleHistogramL1Distance(opalxSpectrum,
                                                                               cainSpectrum);
    EXPECT_LT(l1Distance, 0.03);
}


TEST(TestLinearComptonSpectrum, FiniteBeamJointSpectrumMatchesCainReference) {
    const auto output = runFiniteBeamBenchmark(
        "opalx-linear-compton-90deg-xi029-finite-beam-joint-histogram.csv",
        false,
        true);
    const auto opalxSpectrum = LinearComptonBenchmark::readJointCSV(output);
    const auto cainSpectrum = LinearComptonBenchmark::readJointCSV(referenceFiniteBeamJointSpectrumPath());

    ASSERT_EQ(opalxSpectrum.energyCentersGeV.size(), cainSpectrum.energyCentersGeV.size());
    ASSERT_EQ(opalxSpectrum.thetaCentersRad.size(), cainSpectrum.thetaCentersRad.size());
    EXPECT_NEAR(LinearComptonBenchmark::jointHistogramArea(opalxSpectrum), 1.0, 2.5e-2);
    EXPECT_NEAR(LinearComptonBenchmark::jointHistogramArea(cainSpectrum), 1.0, 1.0e-3);

    EXPECT_NEAR(LinearComptonBenchmark::jointHistogramMeanEnergyGeV(opalxSpectrum),
                LinearComptonBenchmark::jointHistogramMeanEnergyGeV(cainSpectrum),
                LinearComptonBenchmark::jointHistogramMeanEnergyGeV(cainSpectrum) * 2.0e-2);
    EXPECT_NEAR(LinearComptonBenchmark::jointHistogramMeanThetaRad(opalxSpectrum),
                LinearComptonBenchmark::jointHistogramMeanThetaRad(cainSpectrum),
                LinearComptonBenchmark::jointHistogramMeanThetaRad(cainSpectrum) * 5.0e-2);
    EXPECT_LT(LinearComptonBenchmark::jointHistogramL1Distance(opalxSpectrum, cainSpectrum), 0.08);
}

/**
 * This variant adds a small relative energy spread to the same finite-emittance
 * beam. The tolerance is kept at the same scale as the zero-spread finite-beam
 * case because the stored CAIN reference has high enough statistics to keep the
 * broadened distributions stable.
 */
TEST(TestLinearComptonSpectrum, FiniteBeamEnergySpreadSpectrumMatchesCainReference) {
    const auto output = runFiniteBeamBenchmark(
        "opalx-linear-compton-90deg-xi029-finite-beam-energy-spread-histogram.csv",
        false,
        false,
        "--beam-relative-energy-spread 0.001");
    const auto opalxSpectrum = LinearComptonBenchmark::readSpectrumCSV(output);
    const auto cainSpectrum = LinearComptonBenchmark::readSpectrumCSV(
        referenceFiniteBeamEnergySpreadSpectrumPath());

    ASSERT_EQ(opalxSpectrum.centersGeV.size(), cainSpectrum.centersGeV.size());
    EXPECT_NEAR(LinearComptonBenchmark::histogramArea(opalxSpectrum), 1.0, 2.0e-2);
    EXPECT_NEAR(LinearComptonBenchmark::histogramArea(cainSpectrum), 1.0, 5.0e-4);

    const double opalxMean = LinearComptonBenchmark::histogramMeanEnergyGeV(opalxSpectrum);
    const double cainMean = LinearComptonBenchmark::histogramMeanEnergyGeV(cainSpectrum);
    EXPECT_NEAR(opalxMean, cainMean, cainMean * 2.0e-2);

    const double l1Distance = LinearComptonBenchmark::histogramL1Distance(opalxSpectrum,
                                                                          cainSpectrum);
    EXPECT_LT(l1Distance, 0.05);
}

TEST(TestLinearComptonSpectrum, FiniteBeamEnergySpreadAngularSpectrumMatchesCainReference) {
    const auto output = runFiniteBeamBenchmark(
        "opalx-linear-compton-90deg-xi029-finite-beam-energy-spread-theta-histogram.csv",
        true,
        false,
        "--beam-relative-energy-spread 0.001");
    const auto opalxSpectrum = LinearComptonBenchmark::readAngleCSV(output);
    const auto cainSpectrum = LinearComptonBenchmark::readAngleCSV(
        referenceFiniteBeamEnergySpreadAngularSpectrumPath());

    ASSERT_EQ(opalxSpectrum.centersRad.size(), cainSpectrum.centersRad.size());
    EXPECT_NEAR(LinearComptonBenchmark::angleHistogramArea(opalxSpectrum), 1.0, 2.5e-2);
    EXPECT_NEAR(LinearComptonBenchmark::angleHistogramArea(cainSpectrum), 1.0, 5.0e-4);

    const double opalxMean = LinearComptonBenchmark::angleHistogramMeanRad(opalxSpectrum);
    const double cainMean = LinearComptonBenchmark::angleHistogramMeanRad(cainSpectrum);
    EXPECT_NEAR(opalxMean, cainMean, cainMean * 2.0e-2);

    const double l1Distance = LinearComptonBenchmark::angleHistogramL1Distance(opalxSpectrum,
                                                                               cainSpectrum);
    EXPECT_LT(l1Distance, 0.03);
}

TEST(TestLinearComptonSpectrum, FiniteBeamEnergySpreadJointSpectrumMatchesCainReference) {
    const auto output = runFiniteBeamBenchmark(
        "opalx-linear-compton-90deg-xi029-finite-beam-energy-spread-joint-histogram.csv",
        false,
        true,
        "--beam-relative-energy-spread 0.001");
    const auto opalxSpectrum = LinearComptonBenchmark::readJointCSV(output);
    const auto cainSpectrum = LinearComptonBenchmark::readJointCSV(
        referenceFiniteBeamEnergySpreadJointSpectrumPath());

    ASSERT_EQ(opalxSpectrum.energyCentersGeV.size(), cainSpectrum.energyCentersGeV.size());
    ASSERT_EQ(opalxSpectrum.thetaCentersRad.size(), cainSpectrum.thetaCentersRad.size());
    EXPECT_NEAR(LinearComptonBenchmark::jointHistogramArea(opalxSpectrum), 1.0, 2.5e-2);
    EXPECT_NEAR(LinearComptonBenchmark::jointHistogramArea(cainSpectrum), 1.0, 1.5e-3);

    EXPECT_NEAR(LinearComptonBenchmark::jointHistogramMeanEnergyGeV(opalxSpectrum),
                LinearComptonBenchmark::jointHistogramMeanEnergyGeV(cainSpectrum),
                LinearComptonBenchmark::jointHistogramMeanEnergyGeV(cainSpectrum) * 2.0e-2);
    EXPECT_NEAR(LinearComptonBenchmark::jointHistogramMeanThetaRad(opalxSpectrum),
                LinearComptonBenchmark::jointHistogramMeanThetaRad(cainSpectrum),
                LinearComptonBenchmark::jointHistogramMeanThetaRad(cainSpectrum) * 5.0e-2);
    EXPECT_LT(LinearComptonBenchmark::jointHistogramL1Distance(opalxSpectrum, cainSpectrum), 0.08);
}

TEST(TestLinearComptonSpectrum, FiniteBeamOverlapSpectrumMatchesCainReference) {
    const auto output = runFiniteBeamBenchmark(
        "opalx-linear-compton-90deg-xi029-finite-beam-overlap-histogram.csv",
        false,
        false,
        "--overlap-weighting --beam-sigma-longitudinal 0.000299792 --laser-rayleigh 12.5 --laser-sigma-t 0.000299792");
    const auto opalxSpectrum = LinearComptonBenchmark::readSpectrumCSV(output);
    const auto cainSpectrum = LinearComptonBenchmark::readSpectrumCSV(referenceFiniteBeamOverlapSpectrumPath());

    ASSERT_EQ(opalxSpectrum.centersGeV.size(), cainSpectrum.centersGeV.size());
    EXPECT_NEAR(LinearComptonBenchmark::histogramArea(opalxSpectrum), 1.0, 2.0e-2);
    EXPECT_NEAR(LinearComptonBenchmark::histogramMeanEnergyGeV(opalxSpectrum),
                LinearComptonBenchmark::histogramMeanEnergyGeV(cainSpectrum),
                LinearComptonBenchmark::histogramMeanEnergyGeV(cainSpectrum) * 2.0e-2);
    EXPECT_LT(LinearComptonBenchmark::histogramL1Distance(opalxSpectrum, cainSpectrum), 0.06);
}

TEST(TestLinearComptonSpectrum, FiniteBeamOverlapAngularSpectrumMatchesCainReference) {
    const auto output = runFiniteBeamBenchmark(
        "opalx-linear-compton-90deg-xi029-finite-beam-overlap-theta-histogram.csv",
        true,
        false,
        "--overlap-weighting --beam-sigma-longitudinal 0.000299792 --laser-rayleigh 12.5 --laser-sigma-t 0.000299792");
    const auto opalxSpectrum = LinearComptonBenchmark::readAngleCSV(output);
    const auto cainSpectrum = LinearComptonBenchmark::readAngleCSV(referenceFiniteBeamOverlapAngularSpectrumPath());

    ASSERT_EQ(opalxSpectrum.centersRad.size(), cainSpectrum.centersRad.size());
    EXPECT_NEAR(LinearComptonBenchmark::angleHistogramArea(opalxSpectrum), 1.0, 2.5e-2);
    EXPECT_NEAR(LinearComptonBenchmark::angleHistogramMeanRad(opalxSpectrum),
                LinearComptonBenchmark::angleHistogramMeanRad(cainSpectrum),
                LinearComptonBenchmark::angleHistogramMeanRad(cainSpectrum) * 2.0e-2);
    EXPECT_LT(LinearComptonBenchmark::angleHistogramL1Distance(opalxSpectrum, cainSpectrum), 0.05);
}
