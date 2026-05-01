/**
 * \file TestLinearCompton.cpp
 * \brief Unit tests for the unpolarized linear Compton helper.
 *
 * The implemented helper is intentionally limited to single-electron,
 * unpolarized Klein-Nishina physics. The tests below validate four things:
 *
 * \test TestLinearCompton::ThomsonLimitRecoveredAtLowEnergy
 * Verifies that the exact Klein-Nishina total cross section approaches the
 * Thomson limit for @f$\omega_1 \ll m_e@f$.
 *
 * \test TestLinearCompton::EnergySupportMatchesComptonEndpoints
 * Verifies the exact electron-rest-frame energy support
 * @f$\omega_2 \in [\omega_1/(1+2\omega_1/m_e), \omega_1]@f$ and checks the
 * recoil formula at the forward and backward endpoints.
 *
 * \test TestLinearCompton::DifferentialSpectrumIntegratesToTotalCrossSection
 * Numerically integrates @f$d\sigma/d\omega_2@f$ over the exact support and
 * compares it against the analytic Klein-Nishina total cross section.
 *
 * \test TestLinearCompton::NinetyDegreeLabForwardBenchmarkMatchesExactFormula
 * Reproduces the same 90 degree crossing benchmark used by the laser element:
 * a 1 GeV electron moving along @f$+\hat z@f$ and a 1030 nm laser propagating
 * along @f$+\hat x@f$. The forward scattered photon energy is compared against
 * the exact stable laboratory-frame formula
 * @f[
 *   \omega_\gamma' =
 *   \frac{E_e \omega_L}{\omega_L + m_e^2 / (E_e + p_e)}.
 * @f]
 */

#include "Physics/LinearCompton.h"
#include "Physics/Physics.h"
#include "Utilities/Options.h"
#include "gtest/gtest.h"

#include <cmath>

namespace {
    double integrateSpectrumCompositeSimpson(double incomingPhotonEnergyERFGeV,
                                             double minEnergy,
                                             double maxEnergy,
                                             int panels) {
        const double h = (maxEnergy - minEnergy) / panels;
        double sum = Physics::LinearCompton::differentialCrossSectionOmegaERF(incomingPhotonEnergyERFGeV,
                                                                               minEnergy)
            + Physics::LinearCompton::differentialCrossSectionOmegaERF(incomingPhotonEnergyERFGeV,
                                                                       maxEnergy);

        for (int i = 1; i < panels; ++i) {
            const double energy = minEnergy + i * h;
            const double weight = (i % 2 == 0) ? 2.0 : 4.0;
            sum += weight * Physics::LinearCompton::differentialCrossSectionOmegaERF(
                incomingPhotonEnergyERFGeV,
                energy);
        }

        return sum * h / 3.0;
    }
}

TEST(TestLinearCompton, ThomsonLimitRecoveredAtLowEnergy) {
    const double incomingPhotonEnergyERFGeV = 1.0e-9;
    const double sigmaKN = Physics::LinearCompton::totalCrossSection(incomingPhotonEnergyERFGeV);
    const double sigmaThomson = Physics::LinearCompton::thomsonCrossSection();

    EXPECT_NEAR(sigmaKN / sigmaThomson, 1.0, 1.0e-5);
}

TEST(TestLinearCompton, EnergySupportMatchesComptonEndpoints) {
    const double incomingPhotonEnergyERFGeV = 2.5e-4;
    const double minEnergy = Physics::LinearCompton::scatteredPhotonEnergyMinERFGeV(
        incomingPhotonEnergyERFGeV);
    const double maxEnergy = Physics::LinearCompton::scatteredPhotonEnergyMaxERFGeV(
        incomingPhotonEnergyERFGeV);

    EXPECT_DOUBLE_EQ(maxEnergy, incomingPhotonEnergyERFGeV);
    EXPECT_NEAR(Physics::LinearCompton::scatteredPhotonEnergyERFGeV(incomingPhotonEnergyERFGeV, 1.0),
                maxEnergy,
                1.0e-15);
    EXPECT_NEAR(Physics::LinearCompton::scatteredPhotonEnergyERFGeV(incomingPhotonEnergyERFGeV, -1.0),
                minEnergy,
                1.0e-15);
    EXPECT_LT(minEnergy, maxEnergy);
}

TEST(TestLinearCompton, DifferentialSpectrumIntegratesToTotalCrossSection) {
    const double incomingPhotonEnergyERFGeV = 2.0e-4;
    const double minEnergy = Physics::LinearCompton::scatteredPhotonEnergyMinERFGeV(
        incomingPhotonEnergyERFGeV);
    const double maxEnergy = Physics::LinearCompton::scatteredPhotonEnergyMaxERFGeV(
        incomingPhotonEnergyERFGeV);
    const double integral = integrateSpectrumCompositeSimpson(incomingPhotonEnergyERFGeV,
                                                              minEnergy,
                                                              maxEnergy,
                                                              20000);

    const double sigmaKN = Physics::LinearCompton::totalCrossSection(incomingPhotonEnergyERFGeV);
    EXPECT_NEAR(integral, sigmaKN, sigmaKN * 2.0e-6);
}

TEST(TestLinearCompton, NinetyDegreeLabForwardBenchmarkMatchesExactFormula) {
    const double electronTotalEnergyGeV = 1.0;
    const double wavelength_m = 1.03e-6;
    const double laserPhotonEnergyGeV = Physics::LinearCompton::photonEnergyFromWavelengthGeV(wavelength_m);

    Vector_t<double, 3> beamDirection(0.0);
    beamDirection(2) = 1.0;
    Vector_t<double, 3> laserDirection(0.0);
    laserDirection(0) = 1.0;

    const double electronMomentumGeV = std::sqrt(
        electronTotalEnergyGeV * electronTotalEnergyGeV - Physics::m_e * Physics::m_e);
    const double expectedInvariantX =
        2.0 * laserPhotonEnergyGeV * electronTotalEnergyGeV / (Physics::m_e * Physics::m_e);
    const double stableEnergyMinusMomentum =
        Physics::m_e * Physics::m_e / (electronTotalEnergyGeV + electronMomentumGeV);
    const double expectedForwardPhotonEnergyGeV =
        laserPhotonEnergyGeV * electronTotalEnergyGeV
        / (stableEnergyMinusMomentum + laserPhotonEnergyGeV);

    EXPECT_NEAR(Physics::LinearCompton::invariantX(electronTotalEnergyGeV,
                                                   laserPhotonEnergyGeV,
                                                   beamDirection,
                                                   laserDirection),
                expectedInvariantX,
                expectedInvariantX * 1.0e-12);
    EXPECT_NEAR(Physics::LinearCompton::labForwardPhotonEnergyGeV(electronTotalEnergyGeV,
                                                                  laserPhotonEnergyGeV,
                                                                  beamDirection,
                                                                  laserDirection),
                expectedForwardPhotonEnergyGeV,
                expectedForwardPhotonEnergyGeV * 1.0e-12);
}


TEST(TestLinearCompton, ExplicitLabPhotonEnergyMatchesForwardHelper) {
    const double electronTotalEnergyGeV = 1.0;
    const double wavelength_m = 1.03e-6;
    const double laserPhotonEnergyGeV = Physics::LinearCompton::photonEnergyFromWavelengthGeV(wavelength_m);

    Vector_t<double, 3> beamDirection(0.0);
    beamDirection(2) = 1.0;
    Vector_t<double, 3> laserDirection(0.0);
    laserDirection(0) = 1.0;

    const double scatteringCosine =
        Physics::LinearCompton::restFrameScatteringCosineForLabForwardPhoton(electronTotalEnergyGeV,
                                                                             beamDirection,
                                                                             laserDirection);

    EXPECT_NEAR(Physics::LinearCompton::labPhotonEnergyGeV(electronTotalEnergyGeV,
                                                           laserPhotonEnergyGeV,
                                                           beamDirection,
                                                           laserDirection,
                                                           scatteringCosine,
                                                           0.0),
                Physics::LinearCompton::labForwardPhotonEnergyGeV(electronTotalEnergyGeV,
                                                                  laserPhotonEnergyGeV,
                                                                  beamDirection,
                                                                  laserDirection),
                1.0e-12);
}


TEST(TestLinearCompton, SampledEventsRespectKinematicBounds) {
    const int previousSeed = Options::seed;
    Options::seed = 20260403;

    const double electronTotalEnergyGeV = 1.0;
    const double wavelength_m = 1.03e-6;
    const double laserPhotonEnergyGeV = Physics::LinearCompton::photonEnergyFromWavelengthGeV(wavelength_m);

    Vector_t<double, 3> beamDirection(0.0);
    beamDirection(2) = 1.0;
    Vector_t<double, 3> laserDirection(0.0);
    laserDirection(0) = 1.0;

    const auto kernel = Physics::LinearCompton::makeSamplingKernel(electronTotalEnergyGeV,
                                                                   laserPhotonEnergyGeV,
                                                                   beamDirection,
                                                                   laserDirection);
    auto engine = Physics::LinearCompton::makeHostRandomEngine();
    const auto event = Physics::LinearCompton::sampleEvent(kernel, engine);

    const double minERF = Physics::LinearCompton::scatteredPhotonEnergyMinERFGeV(
        kernel.incomingPhotonEnergyERFGeV);
    const double maxERF = Physics::LinearCompton::scatteredPhotonEnergyMaxERFGeV(
        kernel.incomingPhotonEnergyERFGeV);
    const double directionNorm = std::sqrt(dot(event.scatteredPhotonDirectionLab,
                                               event.scatteredPhotonDirectionLab));

    EXPECT_GE(event.scatteringCosineERF, -1.0);
    EXPECT_LE(event.scatteringCosineERF, 1.0);
    EXPECT_GE(event.azimuthERF, 0.0);
    EXPECT_LT(event.azimuthERF, Physics::two_pi);
    EXPECT_GE(event.scatteredPhotonEnergyERFGeV, minERF);
    EXPECT_LE(event.scatteredPhotonEnergyERFGeV, maxERF);
    EXPECT_GT(event.scatteredPhotonEnergyLabGeV, 0.0);
    EXPECT_NEAR(directionNorm, 1.0, 1.0e-12);

    Options::seed = previousSeed;
}

TEST(TestLinearCompton, FixedSeedProducesDeterministicSampleSequence) {
    const int previousSeed = Options::seed;
    Options::seed = 424242;

    const double electronTotalEnergyGeV = 1.0;
    const double wavelength_m = 1.03e-6;
    const double laserPhotonEnergyGeV = Physics::LinearCompton::photonEnergyFromWavelengthGeV(wavelength_m);

    Vector_t<double, 3> beamDirection(0.0);
    beamDirection(2) = 1.0;
    Vector_t<double, 3> laserDirection(0.0);
    laserDirection(0) = 1.0;

    const auto kernel = Physics::LinearCompton::makeSamplingKernel(electronTotalEnergyGeV,
                                                                   laserPhotonEnergyGeV,
                                                                   beamDirection,
                                                                   laserDirection);
    auto engine1 = Physics::LinearCompton::makeHostRandomEngine();
    auto engine2 = Physics::LinearCompton::makeHostRandomEngine();

    for (int i = 0; i < 8; ++i) {
        const auto event1 = Physics::LinearCompton::sampleEvent(kernel, engine1);
        const auto event2 = Physics::LinearCompton::sampleEvent(kernel, engine2);
        EXPECT_DOUBLE_EQ(event1.scatteringCosineERF, event2.scatteringCosineERF);
        EXPECT_DOUBLE_EQ(event1.azimuthERF, event2.azimuthERF);
        EXPECT_DOUBLE_EQ(event1.scatteredPhotonEnergyERFGeV, event2.scatteredPhotonEnergyERFGeV);
        EXPECT_DOUBLE_EQ(event1.scatteredPhotonEnergyLabGeV, event2.scatteredPhotonEnergyLabGeV);
        for (int axis = 0; axis < 3; ++axis) {
            EXPECT_DOUBLE_EQ(event1.scatteredPhotonDirectionLab(axis),
                             event2.scatteredPhotonDirectionLab(axis));
        }
    }

    Options::seed = previousSeed;
}
