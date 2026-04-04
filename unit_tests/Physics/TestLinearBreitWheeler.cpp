/**
 * \file TestLinearBreitWheeler.cpp
 * \brief Unit tests for the first host-only linear Breit-Wheeler helper.
 *
 * The first OPALX implementation is intentionally limited to unpolarized,
 * fixed-geometry, host-side event generation. The tests below validate:
 *
 * - the two-photon invariant and threshold,
 * - the standard total linear Breit-Wheeler cross section,
 * - event-level energy-momentum conservation,
 * - deterministic sampling when the host RNG is seeded from `Options::seed`.
 */

#include "Physics/LinearBreitWheeler.h"
#include "Physics/Physics.h"
#include "Utilities/Options.h"
#include "gtest/gtest.h"

#include <cmath>

namespace {
    Vector_t<double, 3> makeDirection(double x, double y, double z) {
        Vector_t<double, 3> value(0.0);
        value(0) = x;
        value(1) = y;
        value(2) = z;
        return value;
    }
}

TEST(TestLinearBreitWheeler, HeadOnThresholdMatchesAnalyticCondition) {
    const double wavelength_m = 1.0e-9;
    const double laserPhotonEnergyGeV = Physics::LinearBreitWheeler::photonEnergyFromWavelengthGeV(wavelength_m);
    const double expectedThresholdPhotonEnergyGeV = Physics::m_e * Physics::m_e / laserPhotonEnergyGeV;

    const auto n1 = makeDirection(0.0, 0.0, 1.0);
    const auto n2 = makeDirection(0.0, 0.0, -1.0);
    const double thresholdS = Physics::LinearBreitWheeler::thresholdInvariantSGeV2();
    const double computedS = Physics::LinearBreitWheeler::invariantSGeV2(expectedThresholdPhotonEnergyGeV,
                                                                         laserPhotonEnergyGeV,
                                                                         n1,
                                                                         n2);

    EXPECT_NEAR(computedS, thresholdS, thresholdS * 1.0e-12);
    EXPECT_TRUE(Physics::LinearBreitWheeler::isAboveThreshold(computedS));
}

TEST(TestLinearBreitWheeler, TotalCrossSectionVanishesBelowThreshold) {
    const double belowThresholdS = 0.99 * Physics::LinearBreitWheeler::thresholdInvariantSGeV2();
    EXPECT_DOUBLE_EQ(Physics::LinearBreitWheeler::totalCrossSection(belowThresholdS), 0.0);
}

TEST(TestLinearBreitWheeler, TotalCrossSectionIsPositiveAboveThreshold) {
    const double aboveThresholdS = 16.0 * Physics::m_e * Physics::m_e;
    const double sigma = Physics::LinearBreitWheeler::totalCrossSection(aboveThresholdS);
    EXPECT_GT(sigma, 0.0);
}

TEST(TestLinearBreitWheeler, SampledEventConservesFourMomentumForHeadOnGeometry) {
    const int previousSeed = Options::seed;
    Options::seed = 20260404;

    const double wavelength_m = 1.0e-9;
    const double laserPhotonEnergyGeV = Physics::LinearBreitWheeler::photonEnergyFromWavelengthGeV(wavelength_m);
    const double highEnergyPhotonEnergyGeV = 0.5;
    const auto highEnergyDirection = makeDirection(0.0, 0.0, 1.0);
    const auto laserDirection = makeDirection(0.0, 0.0, -1.0);

    const auto kernel = Physics::LinearBreitWheeler::makeSamplingKernel(highEnergyPhotonEnergyGeV,
                                                                        laserPhotonEnergyGeV,
                                                                        highEnergyDirection,
                                                                        laserDirection);
    auto engine = Physics::LinearBreitWheeler::makeHostRandomEngine();
    const auto event = Physics::LinearBreitWheeler::sampleEvent(kernel, engine);

    Vector_t<double, 3> incomingMomentum(0.0);
    incomingMomentum += highEnergyPhotonEnergyGeV * highEnergyDirection;
    incomingMomentum += laserPhotonEnergyGeV * laserDirection;

    const Vector_t<double, 3> outgoingMomentum = event.electronMomentumLabGeV + event.positronMomentumLabGeV;
    const double incomingEnergy = highEnergyPhotonEnergyGeV + laserPhotonEnergyGeV;
    const double outgoingEnergy = event.electronEnergyLabGeV + event.positronEnergyLabGeV;

    EXPECT_NEAR(outgoingEnergy, incomingEnergy, incomingEnergy * 1.0e-10);
    EXPECT_NEAR(outgoingMomentum(0), incomingMomentum(0), 1.0e-10);
    EXPECT_NEAR(outgoingMomentum(1), incomingMomentum(1), 1.0e-10);
    EXPECT_NEAR(outgoingMomentum(2), incomingMomentum(2), 1.0e-10);
    EXPECT_GE(event.scatteringCosineCM, -1.0);
    EXPECT_LE(event.scatteringCosineCM, 1.0);
    EXPECT_GE(event.azimuthCM, 0.0);
    EXPECT_LT(event.azimuthCM, Physics::two_pi);

    Options::seed = previousSeed;
}

TEST(TestLinearBreitWheeler, FixedSeedProducesDeterministicSampleSequence) {
    const int previousSeed = Options::seed;
    Options::seed = 424242;

    const double wavelength_m = 1.0e-9;
    const double laserPhotonEnergyGeV = Physics::LinearBreitWheeler::photonEnergyFromWavelengthGeV(wavelength_m);
    const double highEnergyPhotonEnergyGeV = 0.5;
    const auto highEnergyDirection = makeDirection(0.0, 0.0, 1.0);
    const auto laserDirection = makeDirection(0.0, 0.0, -1.0);

    const auto kernel = Physics::LinearBreitWheeler::makeSamplingKernel(highEnergyPhotonEnergyGeV,
                                                                        laserPhotonEnergyGeV,
                                                                        highEnergyDirection,
                                                                        laserDirection);
    auto engine1 = Physics::LinearBreitWheeler::makeHostRandomEngine();
    auto engine2 = Physics::LinearBreitWheeler::makeHostRandomEngine();

    for (int i = 0; i < 8; ++i) {
        const auto event1 = Physics::LinearBreitWheeler::sampleEvent(kernel, engine1);
        const auto event2 = Physics::LinearBreitWheeler::sampleEvent(kernel, engine2);
        EXPECT_DOUBLE_EQ(event1.scatteringCosineCM, event2.scatteringCosineCM);
        EXPECT_DOUBLE_EQ(event1.azimuthCM, event2.azimuthCM);
        EXPECT_DOUBLE_EQ(event1.electronEnergyLabGeV, event2.electronEnergyLabGeV);
        EXPECT_DOUBLE_EQ(event1.positronEnergyLabGeV, event2.positronEnergyLabGeV);
    }

    Options::seed = previousSeed;
}
