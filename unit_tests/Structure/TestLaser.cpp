#include "AbsBeamline/Laser.h"
#include "Attributes/Attributes.h"
#include "Elements/OpalLaser.h"
#include "Physics/Physics.h"
#include "Utilities/LogicalError.h"
#include "Utilities/OpalException.h"

#include "gtest/gtest.h"

#include <cmath>
#include "../KeywordObjects.h"

namespace {
    void setValidLaserAttributes(OpalLaser& laser) {
        Attributes::setReal(laser.itsAttr[OpalLaser::WAVELENGTH], 1.03e-6);
        Attributes::setReal(laser.itsAttr[OpalLaser::PULSEENERGY], 1.0);
        Attributes::setReal(laser.itsAttr[OpalLaser::PULSELENGTH], 2.0e-12);
        Attributes::setReal(laser.itsAttr[OpalLaser::WAISTX], 5.0e-6);
        Attributes::setReal(laser.itsAttr[OpalLaser::WAISTY], 6.0e-6);
        Attributes::setRealArray(laser.itsAttr[OpalLaser::DIR], {0.0, 0.0, -2.0});
    }
}  // namespace

TEST(TestLaser, UpdateStoresValidatedParameters) {
    auto laserObject = opalx::test::makePlacedElement<OpalLaser>("LASER", "LASER");
    OpalLaser& laser = *laserObject;
    setValidLaserAttributes(laser);
    Attributes::setReal(laser.itsAttr[OpalElement::LENGTH], 0.01);
    Attributes::setRealArray(laser.itsAttr[OpalLaser::STOKES], {0.0, 0.0, 1.0});

    std::shared_ptr<ElementBase> built;
    EXPECT_NO_THROW(built = opalx::test::buildFrom(laser));

    auto* element = dynamic_cast<Laser*>(built.get());
    ASSERT_NE(element, nullptr);
    EXPECT_DOUBLE_EQ(element->getGeometry().getElementLength(), 0.01);
    EXPECT_DOUBLE_EQ(element->getWavelength(), 1.03e-6);
    EXPECT_DOUBLE_EQ(element->getPulseEnergy(), 1.0);
    EXPECT_DOUBLE_EQ(element->getPulseLength(), 2.0e-12);
    EXPECT_DOUBLE_EQ(element->getWaistX(), 5.0e-6);
    EXPECT_DOUBLE_EQ(element->getWaistY(), 6.0e-6);
    EXPECT_NEAR(element->getDirection()(0), 0.0, 1.0e-12);
    EXPECT_NEAR(element->getDirection()(1), 0.0, 1.0e-12);
    EXPECT_NEAR(element->getDirection()(2), -1.0, 1.0e-12);
    EXPECT_DOUBLE_EQ(element->getStokes()(2), 1.0);
}

TEST(TestLaser, UpdateDefaultsStokesToZero) {
    auto laserObject = opalx::test::makePlacedElement<OpalLaser>("LASER", "LASER");
    OpalLaser& laser = *laserObject;
    setValidLaserAttributes(laser);

    std::shared_ptr<ElementBase> built;
    EXPECT_NO_THROW(built = opalx::test::buildFrom(laser));

    auto* element = dynamic_cast<Laser*>(built.get());
    ASSERT_NE(element, nullptr);
    EXPECT_DOUBLE_EQ(element->getStokes()(0), 0.0);
    EXPECT_DOUBLE_EQ(element->getStokes()(1), 0.0);
    EXPECT_DOUBLE_EQ(element->getStokes()(2), 0.0);
}

TEST(TestLaser, UpdateRequiresWavelength) {
    auto laserObject = opalx::test::makePlacedElement<OpalLaser>("LASER", "LASER");
    OpalLaser& laser = *laserObject;
    Attributes::setReal(laser.itsAttr[OpalLaser::PULSEENERGY], 1.0);
    Attributes::setReal(laser.itsAttr[OpalLaser::PULSELENGTH], 2.0e-12);
    Attributes::setReal(laser.itsAttr[OpalLaser::WAISTX], 5.0e-6);
    Attributes::setReal(laser.itsAttr[OpalLaser::WAISTY], 6.0e-6);
    Attributes::setRealArray(laser.itsAttr[OpalLaser::DIR], {0.0, 0.0, -2.0});

    EXPECT_THROW(opalx::test::buildFrom(laser), OpalException);
}

TEST(TestLaser, UpdateRejectsNonPositivePulseEnergy) {
    auto laserObject = opalx::test::makePlacedElement<OpalLaser>("LASER", "LASER");
    OpalLaser& laser = *laserObject;
    setValidLaserAttributes(laser);
    Attributes::setReal(laser.itsAttr[OpalLaser::PULSEENERGY], 0.0);

    EXPECT_THROW(opalx::test::buildFrom(laser), OpalException);
}

TEST(TestLaser, UpdateRejectsBadDirectionSize) {
    auto laserObject = opalx::test::makePlacedElement<OpalLaser>("LASER", "LASER");
    OpalLaser& laser = *laserObject;
    setValidLaserAttributes(laser);
    Attributes::setRealArray(laser.itsAttr[OpalLaser::DIR], {0.0, 1.0});

    EXPECT_THROW(opalx::test::buildFrom(laser), OpalException);
}

TEST(TestLaser, UpdateRejectsZeroDirection) {
    auto laserObject = opalx::test::makePlacedElement<OpalLaser>("LASER", "LASER");
    OpalLaser& laser = *laserObject;
    setValidLaserAttributes(laser);
    Attributes::setRealArray(laser.itsAttr[OpalLaser::DIR], {0.0, 0.0, 0.0});

    EXPECT_THROW(opalx::test::buildFrom(laser), OpalException);
}

TEST(TestLaser, UpdateRejectsBadStokesSize) {
    auto laserObject = opalx::test::makePlacedElement<OpalLaser>("LASER", "LASER");
    OpalLaser& laser = *laserObject;
    setValidLaserAttributes(laser);
    Attributes::setRealArray(laser.itsAttr[OpalLaser::STOKES], {0.0, 1.0});

    EXPECT_THROW(opalx::test::buildFrom(laser), OpalException);
}

TEST(TestLaser, UpdateRejectsStokesNormAboveOne) {
    auto laserObject = opalx::test::makePlacedElement<OpalLaser>("LASER", "LASER");
    OpalLaser& laser = *laserObject;
    setValidLaserAttributes(laser);
    Attributes::setRealArray(laser.itsAttr[OpalLaser::STOKES], {1.0, 1.0, 0.0});

    EXPECT_THROW(opalx::test::buildFrom(laser), OpalException);
}

TEST(TestLaser, UpdateRejectsNegativeLength) {
    auto laserObject = opalx::test::makePlacedElement<OpalLaser>("LASER", "LASER");
    OpalLaser& laser = *laserObject;
    setValidLaserAttributes(laser);
    Attributes::setReal(laser.itsAttr[OpalElement::LENGTH], -1.0);

    EXPECT_THROW(opalx::test::buildFrom(laser), OpalException);
}

/**
 * @brief Benchmark the CAIN-inspired linear Compton helper for a 90 degree laser crossing.
 *
 * The configuration is the simplest fixed-angle case used as a physics-facing unit test:
 * an electron with total energy \f$E_e\f$ moves along \f$+\hat z\f$, the laser propagates
 * along \f$+\hat x\f$, and the scattered photon is evaluated in the forward electron
 * direction \f$+\hat z\f$.
 *
 * The exact reference formulas used in this test are
 * \f[
 *   \omega_L = \frac{2\pi \hbar c}{\lambda_L},
 * \f]
 * \f[
 *   p_e = \sqrt{E_e^2 - m_e^2},
 * \f]
 * \f[
 *   x = \frac{2 E_e \omega_L}{m_e^2},
 * \f]
 * and the forward scattered photon energy from exact energy-momentum conservation,
 * written in a numerically stable form,
 * \f[
 *   \omega_\gamma^{\prime}
 *   = \frac{E_e\,\omega_L}{\omega_L + m_e^2/(E_e + p_e)}.
 * \f]
 *
 * The stable denominator avoids catastrophic cancellation in the equivalent expression
 * \f$E_e - p_e + \omega_L\f$ for relativistic electrons.
 */
TEST(TestLaser, LinearComptonForwardPhotonEnergyMatchesExactNinetyDegreeKinematics) {
    auto laserObject = opalx::test::makePlacedElement<OpalLaser>("LASER", "LASER");
    OpalLaser& laser = *laserObject;
    setValidLaserAttributes(laser);
    Attributes::setRealArray(laser.itsAttr[OpalLaser::DIR], {1.0, 0.0, 0.0});

    std::shared_ptr<ElementBase> built;
    ASSERT_NO_THROW(built = opalx::test::buildFrom(laser));

    auto* element = dynamic_cast<Laser*>(built.get());
    ASSERT_NE(element, nullptr);

    Vector_t<double, 3> beamDirection(0.0);
    beamDirection(2) = 1.0;

    const double electronTotalEnergyGeV = 1.0;
    const double laserPhotonEnergyGeV   = element->getPhotonEnergyGeV();
    const double electronMomentumGeV    = std::sqrt(
            electronTotalEnergyGeV * electronTotalEnergyGeV - Physics::m_e * Physics::m_e);
    const double expectedInvariantX =
            2.0 * laserPhotonEnergyGeV * electronTotalEnergyGeV / (Physics::m_e * Physics::m_e);
    const double stableEnergyMinusMomentum =
            Physics::m_e * Physics::m_e / (electronTotalEnergyGeV + electronMomentumGeV);
    const double expectedForwardPhotonEnergyGeV =
            laserPhotonEnergyGeV * electronTotalEnergyGeV
            / (stableEnergyMinusMomentum + laserPhotonEnergyGeV);

    EXPECT_NEAR(
            element->getLinearComptonInvariantX(electronTotalEnergyGeV, beamDirection),
            expectedInvariantX, expectedInvariantX * 1.0e-12);
    EXPECT_NEAR(
            element->getLinearComptonForwardPhotonEnergyGeV(electronTotalEnergyGeV, beamDirection),
            expectedForwardPhotonEnergyGeV, expectedForwardPhotonEnergyGeV * 1.0e-12);
}
