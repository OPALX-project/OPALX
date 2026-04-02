#include "Attributes/Attributes.h"
#include "BeamlineCore/LaserRep.h"
#include "Elements/OpalLaser.h"
#include "Utilities/OpalException.h"

#include "gtest/gtest.h"

namespace {
    void setValidLaserAttributes(OpalLaser& laser) {
        Attributes::setReal(laser.itsAttr[OpalLaser::WAVELENGTH], 1.03e-6);
        Attributes::setReal(laser.itsAttr[OpalLaser::PULSEENERGY], 1.0);
        Attributes::setReal(laser.itsAttr[OpalLaser::PULSELENGTH], 2.0e-12);
        Attributes::setReal(laser.itsAttr[OpalLaser::WAISTX], 5.0e-6);
        Attributes::setReal(laser.itsAttr[OpalLaser::WAISTY], 6.0e-6);
        Attributes::setRealArray(laser.itsAttr[OpalLaser::DIR], {0.0, 0.0, -2.0});
    }
}

TEST(TestLaser, UpdateStoresValidatedParameters) {
    OpalLaser laser;
    setValidLaserAttributes(laser);
    Attributes::setReal(laser.itsAttr[OpalElement::LENGTH], 0.01);
    Attributes::setRealArray(laser.itsAttr[OpalLaser::STOKES], {0.0, 0.0, 1.0});

    EXPECT_NO_THROW(laser.update());

    auto* rep = dynamic_cast<LaserRep*>(laser.getElement());
    ASSERT_NE(rep, nullptr);
    EXPECT_DOUBLE_EQ(rep->getElementLength(), 0.01);
    EXPECT_DOUBLE_EQ(rep->getWavelength(), 1.03e-6);
    EXPECT_DOUBLE_EQ(rep->getPulseEnergy(), 1.0);
    EXPECT_DOUBLE_EQ(rep->getPulseLength(), 2.0e-12);
    EXPECT_DOUBLE_EQ(rep->getWaistX(), 5.0e-6);
    EXPECT_DOUBLE_EQ(rep->getWaistY(), 6.0e-6);
    EXPECT_NEAR(rep->getDirection()(0), 0.0, 1.0e-15);
    EXPECT_NEAR(rep->getDirection()(1), 0.0, 1.0e-15);
    EXPECT_NEAR(rep->getDirection()(2), -1.0, 1.0e-15);
    EXPECT_DOUBLE_EQ(rep->getStokes()(2), 1.0);
}

TEST(TestLaser, UpdateDefaultsStokesToZero) {
    OpalLaser laser;
    setValidLaserAttributes(laser);

    EXPECT_NO_THROW(laser.update());

    auto* rep = dynamic_cast<LaserRep*>(laser.getElement());
    ASSERT_NE(rep, nullptr);
    EXPECT_DOUBLE_EQ(rep->getStokes()(0), 0.0);
    EXPECT_DOUBLE_EQ(rep->getStokes()(1), 0.0);
    EXPECT_DOUBLE_EQ(rep->getStokes()(2), 0.0);
}

TEST(TestLaser, UpdateRequiresWavelength) {
    OpalLaser laser;
    Attributes::setReal(laser.itsAttr[OpalLaser::PULSEENERGY], 1.0);
    Attributes::setReal(laser.itsAttr[OpalLaser::PULSELENGTH], 2.0e-12);
    Attributes::setReal(laser.itsAttr[OpalLaser::WAISTX], 5.0e-6);
    Attributes::setReal(laser.itsAttr[OpalLaser::WAISTY], 6.0e-6);
    Attributes::setRealArray(laser.itsAttr[OpalLaser::DIR], {0.0, 0.0, -2.0});

    EXPECT_THROW(laser.update(), OpalException);
}

TEST(TestLaser, UpdateRejectsNonPositivePulseEnergy) {
    OpalLaser laser;
    setValidLaserAttributes(laser);
    Attributes::setReal(laser.itsAttr[OpalLaser::PULSEENERGY], 0.0);

    EXPECT_THROW(laser.update(), OpalException);
}

TEST(TestLaser, UpdateRejectsBadDirectionSize) {
    OpalLaser laser;
    setValidLaserAttributes(laser);
    Attributes::setRealArray(laser.itsAttr[OpalLaser::DIR], {0.0, 1.0});

    EXPECT_THROW(laser.update(), OpalException);
}

TEST(TestLaser, UpdateRejectsZeroDirection) {
    OpalLaser laser;
    setValidLaserAttributes(laser);
    Attributes::setRealArray(laser.itsAttr[OpalLaser::DIR], {0.0, 0.0, 0.0});

    EXPECT_THROW(laser.update(), OpalException);
}

TEST(TestLaser, UpdateRejectsBadStokesSize) {
    OpalLaser laser;
    setValidLaserAttributes(laser);
    Attributes::setRealArray(laser.itsAttr[OpalLaser::STOKES], {0.0, 1.0});

    EXPECT_THROW(laser.update(), OpalException);
}

TEST(TestLaser, UpdateRejectsStokesNormAboveOne) {
    OpalLaser laser;
    setValidLaserAttributes(laser);
    Attributes::setRealArray(laser.itsAttr[OpalLaser::STOKES], {1.0, 1.0, 0.0});

    EXPECT_THROW(laser.update(), OpalException);
}

TEST(TestLaser, UpdateRejectsNegativeLength) {
    OpalLaser laser;
    setValidLaserAttributes(laser);
    Attributes::setReal(laser.itsAttr[OpalElement::LENGTH], -1.0);

    EXPECT_THROW(laser.update(), OpalException);
}
