//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#include "AbsBeamline/MultipoleT.h"
#include "Attributes/Attributes.h"
#include "Elements/OpalMultipoleT.h"
#include "gtest/gtest.h"

class TestOpalMultipoleT : public testing::Test {
public:
    TestOpalMultipoleT() = default;

    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
        // Many OPAL writers assume `gmsg` is initialized (see SDDSWriter/StatWriter).
        // Unit tests normally don't set this up via Main().
        gmsg = new Inform(nullptr, -1);
    }
    static void TearDownTestSuite() {
        delete gmsg;
        gmsg = nullptr;
        ippl::finalize();
    }

};

// Does the user interface affect the correct magnet configuration
TEST_F(TestOpalMultipoleT, UserInterface) {
    // Make the UI
    OpalMultipoleT ui;
    // Set the attributes
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::LENGTH], 4.1);
    Attributes::setRealArray(ui.itsAttr[OpalMultipoleT::TP], {0.2, 0.3});
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::LFRINGE], 0.5);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::RFRINGE], 0.6);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::HAPERT], 1.0);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::VAPERT], 1.1);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::MAXFORDER], 4.0);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::ROTATION], 0.0);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::EANGLE], 0.01);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::BBLENGTH], 6.0);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::ANGLE], 0.628);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::MAXXORDER], 7.0);
    Attributes::setBool(ui.itsAttr[OpalMultipoleT::VARRADIUS], false);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::ENTRYOFFSET], 0.0);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::MISALIGN_H], 1.3);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::MISALIGN_V], 1.4);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::MISALIGN_S], 1.5);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::MISALIGN_ROLL], 1.6);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::MISALIGN_YAW], 1.7);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::MISALIGN_PITCH], 1.8);
    // Update the magnet
    EXPECT_NO_THROW(ui.update());
    // Check the values
    auto* myMagnet = dynamic_cast<MultipoleT*>(ui.getElement());
    EXPECT_TRUE(myMagnet);
    EXPECT_NEAR(myMagnet->getElementLength(), 4.1, 1e-6);
    auto tp = myMagnet->getTransProfile();
    EXPECT_EQ(tp.size(), MultipoleTConfig::NumPoles);
    EXPECT_NEAR(tp[0], 0.2, 1e-6);
    EXPECT_NEAR(tp[1], 0.3, 1e-6);
    auto [s0, left, right] = myMagnet->getFringeField();
    EXPECT_NEAR(s0, 4.1 / 2.0, 1e-6);
    EXPECT_NEAR(left, 0.5, 1e-6);
    EXPECT_NEAR(right, 0.6, 1e-6);
    auto [vertical, horizontal] = myMagnet->getAperture();
    EXPECT_NEAR(vertical, 1.1, 1e-6);
    EXPECT_NEAR(horizontal, 1.0, 1e-6);
    EXPECT_NEAR(myMagnet->getMaxFOrder(), 4.0, 1e-6);
    EXPECT_NEAR(myMagnet->getRotation(), 0.0, 1e-6);
    EXPECT_NEAR(myMagnet->getEntranceAngle(), 0.01, 1e-6);
    EXPECT_NEAR(myMagnet->getBoundingBoxLength(), 6.0, 1e-6);
    EXPECT_NEAR(myMagnet->getBendAngle(), 0.628, 1e-6);
    EXPECT_NEAR(myMagnet->getMaxXOrder(), 7.0, 1e-6);
    EXPECT_FALSE(myMagnet->getVariableRadius());
    EXPECT_NEAR(myMagnet->getEntryOffset(), 0.0, 1e-6);
    // Check rotation (only works for straight magnets)
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::ANGLE], 0.0);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::ROTATION], 0.1);
    EXPECT_NO_THROW(ui.update());
    EXPECT_NEAR(myMagnet->getRotation(), 0.1, 1e-6);
    // Check entry offset (only works for var radius magnets)
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::ROTATION], 0.0);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::ANGLE], 0.5);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::ENTRYOFFSET], 1.2);
    Attributes::setBool(ui.itsAttr[OpalMultipoleT::VARRADIUS], true);
    EXPECT_ANY_THROW(ui.update());   // Currently throws 'cos var radius not yet supported
    // Check time dependency
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::ENTRYOFFSET], 0.0);
    Attributes::setBool(ui.itsAttr[OpalMultipoleT::VARRADIUS], false);
    Attributes::setUpperCaseString(ui.itsAttr[OpalMultipoleT::SCALING_MODEL], "Scaling");
    EXPECT_NO_THROW(ui.update());
    EXPECT_EQ(myMagnet->getScalingName(), "SCALING");
}

// Does the user interface clone have the correct magnet configuration
TEST_F(TestOpalMultipoleT, UserInterfaceClone) {
    // Make the UI
    OpalMultipoleT ui;
    // Set the attributes
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::LENGTH], 4.1);
    Attributes::setRealArray(ui.itsAttr[OpalMultipoleT::TP], {0.2, 0.3});
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::LFRINGE], 0.5);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::RFRINGE], 0.6);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::HAPERT], 1.0);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::VAPERT], 1.1);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::MAXFORDER], 4.0);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::ROTATION], 0.0);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::EANGLE], 0.01);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::BBLENGTH], 6.0);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::ANGLE], 0.628);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::MAXXORDER], 7.0);
    Attributes::setBool(ui.itsAttr[OpalMultipoleT::VARRADIUS], false);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::ENTRYOFFSET], 0.0);
    Attributes::setUpperCaseString(ui.itsAttr[OpalMultipoleT::SCALING_MODEL], "Scaling");
    // Make the clone
    std::unique_ptr<OpalMultipoleT> uiClone{ui.clone("Clone")};
    // Update the magnet
    EXPECT_NO_THROW(uiClone->update());
    // Check the values
    auto* myMagnet = dynamic_cast<MultipoleT*>(uiClone->getElement());
    EXPECT_TRUE(myMagnet);
    EXPECT_NEAR(myMagnet->getElementLength(), 4.1, 1e-6);
    auto tp = myMagnet->getTransProfile();
    EXPECT_EQ(tp.size(), 6);
    EXPECT_NEAR(tp[0], 0.2, 1e-6);
    EXPECT_NEAR(tp[1], 0.3, 1e-6);
    EXPECT_NEAR(tp[2], 0.0, 1e-6);
    EXPECT_NEAR(tp[3], 0.0, 1e-6);
    EXPECT_NEAR(tp[4], 0.0, 1e-6);
    EXPECT_NEAR(tp[5], 0.0, 1e-6);
    auto [s0, left, right] = myMagnet->getFringeField();
    EXPECT_NEAR(s0, 4.1 / 2.0, 1e-6);
    EXPECT_NEAR(left, 0.5, 1e-6);
    EXPECT_NEAR(right, 0.6, 1e-6);
    auto [vertical, horizontal] = myMagnet->getAperture();
    EXPECT_NEAR(vertical, 1.1, 1e-6);
    EXPECT_NEAR(horizontal, 1.0, 1e-6);
    EXPECT_NEAR(myMagnet->getMaxFOrder(), 4.0, 1e-6);
    EXPECT_NEAR(myMagnet->getRotation(), 0.0, 1e-6);
    EXPECT_NEAR(myMagnet->getEntranceAngle(), 0.01, 1e-6);
    EXPECT_NEAR(myMagnet->getBoundingBoxLength(), 6.0, 1e-6);
    EXPECT_NEAR(myMagnet->getBendAngle(), 0.628, 1e-6);
    EXPECT_NEAR(myMagnet->getMaxXOrder(), 7.0, 1e-6);
    EXPECT_FALSE(myMagnet->getVariableRadius());
    EXPECT_NEAR(myMagnet->getEntryOffset(), 0.0, 1e-6);
    EXPECT_EQ(myMagnet->getScalingName(), "SCALING");
}

// Check the various configuration validation tests
TEST_F(TestOpalMultipoleT, UserInterfaceSanityCheck) {
    // Make the UI
    OpalMultipoleT ui;
    // Set the attributes
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::LENGTH], 4.1);
    Attributes::setRealArray(ui.itsAttr[OpalMultipoleT::TP], {0.2, 0.3});
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::LFRINGE], 0.5);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::RFRINGE], 0.6);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::HAPERT], 1.0);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::VAPERT], 1.1);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::ROTATION], 0.0);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::EANGLE], 0.0);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::BBLENGTH], 6.0);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::ANGLE], 0.0);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::MAXXORDER], 7.0);
    Attributes::setBool(ui.itsAttr[OpalMultipoleT::VARRADIUS], false);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::ENTRYOFFSET], 0.0);
    // Try F order of zero
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::MAXFORDER], 0.0);
    EXPECT_ANY_THROW(ui.update());
    // Try F order of 21
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::MAXFORDER], 21.0);
    EXPECT_ANY_THROW(ui.update());
    // Try skew straight magnet
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::MAXFORDER], 4);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::ANGLE], 0.0);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::ROTATION], 0.5);
    EXPECT_NO_THROW(ui.update());
    // Try skew bent magnet
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::ANGLE], 0.5);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::ROTATION], 0.5);
    EXPECT_ANY_THROW(ui.update());
    // Try the variable radius magnet
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::ROTATION], 0.0);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::ANGLE], 0.5);
    Attributes::setBool(ui.itsAttr[OpalMultipoleT::VARRADIUS], true);
    EXPECT_ANY_THROW(ui.update());
    // Try the entry offset with a straight magnet
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::ANGLE], 0.0);
    Attributes::setBool(ui.itsAttr[OpalMultipoleT::VARRADIUS], false);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::ENTRYOFFSET], 1.0);
    EXPECT_ANY_THROW(ui.update());
    // Try the entry offset with a bent magnet
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::ANGLE], 0.5);
    Attributes::setBool(ui.itsAttr[OpalMultipoleT::VARRADIUS], false);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::ENTRYOFFSET], 1.0);
    EXPECT_ANY_THROW(ui.update());
    // Try the entry offset with a variable radius bent magnet
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::ANGLE], 0.5);
    Attributes::setBool(ui.itsAttr[OpalMultipoleT::VARRADIUS], true);
    Attributes::setReal(ui.itsAttr[OpalMultipoleT::ENTRYOFFSET], 1.0);
    EXPECT_ANY_THROW(ui.update());
}

// Does the print API work
TEST_F(TestOpalMultipoleT, Print) {
    // Make the UI
    const OpalMultipoleT ui;
    // And print it
    std::stringstream ss;
    ui.print(ss);
    EXPECT_EQ(ss.str(), "MULTIPOLET;\n");
}
