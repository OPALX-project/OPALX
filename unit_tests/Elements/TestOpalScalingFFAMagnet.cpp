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

#include "AbsBeamline/EndFieldModel/EndFieldModelManager.h"
#include "AbsBeamline/EndFieldModel/Tanh.h"
#include "AbsBeamline/ScalingFFAMagnet.h"
#include "Attributes/Attributes.h"
#include "Elements/OpalScalingFFAMagnet.h"
#include "gtest/gtest.h"

class TestOpalScalingFFAMagnet : public testing::Test {
public:
    TestOpalScalingFFAMagnet() = default;

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
TEST_F(TestOpalScalingFFAMagnet, UserInterface) {
    // Make the UI
    OpalScalingFFAMagnet ui;
    // Set the attributes
    Attributes::setReal(ui.itsAttr[OpalElement::LENGTH], 1);
    Attributes::setReal(ui.itsAttr[OpalScalingFFAMagnet::B0], 2);
    Attributes::setReal(ui.itsAttr[OpalScalingFFAMagnet::R0], 19);
    Attributes::setReal(ui.itsAttr[OpalScalingFFAMagnet::FIELD_INDEX], 4);
    Attributes::setReal(ui.itsAttr[OpalScalingFFAMagnet::TAN_DELTA], 5);
    Attributes::setReal(ui.itsAttr[OpalScalingFFAMagnet::MAX_Y_POWER], 6);
    //Attributes::setString(ui.itsAttr[OpalScalingFFAMagnet::END_FIELD_MODEL], "");
    Attributes::setReal(ui.itsAttr[OpalScalingFFAMagnet::END_LENGTH], 8);
    Attributes::setReal(ui.itsAttr[OpalScalingFFAMagnet::RADIAL_NEG_EXTENT], 10);
    Attributes::setReal(ui.itsAttr[OpalScalingFFAMagnet::RADIAL_POS_EXTENT], 11);
    Attributes::setReal(ui.itsAttr[OpalScalingFFAMagnet::HEIGHT], 12);
    Attributes::setReal(ui.itsAttr[OpalScalingFFAMagnet::LAYOUT_START], 13);
    Attributes::setReal(ui.itsAttr[OpalScalingFFAMagnet::AZIMUTHAL_EXTENT], 14);
    // Update the magnet
    EXPECT_NO_THROW(ui.update());
    // Check the values
    auto* ffa = dynamic_cast<ScalingFFAMagnet*>(ui.getElement());
    ASSERT_TRUE(ffa);
    ffa->setupEndField();
    EXPECT_NEAR(ffa->getEndField()->getCentreLength(), 1.0/ffa->getR0(), 1e-12);
    EXPECT_NEAR(ffa->getDipoleConstant(), 2.0, 1e-12);
    EXPECT_NEAR(ffa->getR0(), 19.0, 1e-12);
    EXPECT_NEAR(ffa->getFieldIndex(), 4, 1e-12);
    EXPECT_NEAR(ffa->getTanDelta(), 5, 1e-12);
    EXPECT_EQ(ffa->getMaxOrder(), 6);
    EXPECT_NEAR(ffa->getEndField()->getEndLength(), 8/ffa->getR0(), 1e-12);
    EXPECT_NEAR(ffa->getRMin(), ffa->getR0()-10, 1e-12);
    EXPECT_NEAR(ffa->getRMax(), ffa->getR0()+11, 1e-12);
    EXPECT_NEAR(ffa->getVerticalExtent()*2.0, 12, 1e-12);
    // phistart is the (MAGNET_START + LENGTH/2)/R0 [radians]
    EXPECT_NEAR(ffa->getPhiStart()*ffa->getR0(), 13+0.5, 1e-12);
    EXPECT_NEAR(ffa->getAzimuthalExtent()*ffa->getR0(), 14, 1e-12);
}

TEST_F(TestOpalScalingFFAMagnet, CentreLength) {
    OpalScalingFFAMagnet ui1;
    Attributes::setReal(ui1.itsAttr[OpalElement::LENGTH], 2.0);
    Attributes::setReal(ui1.itsAttr[OpalScalingFFAMagnet::R0], 4.0);
    EXPECT_NO_THROW(ui1.update());
    auto* ffa1 = dynamic_cast<ScalingFFAMagnet*>(ui1.getElement());
    EXPECT_NEAR(ffa1->getEndField()->getCentreLength(), 2.0/4.0, 1e-12);

    OpalScalingFFAMagnet ui2;
    Attributes::setReal(ui2.itsAttr[OpalScalingFFAMagnet::CENTRE_LENGTH], 3.0);
    Attributes::setReal(ui2.itsAttr[OpalScalingFFAMagnet::R0], 4.0);
    EXPECT_NO_THROW(ui2.update());
    auto* ffa2 = dynamic_cast<ScalingFFAMagnet*>(ui2.getElement());
    EXPECT_NEAR(ffa2->getEndField()->getCentreLength(), 3.0/4.0, 1e-12);

    OpalScalingFFAMagnet ui3;
    Attributes::setReal(ui3.itsAttr[OpalElement::LENGTH], 5.0); // should take this value
    Attributes::setReal(ui3.itsAttr[OpalScalingFFAMagnet::CENTRE_LENGTH], 4.0); // not this value
    Attributes::setReal(ui3.itsAttr[OpalScalingFFAMagnet::R0], 4.0);
    EXPECT_NO_THROW(ui3.update());
    auto* ffa3 = dynamic_cast<ScalingFFAMagnet*>(ui3.getElement());
    EXPECT_NEAR(ffa3->getEndField()->getCentreLength(), 5.0/4.0, 1e-12);
}

TEST_F(TestOpalScalingFFAMagnet, Aperture) {
    OpalScalingFFAMagnet ui;
    Attributes::setString(ui.itsAttr[OpalElement::APERT], "ABC");
    EXPECT_ANY_THROW(ui.update());
}

TEST_F(TestOpalScalingFFAMagnet, EndFieldModel) {
    std::string endName = "TestEndField";
    auto endField = std::make_shared<endfieldmodel::Tanh>(22.0, 3.0, 5);
    endfieldmodel::EndFieldModelManager::getEFMManager()->setEndFieldModel(endName, endField);

    OpalScalingFFAMagnet ui;
    Attributes::setString(ui.itsAttr[OpalScalingFFAMagnet::END_FIELD_MODEL], "TestEndField");
    Attributes::setReal(ui.itsAttr[OpalScalingFFAMagnet::R0], 5.0);
    EXPECT_NO_THROW(ui.update());
    auto* ffa = dynamic_cast<ScalingFFAMagnet*>(ui.getElement());
    ffa->setupEndField();
    EXPECT_NEAR(ffa->getEndField()->getCentreLength()/2, 22.0/ffa->getR0(), 1e-12);
}

TEST_F(TestOpalScalingFFAMagnet, AzimuthalExtentDefault) {
    OpalScalingFFAMagnet ui;
    Attributes::setReal(ui.itsAttr[OpalElement::LENGTH], 1);
    Attributes::setReal(ui.itsAttr[OpalScalingFFAMagnet::R0], 19);
    Attributes::setReal(ui.itsAttr[OpalScalingFFAMagnet::END_LENGTH], 8);
    ui.update();
    auto* ffa = dynamic_cast<ScalingFFAMagnet*>(ui.getElement());
    EXPECT_TRUE(ffa);
    ffa->setupEndField();

    // default azimuthal extent is (endlength*5 + centrelength*0.5)/r0
    EXPECT_NEAR(ffa->getAzimuthalExtent(), (5*8+0.5)/ffa->getR0(), 1e-12);
}
