//
// Unit tests for OpalVariableRFCavity element definition
//
// Copyright (c) 2014, Chris Rogers, STFC Rutherford Appleton Laboratory, Didcot, UK
// All rights reserved.
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#include "AbsBeamline/VariableRFCavity.h"
#include "Algorithms/PolynomialTimeDependence.h"
#include "Attributes/Attributes.h"
#include "Elements/OpalVariableRFCavity.h"
#include "Structure/DataSink.h"
#include "gtest/gtest.h"

class TestOpalVariableRFCavity : public testing::Test {
public:
    TestOpalVariableRFCavity() = default;

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

TEST_F(TestOpalVariableRFCavity, ConstructorDestructor) {
    OpalVariableRFCavity cav1;
    EXPECT_EQ(cav1.getOpalName(), "VARIABLE_RF_CAVITY");
    OpalVariableRFCavity cav2("name", &cav1);
    EXPECT_EQ(cav2.getOpalName(), "name");
    const std::unique_ptr<OpalVariableRFCavity> cav3{cav2.clone()};
    EXPECT_EQ(cav3->getOpalName(), "name");
    const std::unique_ptr<OpalVariableRFCavity> cav4{cav2.clone("other_name")};
    EXPECT_EQ(cav4->getOpalName(), "other_name");
}

TEST_F(TestOpalVariableRFCavity, UserInterface) {
    // Declare some time dependencies
    auto pd1 = std::make_shared<PolynomialTimeDependence>(std::vector{1.0});
    auto pd2 = std::make_shared<PolynomialTimeDependence>(std::vector{2.0});
    auto pd3 = std::make_shared<PolynomialTimeDependence>(std::vector{3.0});
    AbstractTimeDependence::setTimeDependence("PD1", pd1);
    AbstractTimeDependence::setTimeDependence("PD2", pd2);
    AbstractTimeDependence::setTimeDependence("PD3", pd3);
    // Make the UI
    OpalVariableRFCavity ui;
    // Set the attributes
    Attributes::setUpperCaseString(ui.itsAttr[OpalVariableRFCavity::PHASE_MODEL], "pd1");
    Attributes::setUpperCaseString(ui.itsAttr[OpalVariableRFCavity::AMPLITUDE_MODEL], "pd2");
    Attributes::setUpperCaseString(ui.itsAttr[OpalVariableRFCavity::FREQUENCY_MODEL], "pd3");
    Attributes::setReal(ui.itsAttr[OpalVariableRFCavity::WIDTH], 4.0);
    Attributes::setReal(ui.itsAttr[OpalVariableRFCavity::HEIGHT], 5.0);
    Attributes::setReal(ui.itsAttr[OpalVariableRFCavity::LENGTH], 6.0);
    // Update the cavity
    EXPECT_NO_THROW(ui.update());
    // Check the values
    auto* myCavity = dynamic_cast<VariableRFCavity*>(ui.getElement());
    EXPECT_TRUE(myCavity);
    EXPECT_NO_THROW(myCavity->initialiseTimeDependencies());
    EXPECT_DOUBLE_EQ(myCavity->getPhaseModel()->getValue(0.0), 1.0);
    EXPECT_DOUBLE_EQ(myCavity->getAmplitudeModel()->getValue(0.0), 2.0);
    EXPECT_DOUBLE_EQ(myCavity->getFrequencyModel()->getValue(0.0), 3.0);
    EXPECT_DOUBLE_EQ(myCavity->getWidth(), 4.0);
    EXPECT_DOUBLE_EQ(myCavity->getHeight(), 5.0);
    EXPECT_DOUBLE_EQ(myCavity->getLength(), 6.0);
}
