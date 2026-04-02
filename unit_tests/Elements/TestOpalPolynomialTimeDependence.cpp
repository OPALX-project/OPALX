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

#include <sstream>
#include "Algorithms/AbstractTimeDependence.h"
#include "Algorithms/PolynomialTimeDependence.h"
#include "Attributes/Attributes.h"
#include "Elements/OpalPolynomialTimeDependence.h"
#include "gtest/gtest.h"

TEST(TestOpalPolynomialTimeDependence, ConstructorTest) {
    OpalPolynomialTimeDependence dep;
    const OpalPolynomialTimeDependence* dep_clone = dep.clone("new name");
    EXPECT_EQ(dep_clone->getOpalName(), "new name");
}

TEST(TestOpalPolynomialTimeDependence, PrintTest) {
    const OpalPolynomialTimeDependence dep;
    std::stringstream _string;
    dep.print(_string);
    EXPECT_EQ(_string.str(), "POLYNOMIAL_TIME_DEPENDENCE;\n");
}

TEST(TestOpalPolynomialTimeDependence, UpdateTest) {
    OpalPolynomialTimeDependence dependence;
    dependence.setOpalName("SCALE1");
    Attributes::setReal(dependence.itsAttr[OpalPolynomialTimeDependence::P0], 1.0);
    Attributes::setReal(dependence.itsAttr[OpalPolynomialTimeDependence::P1], 2.0);
    Attributes::setReal(dependence.itsAttr[OpalPolynomialTimeDependence::P2], 3.0);
    Attributes::setReal(dependence.itsAttr[OpalPolynomialTimeDependence::P3], 4.0);
    dependence.update();
    const auto dep          = AbstractTimeDependence::getTimeDependence("SCALE1");
    const auto value        = dep->getValue(0.1);
    constexpr auto shouldBe = 1.0 + 2.0 * 0.1 + 3.0 * 0.1 * 0.1 + 4.0 * 0.1 * 0.1 * 0.1;
    EXPECT_DOUBLE_EQ(value, shouldBe);
}

TEST(TestOpalPolynomialTimeDependence, UpdateTest2) {
    OpalPolynomialTimeDependence dependence;
    dependence.setOpalName("SCALE1");
    Attributes::setRealArray(
            dependence.itsAttr[OpalPolynomialTimeDependence::COEFFICIENTS], {1.0, 2.0, 3.0, 4.0});
    dependence.update();
    const auto dep          = AbstractTimeDependence::getTimeDependence("SCALE1");
    const auto value        = dep->getValue(0.1);
    constexpr auto shouldBe = 1.0 + 2.0 * 0.1 + 3.0 * 0.1 * 0.1 + 4.0 * 0.1 * 0.1 * 0.1;
    EXPECT_DOUBLE_EQ(value, shouldBe);
}

TEST(TestOpalPolynomialTimeDependence, UpdateTest3) {
    OpalPolynomialTimeDependence dependence;
    dependence.setOpalName("SCALE1");
    Attributes::setRealArray(
            dependence.itsAttr[OpalPolynomialTimeDependence::COEFFICIENTS], {1.0, 2.0, 3.0, 4.0});
    Attributes::setReal(dependence.itsAttr[OpalPolynomialTimeDependence::P0], 1.0);
    Attributes::setReal(dependence.itsAttr[OpalPolynomialTimeDependence::P1], 2.0);
    Attributes::setReal(dependence.itsAttr[OpalPolynomialTimeDependence::P2], 3.0);
    Attributes::setReal(dependence.itsAttr[OpalPolynomialTimeDependence::P3], 4.0);
    EXPECT_ANY_THROW(dependence.update());
}

TEST(TestOpalPolynomialTimeDependence, UserInterfacePN) {
    // Make the UI
    OpalPolynomialTimeDependence ui;
    // Set the attributes
    Attributes::setReal(ui.itsAttr[OpalPolynomialTimeDependence::P0], 2);
    Attributes::setReal(ui.itsAttr[OpalPolynomialTimeDependence::P1], 3);
    Attributes::setReal(ui.itsAttr[OpalPolynomialTimeDependence::P2], 4);
    Attributes::setReal(ui.itsAttr[OpalPolynomialTimeDependence::P3], 5);
    // Update the object
    EXPECT_NO_THROW(ui.update());
    // Check the values
    auto* myDependency = dynamic_cast<PolynomialTimeDependence*>(
            AbstractTimeDependence::getTimeDependence("POLYNOMIAL_TIME_DEPENDENCE").get());
    EXPECT_TRUE(myDependency);
    EXPECT_NEAR(myDependency->getCoefficients()[0], 2, 1e-10);
    EXPECT_NEAR(myDependency->getCoefficients()[1], 3, 1e-10);
    EXPECT_NEAR(myDependency->getCoefficients()[2], 4, 1e-10);
    EXPECT_NEAR(myDependency->getCoefficients()[3], 5, 1e-10);
    EXPECT_EQ(myDependency->getCoefficients().size(), 4);
}

TEST(TestOpalPolynomialTimeDependence, UserInterfaceCoeffs) {
    // Make the UI
    OpalPolynomialTimeDependence ui;
    // Set the attributes
    Attributes::setRealArray(
            ui.itsAttr[OpalPolynomialTimeDependence::COEFFICIENTS], {2, 3, 4, 5, 6});
    // Update the object
    EXPECT_NO_THROW(ui.update());
    // Check the values
    auto* myDependency = dynamic_cast<PolynomialTimeDependence*>(
            AbstractTimeDependence::getTimeDependence("POLYNOMIAL_TIME_DEPENDENCE").get());
    EXPECT_TRUE(myDependency);
    EXPECT_NEAR(myDependency->getCoefficients()[0], 2, 1e-10);
    EXPECT_NEAR(myDependency->getCoefficients()[1], 3, 1e-10);
    EXPECT_NEAR(myDependency->getCoefficients()[2], 4, 1e-10);
    EXPECT_NEAR(myDependency->getCoefficients()[3], 5, 1e-10);
    EXPECT_NEAR(myDependency->getCoefficients()[4], 6, 1e-10);
    EXPECT_EQ(myDependency->getCoefficients().size(), 5);
}

TEST(TestOpalPolynomialTimeDependence, UserInterfaceBoth) {
    // Make the UI
    OpalPolynomialTimeDependence ui;
    // Set the attributes
    Attributes::setRealArray(
            ui.itsAttr[OpalPolynomialTimeDependence::COEFFICIENTS], {2, 3, 4, 5, 6});
    Attributes::setReal(ui.itsAttr[OpalPolynomialTimeDependence::P0], 2);
    Attributes::setReal(ui.itsAttr[OpalPolynomialTimeDependence::P1], 3);
    Attributes::setReal(ui.itsAttr[OpalPolynomialTimeDependence::P2], 4);
    Attributes::setReal(ui.itsAttr[OpalPolynomialTimeDependence::P3], 5);
    // Update the object
    EXPECT_THROW(ui.update(), std::invalid_argument);
}
