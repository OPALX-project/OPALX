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
#include "Algorithms/SplineTimeDependence.h"
#include "Attributes/Attributes.h"
#include "Elements/OpalSplineTimeDependence.h"
#include "gtest/gtest.h"

TEST(TestOpalSplineTimeDependence, ConstructorTest) {
    OpalSplineTimeDependence dep;
    const OpalSplineTimeDependence* dep_clone = dep.clone("new name");
    EXPECT_EQ(dep_clone->getOpalName(), "new name");
}

TEST(TestOpalSplineTimeDependence, PrintTest) {
    const OpalSplineTimeDependence dep;
    std::stringstream _string;
    dep.print(_string);
    EXPECT_EQ(_string.str(), "SPLINE_TIME_DEPENDENCE;\n");
}

TEST(TestOpalSplineTimeDependence, UserInterface) {
    // Make the UI
    OpalSplineTimeDependence ui;
    // Set the attributes
    Attributes::setReal(ui.itsAttr[OpalSplineTimeDependence::ORDER], 3);
    Attributes::setRealArray(ui.itsAttr[OpalSplineTimeDependence::TIMES], {1, 2, 3, 4});
    Attributes::setRealArray(ui.itsAttr[OpalSplineTimeDependence::VALUES], {5, 6, 7, 8});
    // Update the object
    EXPECT_NO_THROW(ui.update());
    // Check the values
    auto* myDependency = dynamic_cast<SplineTimeDependence*>(
            AbstractTimeDependence::getTimeDependence("SPLINE_TIME_DEPENDENCE").get());
    EXPECT_TRUE(myDependency);
    EXPECT_EQ(myDependency->getSplineOrder(), 3);
    EXPECT_EQ(myDependency->getTimes(), (std::vector{1.0, 2.0, 3.0, 4.0}));
    EXPECT_EQ(myDependency->getValues(), (std::vector{5.0, 6.0, 7.0, 8.0}));
}

TEST(TestOpalSplineTimeDependence, DefaultValues) {
    // Make the UI
    OpalSplineTimeDependence ui;
    // Set the attributes
    Attributes::setReal(ui.itsAttr[OpalSplineTimeDependence::ORDER], 1);
    // Update the object
    EXPECT_NO_THROW(ui.update());
    // Check the values
    auto* myDependency = dynamic_cast<SplineTimeDependence*>(
            AbstractTimeDependence::getTimeDependence("SPLINE_TIME_DEPENDENCE").get());
    EXPECT_TRUE(myDependency);
    EXPECT_EQ(myDependency->getSplineOrder(), 1);
    EXPECT_EQ(myDependency->getTimes(), (std::vector{0.0, 1.0}));
    EXPECT_EQ(myDependency->getValues(), (std::vector{0.0, 0.0}));
}

TEST(TestOpalSplineTimeDependence, UserInterfaceBadOrder) {
    // Make the UI
    OpalSplineTimeDependence ui;
    // Set the attributes
    Attributes::setReal(ui.itsAttr[OpalSplineTimeDependence::ORDER], 2);
    Attributes::setRealArray(ui.itsAttr[OpalSplineTimeDependence::TIMES], {1, 2, 3, 4});
    Attributes::setRealArray(ui.itsAttr[OpalSplineTimeDependence::VALUES], {5, 6, 7, 8});
    // Update the object
    EXPECT_ANY_THROW(ui.update());
}
