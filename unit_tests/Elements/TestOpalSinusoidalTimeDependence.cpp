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
#include "Algorithms/SinusoidalTimeDependence.h"
#include "Attributes/Attributes.h"
#include "Elements/OpalSinusoidalTimeDependence.h"
#include "gtest/gtest.h"

TEST(TestOpalSinusoidalTimeDependence, UserInterface) {
    // Make the UI
    OpalSinusoidalTimeDependence ui;
    // Set the attributes
    Attributes::setRealArray(ui.itsAttr[OpalSinusoidalTimeDependence::FREQUENCIES], {1, 2});
    Attributes::setRealArray(ui.itsAttr[OpalSinusoidalTimeDependence::PHASE_OFFSETS], {3, 4});
    Attributes::setRealArray(ui.itsAttr[OpalSinusoidalTimeDependence::AMPLITUDES], {5, 6});
    Attributes::setRealArray(ui.itsAttr[OpalSinusoidalTimeDependence::DC_OFFSETS], {7, 8});
    // Update the object
    EXPECT_NO_THROW(ui.update());
    // Check the values
    auto* myDependency = dynamic_cast<SinusoidalTimeDependence*>(
            AbstractTimeDependence::getTimeDependence("SINUSOIDAL_TIME_DEPENDENCE").get());
    EXPECT_TRUE(myDependency);
    EXPECT_TRUE((myDependency->getFrequencies() == std::vector{1.0, 2.0}));
    EXPECT_TRUE((myDependency->getPhases() == std::vector{3.0, 4.0}));
    EXPECT_TRUE((myDependency->getAmplitudes() == std::vector{5.0, 6.0}));
    EXPECT_TRUE((myDependency->getOffsets() == std::vector{7.0, 8.0}));
}

TEST(TestOpalSinusoidalTimeDependence, ConstructorTest) {
    OpalSinusoidalTimeDependence dep;
    const OpalSinusoidalTimeDependence* dep_clone = dep.clone("new name");
    EXPECT_EQ(dep_clone->getOpalName(), "new name");
}

TEST(TestOpalSinusoidalTimeDependence, PrintTest) {
    const OpalSinusoidalTimeDependence dep;
    std::stringstream _string;
    dep.print(_string);
    EXPECT_EQ(_string.str(), "SINUSOIDAL_TIME_DEPENDENCE;\n");
}
