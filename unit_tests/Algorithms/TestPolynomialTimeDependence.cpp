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

#include "Algorithms/PolynomialTimeDependence.h"
#include "Attributes/Attributes.h"
#include "Ippl.h"
#include "Utilities/GeneralClassicException.h"
#include "gtest/gtest.h"

TEST(TestPolynomialTimeDependence, PolynomialTimeDependenceTest) {
    // Check empty polynomial coefficients always returns 0.
    std::vector<double> test;
    PolynomialTimeDependence time_dependence_1(test);
    EXPECT_DOUBLE_EQ(time_dependence_1.getValue(0.1), 0.);

    // Check constant term produces constant
    test.push_back(1.);
    PolynomialTimeDependence time_dependence_2(test);
    EXPECT_DOUBLE_EQ(time_dependence_2.getValue(0.1), 1.);

    // Check cubic terms
    test.push_back(2.);
    test.push_back(3.);
    PolynomialTimeDependence time_dependence_3(test);
    EXPECT_DOUBLE_EQ(time_dependence_3.getValue(0.1), 1.23);

    // Check clone produces same result
    PolynomialTimeDependence* time_dependence_clone = time_dependence_3.clone();
    EXPECT_DOUBLE_EQ(time_dependence_clone->getValue(0.1), 1.23);
    delete time_dependence_clone;
}

TEST(TestPolynomialTimeDependence, TDMapTest) {
    // throw on empty value
    EXPECT_THROW(AbstractTimeDependence::getTimeDependence("name"), GeneralClassicException);
    // set/get time dependence
    PolynomialTimeDependence time_dep(std::vector<double>{});
    const std::shared_ptr<PolynomialTimeDependence> td1(time_dep.clone());
    AbstractTimeDependence::setTimeDependence("td1", td1);
    EXPECT_EQ(AbstractTimeDependence::getTimeDependence("td1"), td1);
    const std::shared_ptr<PolynomialTimeDependence> td2(time_dep.clone());
    AbstractTimeDependence::setTimeDependence("td2", td2);
    EXPECT_EQ(AbstractTimeDependence::getTimeDependence("td2"), td2);
    EXPECT_EQ(AbstractTimeDependence::getTimeDependence("td1"), td1);
    // set time dependence overwriting existing time dependence
    // should overwrite, without memory leak
    const std::shared_ptr<PolynomialTimeDependence> td3(time_dep.clone());
    AbstractTimeDependence::setTimeDependence("td1", td3);
    EXPECT_EQ(AbstractTimeDependence::getTimeDependence("td1"), td3);
}

TEST(TestPolynomialTimeDependence, TDMapNameLookupTest) {
    EXPECT_THROW(AbstractTimeDependence::getName(nullptr), GeneralClassicException);
    PolynomialTimeDependence time_dep(std::vector<double>(1, 1));
    std::shared_ptr<PolynomialTimeDependence> td1(time_dep.clone());
    std::shared_ptr<PolynomialTimeDependence> td2(time_dep.clone());
    std::shared_ptr<PolynomialTimeDependence> td3(time_dep.clone());
    AbstractTimeDependence::setTimeDependence("td1", td1);
    AbstractTimeDependence::setTimeDependence("td2", td2);
    AbstractTimeDependence::setTimeDependence("td3", td2);
    std::string name1 = AbstractTimeDependence::getName(td1);
    EXPECT_EQ(name1, "td1");
    std::string name2 = AbstractTimeDependence::getName(td2);
    EXPECT_TRUE(name2 == "td2" || name2 == "td3");
    EXPECT_THROW(AbstractTimeDependence::getName(td3), GeneralClassicException);
}

TEST(TestPolynomialTimeDependence, Integral) {
    // Check empty polynomial coefficients always returns 0.
    PolynomialTimeDependence time_dependence_1(std::vector<double>{});
    EXPECT_DOUBLE_EQ(time_dependence_1.getIntegral(0.1), 0.0);

    // Check constant term produces constant
    PolynomialTimeDependence time_dependence_2({1.0});
    EXPECT_DOUBLE_EQ(time_dependence_2.getIntegral(0.1), 0.1);

    // Check cubic terms
    PolynomialTimeDependence time_dependence_3({1.0, 2.0, 3.0});
    EXPECT_DOUBLE_EQ(time_dependence_3.getIntegral(0.1), 0.111);
}

TEST(TestPolynomialTimeDependence, Print) {
    int argc    = 0;
    char** argv = nullptr;
    ippl::initialize(argc, argv);
    std::stringstream ss;
    Inform inform("Test", ss);
    inform.setOutputLevel(5);
    PolynomialTimeDependence time_dependence({1.0, 2.0, 3.0});
    inform << time_dependence << endl;
    EXPECT_STREQ(
            "Test> c0= 1.000000e+00 c1= 2.000000e+00 c2= 3.000000e+00 \nTest> \n",
            ss.str().c_str());
    ippl::finalize();
}
