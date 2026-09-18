//
// Unit tests for class Enge
//
// Copyright (c) 2017-2026, Chris Rogers, STFC Rutherford Appleton Laboratory, Didcot, UK
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


#include "AbsBeamline/EndFieldModel/AsymmetricEnge.h"
#include "gtest/gtest.h"

class TestAsymmetricEnge : public testing::Test {
public:
    TestAsymmetricEnge() = default;

    static void SetUpTestSuite() {
        Kokkos::initialize();
    }
    static void TearDownTestSuite() {
        Kokkos::finalize();
    }
};

TEST_F(TestAsymmetricEnge, ConstructorTest) {
    //x0, lambda, max_index
    endfieldmodel::Enge enge0;
    EXPECT_EQ(enge0.getCentreLength(), 0.0);
    EXPECT_EQ(enge0.getEndLength(), 0.0);

    std::vector<double> a0 = {5, 6};
    std::vector<double> a1 = {4, 8};
    auto enge1 = std::make_shared<endfieldmodel::AsymmetricEnge>(a0, 2, 7, a1, 3, 9);
    EXPECT_EQ(enge1->getX0Start(), 2.0);
    EXPECT_EQ(enge1->getLambdaStart(), 7.0);
    EXPECT_EQ(enge1->getX0End(), 3.0);
    EXPECT_EQ(enge1->getLambdaEnd(), 9.0);

    endfieldmodel::AsymmetricEnge* enge2 = enge1->clone();
    EXPECT_EQ(enge1->getX0Start(), 2.0);
    EXPECT_EQ(enge1->getLambdaStart(), 7.0);
    EXPECT_EQ(enge1->getX0End(), 3.0);
    EXPECT_EQ(enge1->getLambdaEnd(), 9.0);
    delete enge2;
}


TEST_F(TestAsymmetricEnge, DerivativeTest) {
    std::vector<double> xVector = {-2.5, 0.0, 0.0, 3.0, 0.0};
    double dx = 1e-6;
    double x0S = 2.5;
    double lambdaS = 0.1;
    double x0E = 3.0;
    double lambdaE = 0.15;
    std::vector<double> aS = {0.0, 1.0, 0.0};
    std::vector<double> aE = {0.0, 1.0, 0.0};
    endfieldmodel::AsymmetricEnge enge(aS, x0S, lambdaS, aE, x0E, lambdaE);
    endfieldmodel::Enge engeS(aS, x0S, lambdaS);
    endfieldmodel::Enge engeE(aE, x0E, lambdaE);

    for (double x = -10.0; x < 10.1; x += 0.5) {
        std::cerr << x << "          ";
        std::cerr << enge.function(x, 0) << "             ";
        std::cerr << engeS.function(x, 0) << "             ";
        std::cerr << engeE.function(x, 0) << "             ";
        std::cerr << std::endl;
    }

    EXPECT_NEAR(enge.function(-5.0, 0), 0.0, 1e-2);
    EXPECT_NEAR(enge.function(-2.5, 0), 0.5, 1e-6); // x0S from the centre
    EXPECT_NEAR(enge.function(0.0, 0), 1.0, 1e-6);
    EXPECT_NEAR(enge.function(3.0, 0), 0.5, 1e-6); // x0E from the centre
    EXPECT_NEAR(enge.function(5.0, 0), 0.0, 1e-2);

    for (auto x: xVector) {
        for (size_t n = 1; n < 5; ++n) {
            double yP = enge.function(x+dx, n-1);
            double yM = enge.function(x-dx, n-1);
            double dyTest = enge.function(x, n);
            EXPECT_NEAR(dyTest, (yP-yM)/2/dx, 1e-6) << " at x " << x << " n " << n;
        }
    }

}


