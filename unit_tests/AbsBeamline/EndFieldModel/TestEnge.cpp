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


#include "AbsBeamline/EndFieldModel/Enge.h"
#include "AbstractObjects/OpalData.h"
#include "Physics/Physics.h"
#include "gtest/gtest.h"

class TestEnge : public testing::Test {
public:
    TestEnge() {
    }

    static void SetUpTestSuite() {
        Kokkos::initialize();
    }
    static void TearDownTestSuite() {
        Kokkos::finalize();
    }
};

TEST_F(TestEnge, ConstructorTest) {
    //x0, lambda, max_index
    endfieldmodel::Enge enge0;
    EXPECT_EQ(enge0.getX0(), 0.0);
    EXPECT_EQ(enge0.getLambda(), 0.0);
    EXPECT_EQ(enge0.getCoefficients().size(), 0);

    std::vector<double> a = {5, 6};
    auto enge1 = std::make_shared<endfieldmodel::Enge>(a, 2, 7);
    EXPECT_EQ(enge1->getX0(), 2.0);
    EXPECT_EQ(enge1->getCentreLength(), 4.0);
    EXPECT_EQ(enge1->getEndLength(), 7.0);
    EXPECT_EQ(enge1->getCoefficients(), a);

    endfieldmodel::Enge* enge2 = enge1->clone();
    EXPECT_EQ(enge2->getX0(), enge1->getX0());
    EXPECT_EQ(enge2->getCentreLength(), enge1->getCentreLength());
    EXPECT_EQ(enge2->getCoefficients(), enge2->getCoefficients());
    delete enge2;
}

double myEnge(double x, std::vector<double> a, double x0, double lambda) {
    double deltaX = (x-x0)/lambda;
    double xPow = 1.0;
    double p = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        p += a[i]*xPow;
        xPow *= deltaX;
    }
    double enge = 1/(1.0+exp(p));
    return enge;
}

TEST_F(TestEnge, FunctionTest) {
    std::vector<double> zVector = {0.0, 1.0, 2.0, 3.0};
    endfieldmodel::Enge enge({0.0, 1.0}, 0.0, 0.5);
    for (auto z: zVector) {
        EXPECT_NEAR(enge.getEnge(z, 0), myEnge(z, {0.0, 1.0}, 00.0, 0.5), 1e-12);
    }
}

TEST_F(TestEnge, HNTest) {
    std::vector<double> a = {1.0, 2.0, 3.0, 4.0};
    endfieldmodel::Enge enge = endfieldmodel::Enge(a, 10.0, 0.5);
    enge.setMaximumDerivative(11);
    double dx = 1e-6;
    for(size_t i = 0; i < 10; ++i) {
        double dhdxNumerical = (enge.hN(10.0+dx, i)-
                                enge.hN(10.0-dx, i))/2/dx;
        double dhdx = enge.hN(10.0, i+1);
        EXPECT_NEAR(dhdx, dhdxNumerical, 1e-5)
                << " for " << i << "^th derivative";
    }
}

TEST_F(TestEnge, GNTest) {
    std::vector<double> a = {1.0, 2.0, 3.0, 4.0};
    endfieldmodel::Enge enge = endfieldmodel::Enge(a, 1.0, 0.5);
    enge.setMaximumDerivative(11);
    double dx = 1e-6;
    for(size_t i = 0; i < 10; ++i) {
        double dgdxNumerical = (enge.gN(0.1+dx, i)-
                                enge.gN(0.1-dx, i))/2/dx;
        double dgdx = enge.gN(0.1, i+1);
        EXPECT_NEAR(dgdx/dgdxNumerical, 1.0, 1e-5)
                << " for " << i << "^th derivative";
    }
}
