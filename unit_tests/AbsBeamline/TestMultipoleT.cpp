//
// Tests for MultipoleT
//
// Copyright (c) 2026, Jon Thompson, STFC Rutherford Appleton Laboratory, Didcot, UK
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL.  If not, see <https://www.gnu.org/licenses/>.
//
#include "gtest/gtest.h"
#include "AbsBeamline/MultipoleT.h"

class MultipoleTTest : public testing::Test, public MultipoleT {
public:
    MultipoleTTest() : MultipoleT("Test") {}

protected:
    static void SetUpTestSuite() {
        int argc = 0;
        char** argv = nullptr;

        ippl::initialize(argc, argv);
    }

    static void TearDownTestSuite() {
        ippl::finalize();
    }

    void SetUp() override {
        // nothing special
    }

    void TearDown() override {
        // nothing special
    }

};

TEST_F(MultipoleTTest, TransverseDerivatives) {
    setTransProfile({0.0, 1.0, 2.0, 3.0, 4.0, 5.0});
    EXPECT_DOUBLE_EQ(getTransDeriv(1, 1.0), 55.0);
    EXPECT_DOUBLE_EQ(getTransDeriv(2, 1.0), 170.0);
    EXPECT_DOUBLE_EQ(getTransDeriv(3, 1.0), 414.0);
    EXPECT_DOUBLE_EQ(getTransDeriv(4, 1.0), 696.0);
    EXPECT_DOUBLE_EQ(getTransDeriv(5, 1.0), 600.0);
    EXPECT_DOUBLE_EQ(getTransDeriv(6, 1.0), 0.0);
}

TEST_F(MultipoleTTest, FringeDerivatives) {
    setFringeField(2.0, 1.0, 1.0);
    EXPECT_DOUBLE_EQ(getFringeDeriv(0, -2.0), 0.49966464986953352);
    EXPECT_DOUBLE_EQ(getFringeDeriv(0, 2.0), 0.49966464986953352);
    EXPECT_DOUBLE_EQ(getFringeDeriv(1, -2.0), 0.49932952465848707);
    EXPECT_DOUBLE_EQ(getFringeDeriv(1, 2.0), -0.49932952465848707);
    EXPECT_DOUBLE_EQ(getFringeDeriv(2, -2.0), -0.0013400513070529474);
    EXPECT_DOUBLE_EQ(getFringeDeriv(2, 2.0), -0.0013400513070529474);
    EXPECT_DOUBLE_EQ(getFringeDeriv(3, -2.0), -1.0026765069198489);
    EXPECT_DOUBLE_EQ(getFringeDeriv(3, 2.0), 1.0026765069198489);
    EXPECT_DOUBLE_EQ(getFringeDeriv(4, -2.0), -0.0053386419156264964);
    EXPECT_DOUBLE_EQ(getFringeDeriv(4, 2.0), -0.0053386419156264964);
    EXPECT_DOUBLE_EQ(getFringeDeriv(5, -2.0), 7.9893801387861263);
    EXPECT_DOUBLE_EQ(getFringeDeriv(5, 2.0), -7.9893801387861263);
    EXPECT_DOUBLE_EQ(getFringeDeriv(6, -2.0), -0.021010422121266359);
    EXPECT_DOUBLE_EQ(getFringeDeriv(6, 2.0), -0.021010422121266359);
}
