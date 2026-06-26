/**
 * \file TestConstEzField.cpp
 * \brief Unit tests for ConstEzField (homogeneous E_z field).
 *
 * Tests: Efield(P), Efield(P,t), setEz, getEz, scale.
 */

#include <gtest/gtest.h>

#include "Fields/ConstEzField.h"
#include "Fields/EMField.h"

#include <cmath>

namespace {

    class ConstEzFieldTest : public ::testing::Test {
    protected:
        void SetUp() override { field_.setEz(10.0); }

        ConstEzField field_;
    };

    TEST_F(ConstEzFieldTest, GetEzSetEz) {
        EXPECT_DOUBLE_EQ(field_.getEz(), 10.0);
        field_.setEz(-3.0);
        EXPECT_DOUBLE_EQ(field_.getEz(), -3.0);
        field_.setEz(0.0);
        EXPECT_DOUBLE_EQ(field_.getEz(), 0.0);
    }

    TEST_F(ConstEzFieldTest, EfieldP) {
        Point3D P(1.0, 2.0, 3.0);
        EVector E = field_.Efield(P);
        EXPECT_DOUBLE_EQ(E.getEx(), 0.0);
        EXPECT_DOUBLE_EQ(E.getEy(), 0.0);
        EXPECT_DOUBLE_EQ(E.getEz(), 10.0);
    }

    TEST_F(ConstEzFieldTest, EfieldPT) {
        Point3D P(0.0, 0.0, 0.0);
        EVector E = field_.Efield(P, 1.0);
        EXPECT_DOUBLE_EQ(E.getEx(), 0.0);
        EXPECT_DOUBLE_EQ(E.getEy(), 0.0);
        EXPECT_DOUBLE_EQ(E.getEz(), 10.0);
    }

    TEST_F(ConstEzFieldTest, Scale) {
        field_.scale(2.0);
        EXPECT_DOUBLE_EQ(field_.getEz(), 20.0);
        Point3D P(0, 0, 0);
        EVector E = field_.Efield(P);
        EXPECT_DOUBLE_EQ(E.getEz(), 20.0);

        field_.scale(-0.5);
        EXPECT_DOUBLE_EQ(field_.getEz(), -10.0);
        E = field_.Efield(P);
        EXPECT_DOUBLE_EQ(E.getEz(), -10.0);
    }

}  // namespace
