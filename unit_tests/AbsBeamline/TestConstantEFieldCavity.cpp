/**
 * \file TestConstantEFieldCavity.cpp
 * \brief Unit tests for the ConstantEFieldCavity element.
 *
 * Covers: getType, the geometry length, getEx/Ey/Ez and setters, getFieldExtent,
 * applyOnHost(R, t, E, B) for various z-positions and full-vector application.
 */

#include <gtest/gtest.h>

#include "AbsBeamline/ConstantEFieldCavity.h"
#include "AbsBeamline/ElementBase.h"

#include <cmath>
#include <memory>

namespace {

    class ConstantEFieldCavityTest : public ::testing::Test {
    protected:
        void SetUp() override {
            cavity_ = std::make_unique<ConstantEFieldCavity>("TestConstantEFieldCavity");
            cavity_->getGeometry().setElementLength(1.0);
            cavity_->setEz(10.0);
        }

        std::unique_ptr<ConstantEFieldCavity> cavity_;
    };

    // ---------------------------------------------------------------------------
    // Type and geometry
    // ---------------------------------------------------------------------------
    TEST_F(ConstantEFieldCavityTest, GetType) {
        EXPECT_EQ(cavity_->getType(), ElementType::CONSTANTEFIELDCAVITY);
    }

    TEST_F(ConstantEFieldCavityTest, GeometryLength) {
        EXPECT_DOUBLE_EQ(cavity_->getGeometry().getElementLength(), 1.0);
        cavity_->getGeometry().setElementLength(3.5);
        EXPECT_DOUBLE_EQ(cavity_->getGeometry().getElementLength(), 3.5);
    }

    TEST_F(ConstantEFieldCavityTest, FieldEz) {
        EXPECT_DOUBLE_EQ(cavity_->getEz(), 10.0);
        cavity_->setEz(-1.5);
        EXPECT_DOUBLE_EQ(cavity_->getEz(), -1.5);
    }

    TEST_F(ConstantEFieldCavityTest, GetExEyEzSetters) {
        EXPECT_DOUBLE_EQ(cavity_->getEx(), 0.0);
        EXPECT_DOUBLE_EQ(cavity_->getEy(), 0.0);
        EXPECT_DOUBLE_EQ(cavity_->getEz(), 10.0);

        cavity_->setEx(1.0);
        cavity_->setEy(-2.0);
        cavity_->setEz(3.0);

        EXPECT_DOUBLE_EQ(cavity_->getEx(), 1.0);
        EXPECT_DOUBLE_EQ(cavity_->getEy(), -2.0);
        EXPECT_DOUBLE_EQ(cavity_->getEz(), 3.0);
    }

    TEST_F(ConstantEFieldCavityTest, GetDimensions) {
        cavity_->initialise(nullptr);

        // getFieldExtent() is the local-chart field-support interval [0, L] (same frame as
        // isInside's r), independent of where the element sits along the reference path; the
        // s-position is carried by the placement transform, not by the field extent.
        double zBegin = 0.0, zEnd = 0.0;
        cavity_->getFieldExtent(zBegin, zEnd);
        EXPECT_DOUBLE_EQ(zBegin, 0.0);
        EXPECT_DOUBLE_EQ(zEnd, 1.0);
    }

    // ---------------------------------------------------------------------------
    // applyOnHost(R, t, E, B): position-based, no bunch needed
    // ---------------------------------------------------------------------------
    TEST_F(ConstantEFieldCavityTest, ApplyInside) {
        Vector_t<double, 3> R = {0.0, 0.0, 0.5};
        Vector_t<double, 3> P = {0.0, 0.0, 1.0};
        Vector_t<double, 3> E = {1.0, 2.0, 3.0};
        Vector_t<double, 3> B = {0.0, 0.0, 0.0};

        cavity_->applyOnHost(R, 0.0, E, B);

        EXPECT_DOUBLE_EQ(E(0), 1.0);
        EXPECT_DOUBLE_EQ(E(1), 2.0);
        EXPECT_DOUBLE_EQ(E(2), 3.0 + 10.0);
    }

    TEST_F(ConstantEFieldCavityTest, ApplyBeforeElement) {
        Vector_t<double, 3> R = {0.0, 0.0, -0.1};
        Vector_t<double, 3> P = {0.0, 0.0, 1.0};
        Vector_t<double, 3> E = {1.0, 2.0, 3.0};
        Vector_t<double, 3> B = {0.0, 0.0, 0.0};

        cavity_->applyOnHost(R, 0.0, E, B);

        EXPECT_DOUBLE_EQ(E(0), 1.0);
        EXPECT_DOUBLE_EQ(E(1), 2.0);
        EXPECT_DOUBLE_EQ(E(2), 3.0);
    }

    TEST_F(ConstantEFieldCavityTest, ApplyAfterElement) {
        Vector_t<double, 3> R = {0.0, 0.0, 1.5};
        Vector_t<double, 3> P = {0.0, 0.0, 1.0};
        Vector_t<double, 3> E = {1.0, 2.0, 3.0};
        Vector_t<double, 3> B = {0.0, 0.0, 0.0};

        cavity_->applyOnHost(R, 0.0, E, B);

        EXPECT_DOUBLE_EQ(E(0), 1.0);
        EXPECT_DOUBLE_EQ(E(1), 2.0);
        EXPECT_DOUBLE_EQ(E(2), 3.0);
    }

    TEST_F(ConstantEFieldCavityTest, ApplyAtZ0) {
        Vector_t<double, 3> R = {0.0, 0.0, 0.0};
        Vector_t<double, 3> P = {0.0, 0.0, 1.0};
        Vector_t<double, 3> E = {0.0, 0.0, 0.0};
        Vector_t<double, 3> B = {0.0, 0.0, 0.0};

        cavity_->applyOnHost(R, 0.0, E, B);

        EXPECT_DOUBLE_EQ(E(0), 0.0);
        EXPECT_DOUBLE_EQ(E(1), 0.0);
        EXPECT_DOUBLE_EQ(E(2), 10.0);
    }

    TEST_F(ConstantEFieldCavityTest, ApplyAtZLength) {
        // The field window is [0, L): the exit belongs to the next element.
        Vector_t<double, 3> R = {0.0, 0.0, 1.0};
        Vector_t<double, 3> P = {0.0, 0.0, 1.0};
        Vector_t<double, 3> E = {0.0, 0.0, 0.0};
        Vector_t<double, 3> B = {0.0, 0.0, 0.0};

        cavity_->applyOnHost(R, 0.0, E, B);

        EXPECT_DOUBLE_EQ(E(0), 0.0);
        EXPECT_DOUBLE_EQ(E(1), 0.0);
        EXPECT_DOUBLE_EQ(E(2), 0.0);
    }

    TEST_F(ConstantEFieldCavityTest, ApplyAddsFullVector) {
        cavity_->setEx(1.0);
        cavity_->setEy(-2.0);
        cavity_->setEz(3.0);

        Vector_t<double, 3> R = {0.0, 0.0, 0.5};
        Vector_t<double, 3> P = {0.0, 0.0, 1.0};
        Vector_t<double, 3> E = {0.5, 0.5, 0.5};
        Vector_t<double, 3> B = {0.0, 0.0, 0.0};

        cavity_->applyOnHost(R, 0.0, E, B);

        EXPECT_DOUBLE_EQ(E(0), 0.5 + 1.0);
        EXPECT_DOUBLE_EQ(E(1), 0.5 - 2.0);
        EXPECT_DOUBLE_EQ(E(2), 0.5 + 3.0);
    }

}  // namespace
