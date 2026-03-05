/**
 * \file TestConstantEFieldCavity.cpp
 * \brief Unit tests for ConstantEFieldCavity component (base layer).
 *
 * Tests the ConstantEFieldCavity API and apply logic via ConstantEFieldCavityRep (concrete type).
 * Covers: getType, bends, getEx/Ey/Ez and setters, getDimensions,
 * apply(R,P,t,E,B) for various z-positions, full-vector application,
 * and applyToReferenceParticle (inside/outside).
 */

#include <gtest/gtest.h>

#include "AbsBeamline/ElementBase.h"
#include "BeamlineCore/ConstantEFieldCavityRep.h"

#include <cmath>
#include <memory>

namespace {

class ConstantEFieldCavityTest : public ::testing::Test {
protected:
    void SetUp() override {
        rep_ = std::make_unique<ConstantEFieldCavityRep>("TestConstantEFieldCavity");
        rep_->setElementLength(1.0);
        rep_->setEz(10.0);
    }

    std::unique_ptr<ConstantEFieldCavityRep> rep_;
};

// ---------------------------------------------------------------------------
// Type and geometry
// ---------------------------------------------------------------------------
TEST_F(ConstantEFieldCavityTest, GetType) {
    EXPECT_EQ(rep_->getType(), ElementType::CONSTANTEFIELDCAVITY);
}

TEST_F(ConstantEFieldCavityTest, Bends) {
    EXPECT_FALSE(rep_->bends());
}

TEST_F(ConstantEFieldCavityTest, GetExEyEzSetters) {
    EXPECT_DOUBLE_EQ(rep_->getEx(), 0.0);
    EXPECT_DOUBLE_EQ(rep_->getEy(), 0.0);
    EXPECT_DOUBLE_EQ(rep_->getEz(), 10.0);

    rep_->setEx(1.0);
    rep_->setEy(-2.0);
    rep_->setEz(3.0);

    EXPECT_DOUBLE_EQ(rep_->getEx(), 1.0);
    EXPECT_DOUBLE_EQ(rep_->getEy(), -2.0);
    EXPECT_DOUBLE_EQ(rep_->getEz(), 3.0);
}

TEST_F(ConstantEFieldCavityTest, GetDimensions) {
    double startField = 0.0;
    double endField   = 0.0;
    rep_->initialise(nullptr, startField, endField);
    EXPECT_DOUBLE_EQ(endField, 1.0);

    double zBegin = 0.0, zEnd = 0.0;
    rep_->getDimensions(zBegin, zEnd);
    EXPECT_DOUBLE_EQ(zBegin, 0.0);
    EXPECT_DOUBLE_EQ(zEnd, 1.0);

    startField = 2.5;
    endField   = 0.0;
    rep_->initialise(nullptr, startField, endField);
    rep_->getDimensions(zBegin, zEnd);
    EXPECT_DOUBLE_EQ(zBegin, 2.5);
    EXPECT_DOUBLE_EQ(zEnd, 3.5);
}

// ---------------------------------------------------------------------------
// apply(R, P, t, E, B): position-based, no bunch needed
// ---------------------------------------------------------------------------
TEST_F(ConstantEFieldCavityTest, ApplyInside) {
    Vector_t<double, 3> R = {0.0, 0.0, 0.5};
    Vector_t<double, 3> P = {0.0, 0.0, 1.0};
    Vector_t<double, 3> E = {1.0, 2.0, 3.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    bool out = rep_->apply(R, P, 0.0, E, B);

    EXPECT_FALSE(out);
    EXPECT_DOUBLE_EQ(E(0), 1.0);
    EXPECT_DOUBLE_EQ(E(1), 2.0);
    EXPECT_DOUBLE_EQ(E(2), 3.0 + 10.0);
}

TEST_F(ConstantEFieldCavityTest, ApplyBeforeElement) {
    Vector_t<double, 3> R = {0.0, 0.0, -0.1};
    Vector_t<double, 3> P = {0.0, 0.0, 1.0};
    Vector_t<double, 3> E = {1.0, 2.0, 3.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    bool out = rep_->apply(R, P, 0.0, E, B);

    EXPECT_FALSE(out);
    EXPECT_DOUBLE_EQ(E(0), 1.0);
    EXPECT_DOUBLE_EQ(E(1), 2.0);
    EXPECT_DOUBLE_EQ(E(2), 3.0);
}

TEST_F(ConstantEFieldCavityTest, ApplyAfterElement) {
    Vector_t<double, 3> R = {0.0, 0.0, 1.5};
    Vector_t<double, 3> P = {0.0, 0.0, 1.0};
    Vector_t<double, 3> E = {1.0, 2.0, 3.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    bool out = rep_->apply(R, P, 0.0, E, B);

    EXPECT_FALSE(out);
    EXPECT_DOUBLE_EQ(E(0), 1.0);
    EXPECT_DOUBLE_EQ(E(1), 2.0);
    EXPECT_DOUBLE_EQ(E(2), 3.0);
}

TEST_F(ConstantEFieldCavityTest, ApplyAtZ0) {
    Vector_t<double, 3> R = {0.0, 0.0, 0.0};
    Vector_t<double, 3> P = {0.0, 0.0, 1.0};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    bool out = rep_->apply(R, P, 0.0, E, B);

    EXPECT_FALSE(out);
    EXPECT_DOUBLE_EQ(E(0), 0.0);
    EXPECT_DOUBLE_EQ(E(1), 0.0);
    EXPECT_DOUBLE_EQ(E(2), 10.0);
}

TEST_F(ConstantEFieldCavityTest, ApplyAtZLength) {
    Vector_t<double, 3> R = {0.0, 0.0, 1.0};
    Vector_t<double, 3> P = {0.0, 0.0, 1.0};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    bool out = rep_->apply(R, P, 0.0, E, B);

    EXPECT_FALSE(out);
    EXPECT_DOUBLE_EQ(E(0), 0.0);
    EXPECT_DOUBLE_EQ(E(1), 0.0);
    EXPECT_DOUBLE_EQ(E(2), 10.0);
}

TEST_F(ConstantEFieldCavityTest, ApplyAddsFullVector) {
    rep_->setEx(1.0);
    rep_->setEy(-2.0);
    rep_->setEz(3.0);

    Vector_t<double, 3> R = {0.0, 0.0, 0.5};
    Vector_t<double, 3> P = {0.0, 0.0, 1.0};
    Vector_t<double, 3> E = {0.5, 0.5, 0.5};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    bool out = rep_->apply(R, P, 0.0, E, B);

    EXPECT_FALSE(out);
    EXPECT_DOUBLE_EQ(E(0), 0.5 + 1.0);
    EXPECT_DOUBLE_EQ(E(1), 0.5 - 2.0);
    EXPECT_DOUBLE_EQ(E(2), 0.5 + 3.0);
}

// ---------------------------------------------------------------------------
// applyToReferenceParticle
// ---------------------------------------------------------------------------
TEST_F(ConstantEFieldCavityTest, ApplyToReferenceParticleInside) {
    Vector_t<double, 3> R = {0.0, 0.0, 0.5};
    Vector_t<double, 3> P = {0.0, 0.0, 1.0};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    bool out = rep_->applyToReferenceParticle(R, P, 0.0, E, B);

    EXPECT_FALSE(out);
    EXPECT_DOUBLE_EQ(E(0), 0.0);
    EXPECT_DOUBLE_EQ(E(1), 0.0);
    EXPECT_DOUBLE_EQ(E(2), 10.0);
}

TEST_F(ConstantEFieldCavityTest, ApplyToReferenceParticleOutside) {
    Vector_t<double, 3> R = {0.0, 0.0, -0.1};
    Vector_t<double, 3> P = {0.0, 0.0, 1.0};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    bool out = rep_->applyToReferenceParticle(R, P, 0.0, E, B);

    EXPECT_FALSE(out);
    EXPECT_DOUBLE_EQ(E(0), 0.0);
    EXPECT_DOUBLE_EQ(E(1), 0.0);
    EXPECT_DOUBLE_EQ(E(2), 0.0);
}

}  // namespace

