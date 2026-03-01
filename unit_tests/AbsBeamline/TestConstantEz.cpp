/**
 * \file TestConstantEz.cpp
 * \brief Unit tests for ConstantEz component (base layer).
 *
 * Tests the ConstantEz API and apply logic via ConstantEzRep (concrete type).
 * Covers: getType, bends, getEz/setEz, getDimensions, apply(R,P,t,E,B),
 * applyToReferenceParticle (inside/outside/boundary).
 */

#include <gtest/gtest.h>

#include "AbsBeamline/ConstantEz.h"
#include "AbsBeamline/ElementBase.h"
#include "BeamlineCore/ConstantEzRep.h"

#include <cmath>
#include <memory>

namespace {

class ConstantEzTest : public ::testing::Test {
protected:
    void SetUp() override {
        rep_ = std::make_unique<ConstantEzRep>("TestConstantEz");
        rep_->setElementLength(1.0);
        rep_->setEz(10.0);
    }

    std::unique_ptr<ConstantEzRep> rep_;
};

// ---------------------------------------------------------------------------
// Type and geometry
// ---------------------------------------------------------------------------
TEST_F(ConstantEzTest, GetType) {
    EXPECT_EQ(rep_->getType(), ElementType::CONSTANTEZ);
}

TEST_F(ConstantEzTest, Bends) {
    EXPECT_FALSE(rep_->bends());
}

TEST_F(ConstantEzTest, GetEzSetEz) {
    EXPECT_DOUBLE_EQ(rep_->getEz(), 10.0);
    rep_->setEz(-5.0);
    EXPECT_DOUBLE_EQ(rep_->getEz(), -5.0);
    rep_->setEz(0.0);
    EXPECT_DOUBLE_EQ(rep_->getEz(), 0.0);
}

TEST_F(ConstantEzTest, GetDimensions) {
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
TEST_F(ConstantEzTest, ApplyInside) {
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

TEST_F(ConstantEzTest, ApplyBeforeElement) {
    Vector_t<double, 3> R = {0.0, 0.0, -0.1};
    Vector_t<double, 3> P = {0.0, 0.0, 1.0};
    Vector_t<double, 3> E = {1.0, 2.0, 3.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    bool out = rep_->apply(R, P, 0.0, E, B);

    EXPECT_FALSE(out);
    EXPECT_DOUBLE_EQ(E(2), 3.0);
}

TEST_F(ConstantEzTest, ApplyAfterElement) {
    Vector_t<double, 3> R = {0.0, 0.0, 1.5};
    Vector_t<double, 3> P = {0.0, 0.0, 1.0};
    Vector_t<double, 3> E = {1.0, 2.0, 3.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    bool out = rep_->apply(R, P, 0.0, E, B);

    EXPECT_FALSE(out);
    EXPECT_DOUBLE_EQ(E(2), 3.0);
}

TEST_F(ConstantEzTest, ApplyAtZ0) {
    Vector_t<double, 3> R = {0.0, 0.0, 0.0};
    Vector_t<double, 3> P = {0.0, 0.0, 1.0};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    bool out = rep_->apply(R, P, 0.0, E, B);

    EXPECT_FALSE(out);
    EXPECT_DOUBLE_EQ(E(2), 10.0);
}

TEST_F(ConstantEzTest, ApplyAtZLength) {
    Vector_t<double, 3> R = {0.0, 0.0, 1.0};
    Vector_t<double, 3> P = {0.0, 0.0, 1.0};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    bool out = rep_->apply(R, P, 0.0, E, B);

    EXPECT_FALSE(out);
    EXPECT_DOUBLE_EQ(E(2), 10.0);
}

// ---------------------------------------------------------------------------
// applyToReferenceParticle
// ---------------------------------------------------------------------------
TEST_F(ConstantEzTest, ApplyToReferenceParticleInside) {
    Vector_t<double, 3> R = {0.0, 0.0, 0.5};
    Vector_t<double, 3> P = {0.0, 0.0, 1.0};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    bool out = rep_->applyToReferenceParticle(R, P, 0.0, E, B);

    EXPECT_FALSE(out);
    EXPECT_DOUBLE_EQ(E(2), 10.0);
}

TEST_F(ConstantEzTest, ApplyToReferenceParticleOutside) {
    Vector_t<double, 3> R = {0.0, 0.0, -0.1};
    Vector_t<double, 3> P = {0.0, 0.0, 1.0};
    Vector_t<double, 3> E = {0.0, 0.0, 0.0};
    Vector_t<double, 3> B = {0.0, 0.0, 0.0};

    bool out = rep_->applyToReferenceParticle(R, P, 0.0, E, B);

    EXPECT_FALSE(out);
    EXPECT_DOUBLE_EQ(E(2), 0.0);
}

}  // namespace
