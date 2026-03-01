/**
 * \file TestConstantEzRep.cpp
 * \brief Unit tests for ConstantEzRep (BeamlineCore layer).
 *
 * Tests: geometry, field, setElementLength/setEz, clone, getChannel("L")/"EZ").
 */

#include <gtest/gtest.h>

#include "BeamlineCore/ConstantEzRep.h"
#include "Channels/IndirectChannel.h"

#include <cmath>
#include <memory>

namespace {

class ConstantEzRepTest : public ::testing::Test {
protected:
    void SetUp() override {
        rep_ = std::make_unique<ConstantEzRep>("TestRep");
        rep_->setElementLength(2.0);
        rep_->setEz(5.0);
    }

    std::unique_ptr<ConstantEzRep> rep_;
};

// ---------------------------------------------------------------------------
// Geometry and field
// ---------------------------------------------------------------------------
TEST_F(ConstantEzRepTest, GeometryLength) {
    EXPECT_DOUBLE_EQ(rep_->getGeometry().getElementLength(), 2.0);
    rep_->setElementLength(3.5);
    EXPECT_DOUBLE_EQ(rep_->getGeometry().getElementLength(), 3.5);
    EXPECT_DOUBLE_EQ(rep_->getElementLength(), 3.5);
}

TEST_F(ConstantEzRepTest, FieldEz) {
    EXPECT_DOUBLE_EQ(rep_->getField().getEz(), 5.0);
    rep_->setEz(-1.5);
    EXPECT_DOUBLE_EQ(rep_->getField().getEz(), -1.5);
    EXPECT_DOUBLE_EQ(rep_->getEz(), -1.5);
}

TEST_F(ConstantEzRepTest, SetEzSyncsBaseAndField) {
    rep_->setEz(100.0);
    EXPECT_DOUBLE_EQ(rep_->getEz(), 100.0);
    EXPECT_DOUBLE_EQ(rep_->getField().getEz(), 100.0);
}

// ---------------------------------------------------------------------------
// Clone
// ---------------------------------------------------------------------------
TEST_F(ConstantEzRepTest, Clone) {
    std::unique_ptr<ElementBase> copy(rep_->clone());
    ASSERT_NE(copy.get(), nullptr);
    ConstantEzRep* repCopy = dynamic_cast<ConstantEzRep*>(copy.get());
    ASSERT_NE(repCopy, nullptr);
    EXPECT_DOUBLE_EQ(repCopy->getElementLength(), 2.0);
    EXPECT_DOUBLE_EQ(repCopy->getEz(), 5.0);
    EXPECT_DOUBLE_EQ(repCopy->getField().getEz(), 5.0);
}

// ---------------------------------------------------------------------------
// Channels L and EZ
// ---------------------------------------------------------------------------
TEST_F(ConstantEzRepTest, ChannelL) {
    Channel* ch = rep_->getChannel("L", false);
    ASSERT_NE(ch, nullptr);
    EXPECT_TRUE(ch->set(1.5));
    EXPECT_DOUBLE_EQ(rep_->getElementLength(), 1.5);
    EXPECT_DOUBLE_EQ(rep_->getGeometry().getElementLength(), 1.5);
    double val = 0.0;
    EXPECT_TRUE(ch->get(val));
    EXPECT_DOUBLE_EQ(val, 1.5);
    delete ch;
}

TEST_F(ConstantEzRepTest, ChannelEZ) {
    Channel* ch = rep_->getChannel("EZ", false);
    ASSERT_NE(ch, nullptr);
    EXPECT_TRUE(ch->set(7.0));
    EXPECT_DOUBLE_EQ(rep_->getEz(), 7.0);
    EXPECT_DOUBLE_EQ(rep_->getField().getEz(), 7.0);
    double val = 0.0;
    EXPECT_TRUE(ch->get(val));
    EXPECT_DOUBLE_EQ(val, 7.0);
    delete ch;
}

TEST_F(ConstantEzRepTest, ChannelUnknownReturnsNull) {
    Channel* ch = rep_->getChannel("UNKNOWN", false);
    EXPECT_EQ(ch, nullptr);
}

}  // namespace
