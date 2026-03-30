#include <gtest/gtest.h>
#include "PartBunch/BunchStateHandler.h"

TEST(BunchStateHandlerTest, DefaultState) {
    BunchStateHandler h;
    EXPECT_FALSE(h.isUnitlessPositions());
    EXPECT_TRUE(h.isMomentsDirty());
    EXPECT_FALSE(h.isEmittingNow());
    EXPECT_TRUE(h.isFirstRepartition());
}

TEST(BunchStateHandlerTest, UnitlessPositionsToggle) {
    BunchStateHandler h;
    h.setUnitlessPositions(true);
    EXPECT_TRUE(h.isUnitlessPositions());
    h.setUnitlessPositions(false);
    EXPECT_FALSE(h.isUnitlessPositions());
}

TEST(BunchStateHandlerTest, MomentsDirtyLifecycle) {
    BunchStateHandler h;

    // Starts dirty
    EXPECT_TRUE(h.isMomentsDirty());

    // Clear it (simulates moments recomputed)
    h.clearMomentsDirty();
    EXPECT_FALSE(h.isMomentsDirty());

    // Mark dirty again (simulates push/kick/emit)
    h.markMomentsDirty();
    EXPECT_TRUE(h.isMomentsDirty());

    // Clear again
    h.clearMomentsDirty();
    EXPECT_FALSE(h.isMomentsDirty());

    // Multiple marks are idempotent
    h.markMomentsDirty();
    h.markMomentsDirty();
    EXPECT_TRUE(h.isMomentsDirty());
}

TEST(BunchStateHandlerTest, EmittingNow) {
    BunchStateHandler h;
    EXPECT_FALSE(h.isEmittingNow());

    h.setEmittingNow(true);
    EXPECT_TRUE(h.isEmittingNow());

    h.setEmittingNow(false);
    EXPECT_FALSE(h.isEmittingNow());
}

TEST(BunchStateHandlerTest, FirstRepartition) {
    BunchStateHandler h;
    EXPECT_TRUE(h.isFirstRepartition());

    h.setFirstRepartition(false);
    EXPECT_FALSE(h.isFirstRepartition());

    h.setFirstRepartition(true);
    EXPECT_TRUE(h.isFirstRepartition());
}
