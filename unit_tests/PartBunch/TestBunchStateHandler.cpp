#include <gtest/gtest.h>

#include "Ippl.h"
#include "PartBunch/BunchStateHandler.h"
#include "Utilities/Options.h"

class BunchStateHandlerTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
    }

    static void TearDownTestSuite() { ippl::finalize(); }

    // Each test starts with the sync flag off; tests that need it turn it on.
    void SetUp() override { Options::aggressiveStateSync = false; }
    void TearDown() override { Options::aggressiveStateSync = false; }
};

TEST_F(BunchStateHandlerTest, RegisterContainerReturnsFreshSlot) {
    BunchStateHandler h;
    auto s = h.registerContainer();
    ASSERT_TRUE(s);
    EXPECT_FALSE(s->unitlessPositions);
    EXPECT_TRUE(s->momentsDirty);  // default is dirty so first compute runs
}

TEST_F(BunchStateHandlerTest, BunchWideFirstRepartitionDefault) {
    BunchStateHandler h;
    EXPECT_TRUE(h.isFirstRepartition());
}

TEST_F(BunchStateHandlerTest, UnitlessPositionsToggle) {
    BunchStateHandler h;
    auto s = h.registerContainer();

    s->setUnitlessPositions(true);
    EXPECT_TRUE(s->unitlessPositions);
    s->setUnitlessPositions(false);
    EXPECT_FALSE(s->unitlessPositions);
}

TEST_F(BunchStateHandlerTest, MomentsDirtyLifecycle) {
    BunchStateHandler h;
    auto s = h.registerContainer();

    // Starts dirty
    EXPECT_TRUE(s->momentsDirty);

    s->clearMomentsDirty();
    EXPECT_FALSE(s->momentsDirty);

    s->markMomentsDirty();
    EXPECT_TRUE(s->momentsDirty);

    s->clearMomentsDirty();
    EXPECT_FALSE(s->momentsDirty);

    // Multiple marks are idempotent
    s->markMomentsDirty();
    s->markMomentsDirty();
    EXPECT_TRUE(s->momentsDirty);
}

TEST_F(BunchStateHandlerTest, FirstRepartition) {
    BunchStateHandler h;
    EXPECT_TRUE(h.isFirstRepartition());

    h.setFirstRepartition(false);
    EXPECT_FALSE(h.isFirstRepartition());

    h.setFirstRepartition(true);
    EXPECT_TRUE(h.isFirstRepartition());
}

TEST_F(BunchStateHandlerTest, PerContainerSlotsAreIndependent) {
    BunchStateHandler h;
    auto a = h.registerContainer();
    auto b = h.registerContainer();
    ASSERT_NE(a.get(), b.get());

    a->clearMomentsDirty();
    b->clearMomentsDirty();

    // Marking only `a` dirty must not affect `b`.
    a->markMomentsDirty();
    EXPECT_TRUE(a->momentsDirty);
    EXPECT_FALSE(b->momentsDirty);

    // And vice-versa.
    a->clearMomentsDirty();
    b->markMomentsDirty();
    EXPECT_FALSE(a->momentsDirty);
    EXPECT_TRUE(b->momentsDirty);
}

TEST_F(BunchStateHandlerTest, SlotLifetimeFollowsCallerOwnership) {
    BunchStateHandler h;

    // Keep a weak_ptr to observe when the slot dies.
    std::weak_ptr<BunchStateHandler::ContainerState> observer;
    {
        auto s = h.registerContainer();
        observer = s;
        EXPECT_FALSE(observer.expired());
    }
    // Last strong ref released: slot is destroyed, handler's weak_ptr expires.
    EXPECT_TRUE(observer.expired());
}

// -----------------------------------------------------------------------------
// Aggressive-sync tests: require >1 MPI rank to actually exercise divergence.
// With a single rank the allreduce is a no-op and proves nothing, so skip.
// -----------------------------------------------------------------------------

TEST_F(BunchStateHandlerTest, AggressiveSync_MomentsDirty_OrAcrossRanks) {
    if (ippl::Comm->size() < 2) {
        GTEST_SKIP() << "Requires at least 2 MPI ranks to exercise divergent setter calls.";
    }

    Options::aggressiveStateSync = true;
    BunchStateHandler h;
    auto s = h.registerContainer();
    s->clearMomentsDirty();  // sync'd clear: all ranks agree false -> false
    ASSERT_FALSE(s->momentsDirty);

    // Only rank 0 marks dirty; logical_or should propagate to every rank.
    if (ippl::Comm->rank() == 0) {
        s->markMomentsDirty();
    } else {
        // Non-zero ranks must also call a setter so they participate in the
        // collective. Clearing locally; the OR with rank 0's "true" wins.
        s->clearMomentsDirty();
    }
    EXPECT_TRUE(s->momentsDirty);
}

TEST_F(BunchStateHandlerTest, AggressiveSync_ClearOnlyWhenAllAgree) {
    if (ippl::Comm->size() < 2) {
        GTEST_SKIP() << "Requires at least 2 MPI ranks to exercise divergent setter calls.";
    }

    Options::aggressiveStateSync = true;
    BunchStateHandler h;
    auto s = h.registerContainer();
    s->markMomentsDirty();
    ASSERT_TRUE(s->momentsDirty);

    // rank 0 -> clearMomentsDirty passes `false`
    // rank i -> markMomentsDirty  passes `true`
    // => OR == true => every rank stays dirty, including rank 0.
    if (ippl::Comm->rank() == 0) {
        s->clearMomentsDirty();
    } else {
        s->markMomentsDirty();
    }
    EXPECT_TRUE(s->momentsDirty);
}

TEST_F(BunchStateHandlerTest, AggressiveSync_UnitlessPositionsConverge) {
    if (ippl::Comm->size() < 2) {
        GTEST_SKIP() << "Requires at least 2 MPI ranks to exercise divergent setter calls.";
    }

    Options::aggressiveStateSync = true;
    BunchStateHandler h;
    auto s = h.registerContainer();
    EXPECT_FALSE(s->unitlessPositions);

    // Only last rank flips to unitless; all ranks should converge to true.
    const bool isLast = (ippl::Comm->rank() == ippl::Comm->size() - 1);
    s->setUnitlessPositions(isLast);
    EXPECT_TRUE(s->unitlessPositions);
}

TEST_F(BunchStateHandlerTest, AggressiveSync_FirstRepartitionConverge) {
    if (ippl::Comm->size() < 2) {
        GTEST_SKIP() << "Requires at least 2 MPI ranks to exercise divergent setter calls.";
    }

    Options::aggressiveStateSync = true;
    BunchStateHandler h;
    ASSERT_TRUE(h.isFirstRepartition());

    // Every rank except rank 0 tries to mark first-repartition done; rank 0
    // still reports true. OR wins => first-repartition remains true everywhere
    // (the conservative choice: do not skip first-time init until all ranks
    // agree it is done).
    h.setFirstRepartition(ippl::Comm->rank() != 0 ? false : true);
    EXPECT_TRUE(h.isFirstRepartition());
}

TEST_F(BunchStateHandlerTest, AggressiveSync_OffDoesNotSynchronize) {
    if (ippl::Comm->size() < 2) {
        GTEST_SKIP() << "Requires at least 2 MPI ranks to distinguish synced from unsynced.";
    }

    Options::aggressiveStateSync = false;
    BunchStateHandler h;
    auto s = h.registerContainer();
    s->clearMomentsDirty();
    ASSERT_FALSE(s->momentsDirty);

    // Only rank 0 marks dirty; with sync off, other ranks must stay clean.
    if (ippl::Comm->rank() == 0) {
        s->markMomentsDirty();
        EXPECT_TRUE(s->momentsDirty);
    } else {
        EXPECT_FALSE(s->momentsDirty);
    }
}
