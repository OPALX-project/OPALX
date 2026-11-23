#include <gtest/gtest.h>

#include "Ippl.h"
#include "PartBunch/BunchStateHandler.h"

#include "Utilities/OpalException.h"

#include <array>
#include <limits>

class BunchStateHandlerTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
    }

    static void TearDownTestSuite() { ippl::finalize(); }
};

TEST_F(BunchStateHandlerTest, RegisterContainerReturnsFreshSlot) {
    BunchStateHandler h;
    auto s = h.registerContainer();
    ASSERT_TRUE(s);
    EXPECT_FALSE(s->unitlessPositions);
}

TEST_F(BunchStateHandlerTest, UnitlessPositionsToggle) {
    BunchStateHandler h;
    auto s = h.registerContainer();

    s->setUnitlessPositions(true);
    EXPECT_TRUE(s->unitlessPositions);
    s->setUnitlessPositions(false);
    EXPECT_FALSE(s->unitlessPositions);
}

TEST_F(BunchStateHandlerTest, PerContainerSlotsAreIndependent) {
    BunchStateHandler h;
    auto a = h.registerContainer();
    auto b = h.registerContainer();
    ASSERT_NE(a.get(), b.get());

    // Switching only `a` must not affect `b`.
    a->setUnitlessPositions(true);
    EXPECT_TRUE(a->unitlessPositions);
    EXPECT_FALSE(b->unitlessPositions);

    // And vice-versa.
    a->setUnitlessPositions(false);
    b->setUnitlessPositions(true);
    EXPECT_FALSE(a->unitlessPositions);
    EXPECT_TRUE(b->unitlessPositions);
}

TEST_F(BunchStateHandlerTest, SlotLifetimeFollowsCallerOwnership) {
    BunchStateHandler h;

    // Keep a weak_ptr to observe when the slot dies.
    std::weak_ptr<BunchStateHandler::ContainerState> observer;
    {
        auto s   = h.registerContainer();
        observer = s;
        EXPECT_FALSE(observer.expired());
    }
    // Last strong ref released: slot is destroyed, handler's weak_ptr expires.
    EXPECT_TRUE(observer.expired());
}

TEST_F(BunchStateHandlerTest, FixedCartesianDomainValidatesAndPersistsExactBounds) {
    BunchStateHandler h;
    const std::array<double, 3> lower{-1.0, -2.0, -3.0};
    const std::array<double, 3> upper{1.0, 2.0, 3.0};

    EXPECT_FALSE(h.fixedCartesianDomain().has_value());
    h.setFixedCartesianDomain(lower, upper);
    ASSERT_TRUE(h.fixedCartesianDomain().has_value());
    EXPECT_EQ(h.fixedCartesianDomain()->lower, lower);
    EXPECT_EQ(h.fixedCartesianDomain()->upper, upper);
    EXPECT_NO_THROW(h.setFixedCartesianDomain(lower, upper));

    EXPECT_THROW(
            h.setFixedCartesianDomain(lower, std::array<double, 3>{1.0, 2.0, 4.0}), OpalException);
    h.clearFixedCartesianDomain();
    EXPECT_FALSE(h.fixedCartesianDomain().has_value());
    EXPECT_NO_THROW(h.setFixedCartesianDomain(lower, std::array<double, 3>{1.0, 2.0, 4.0}));
}

TEST_F(BunchStateHandlerTest, FixedCartesianDomainRejectsInvalidBounds) {
    BunchStateHandler h;
    const std::array<double, 3> lower{-1.0, -2.0, -3.0};
    const std::array<double, 3> upper{1.0, 2.0, 3.0};

    auto equal = upper;
    equal[1]   = lower[1];
    EXPECT_THROW(h.setFixedCartesianDomain(lower, equal), OpalException);

    auto reversed = upper;
    reversed[2]   = lower[2] - 1.0;
    EXPECT_THROW(h.setFixedCartesianDomain(lower, reversed), OpalException);

    auto infinite = upper;
    infinite[0]   = std::numeric_limits<double>::infinity();
    EXPECT_THROW(h.setFixedCartesianDomain(lower, infinite), OpalException);

    auto nan = lower;
    nan[0]   = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(h.setFixedCartesianDomain(nan, upper), OpalException);
    EXPECT_FALSE(h.fixedCartesianDomain().has_value());
}

TEST_F(BunchStateHandlerTest, FixedCartesianDomainIsIndependentOfContainerFlags) {
    BunchStateHandler h;
    auto first  = h.registerContainer();
    auto second = h.registerContainer();
    second->setUnitlessPositions(true);

    h.setFixedCartesianDomain({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0});
    EXPECT_FALSE(first->unitlessPositions);
    EXPECT_TRUE(second->unitlessPositions);

    h.clearFixedCartesianDomain();
    EXPECT_FALSE(first->unitlessPositions);
    EXPECT_TRUE(second->unitlessPositions);
}

// -----------------------------------------------------------------------------
// Aggressive-sync tests: require >1 MPI rank to actually exercise divergence.
// With a single rank the allreduce is a no-op and proves nothing, so skip.
// -----------------------------------------------------------------------------

TEST_F(BunchStateHandlerTest, AggressiveSync_UnitlessPositionsConverge) {
    if (ippl::Comm->size() < 2) {
        GTEST_SKIP() << "Requires at least 2 MPI ranks to exercise divergent setter calls.";
    }

    BunchStateHandler h(true);
    auto s = h.registerContainer();
    EXPECT_FALSE(s->unitlessPositions);

    // Only last rank flips to unitless; all ranks should converge to true.
    const bool isLast = (ippl::Comm->rank() == ippl::Comm->size() - 1);
    s->setUnitlessPositions(isLast);
    EXPECT_TRUE(s->unitlessPositions);
}

TEST_F(BunchStateHandlerTest, AggressiveSync_OffDoesNotSynchronize) {
    if (ippl::Comm->size() < 2) {
        GTEST_SKIP() << "Requires at least 2 MPI ranks to distinguish synced from unsynced.";
    }

    BunchStateHandler h(false);
    auto s = h.registerContainer();

    // Only rank 0 switches; with sync off, other ranks keep physical positions.
    const bool isFirst = ippl::Comm->rank() == 0;
    s->setUnitlessPositions(isFirst);
    EXPECT_EQ(s->unitlessPositions, isFirst);
}
