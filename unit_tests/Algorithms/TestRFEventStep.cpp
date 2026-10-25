// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#include "Algorithms/TrackReferenceStep.h"
#include "Ippl.h"
#include "gtest/gtest.h"

namespace {
    using Vector               = rf_event::Vector;
    constexpr double startTime = 2e-9, nominalStep = 1e-9, mass = 938272088.16;

    // Prescribed SI worldline isolates event scheduling from magnetic physics.
    // Noncommuting impulses make incorrect chronological ordering observable.
    struct DevicePolicy {
        size_t number     = 2;
        double plane[2]   = {.003, .001};
        double minimum[2] = {-1, -1}, maximum[2] = {1, 1};
        double arrival[2]       = {0, 0};
        unsigned long long mask = 0;
        int acceptedCount = 0, order = 0, magneticCalls = 0;
        int stopAfter = 0, failedKick = -1;
        bool failBisection = false;

        KOKKOS_INLINE_FUNCTION size_t count() const { return number; }
        KOKKOS_INLINE_FUNCTION bool magnetic(Vector& r, Vector&, double h, double, double) {
            ++magneticCalls;
            if (failBisection && h > .4e-9 && h < .6e-9) return false;
            r(2) += 1e7 * h;
            return true;
        }
        KOKKOS_INLINE_FUNCTION double distance(size_t j, const Vector& r) const {
            return r(2) - plane[j];
        }
        KOKKOS_INLINE_FUNCTION bool supports(size_t j, const Vector& r) const {
            return r(0) >= minimum[j] && r(0) <= maximum[j];
        }
        KOKKOS_INLINE_FUNCTION bool kick(size_t j, const Vector&, Vector& p, double t, double) {
            if (int(j) == failedKick) return false;
            p(0)       = (j + 2) * p(0) + j + 1;
            arrival[j] = t;
            return true;
        }
        KOKKOS_INLINE_FUNCTION bool visited(size_t j) const { return mask & (1ULL << j); }
        KOKKOS_INLINE_FUNCTION void markVisited(size_t j) { mask |= 1ULL << j; }
        KOKKOS_INLINE_FUNCTION bool accepted(size_t j, const Vector&, const Vector&, double) {
            ++acceptedCount;
            order = 10 * order + j + 1;
            return stopAfter && acceptedCount == stopAfter;
        }
    };

    struct Snapshot {
        Vector r = Vector(0), p = Vector(1, 0, 1);
        DevicePolicy policy;
        rf_event::Result result;
        rf_event::Diagnostic diagnostic;
    };

    Snapshot onDevice(DevicePolicy policy) {
        Kokkos::View<Snapshot*> result("RF event scheduling result", 1);
        Kokkos::parallel_for(
                "RF event scheduling", 1, KOKKOS_LAMBDA(size_t i) {
                    Snapshot state;
                    state.policy = policy;
                    state.result = rf_event::advance(
                            state.r, state.p, startTime, nominalStep, mass, state.policy,
                            &state.diagnostic);
                    result(i) = state;
                });
        return Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), result)(0);
    }

    Snapshot onHost(DevicePolicy policy) {
        Snapshot state;
        state.policy = policy;
        auto& test   = state.policy;
        rf_event::HostPolicy host;
        host.number          = test.number;
        host.magneticAdvance = [&](Vector& r, Vector& p, double h, double t, double m) {
            return test.magnetic(r, p, h, t, m);
        };
        host.planeDistance = [&](size_t j, const Vector& r) {
            return test.distance(j, r);
        };
        host.inProfile = [&](size_t j, const Vector& r) {
            return test.supports(j, r);
        };
        host.applyKick = [&](size_t j, const Vector& r, Vector& p, double t, double m) {
            return test.kick(j, r, p, t, m);
        };
        host.wasVisited = [&](size_t j) {
            return test.visited(j);
        };
        host.mark = [&](size_t j) {
            test.markVisited(j);
        };
        host.onAccepted = [&](size_t j, const Vector& r, const Vector& p, double t) {
            return test.accepted(j, r, p, t);
        };
        state.result = rf_event::advance(
                state.r, state.p, startTime, nominalStep, mass, host, &state.diagnostic);
        return state;
    }

    void expectSameExecution(const Snapshot& host, const Snapshot& device) {
        EXPECT_EQ(host.result.success, device.result.success);
        EXPECT_DOUBLE_EQ(host.result.duration, device.result.duration);
        EXPECT_EQ(host.policy.acceptedCount, device.policy.acceptedCount);
        EXPECT_EQ(host.policy.order, device.policy.order);
        EXPECT_EQ(host.policy.mask, device.policy.mask);
        EXPECT_EQ(host.policy.magneticCalls, device.policy.magneticCalls);
        EXPECT_EQ(host.diagnostic.failure, device.diagnostic.failure);
        EXPECT_EQ(host.diagnostic.gap, device.diagnostic.gap);
        // Event bisection has a 2^-40 step resolution; these tolerances are below
        // 1e-14 m in space and 3e-21 s in time for the prescribed worldline.
        for (unsigned d = 0; d < 3; ++d) {
            EXPECT_NEAR(host.r(d), device.r(d), 1e-14);
            EXPECT_DOUBLE_EQ(host.p(d), device.p(d));
        }
        for (unsigned j = 0; j < 2; ++j)
            EXPECT_NEAR(host.policy.arrival[j], device.policy.arrival[j], 3e-21);
    }
}  // namespace

class RFEventStepTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
    }
    static void TearDownTestSuite() { ippl::finalize(); }
};

TEST_F(RFEventStepTest, OrdersDistinctGapImpulsesByArrivalRatherThanDescriptorIndex) {
    const DevicePolicy policy;
    const auto state = onDevice(policy);
    ASSERT_TRUE(state.result.success);
    EXPECT_DOUBLE_EQ(state.result.duration, nominalStep);
    EXPECT_EQ(state.policy.order, 21);
    EXPECT_EQ(state.policy.acceptedCount, 2);
    EXPECT_EQ(state.policy.mask, 3);
    EXPECT_DOUBLE_EQ(state.p(0), 11);
    EXPECT_NEAR(state.r(2), .01, 1e-14);
    EXPECT_NEAR(state.policy.arrival[1], startTime + 1e-10, 3e-21);
    EXPECT_NEAR(state.policy.arrival[0], startTime + 3e-10, 3e-21);
    EXPECT_GT(state.policy.magneticCalls, state.policy.acceptedCount);
    expectSameExecution(onHost(policy), state);
}

TEST_F(RFEventStepTest, EqualArrivalSelectsLaterGapAndDoesNotRepeatAtItsPlane) {
    DevicePolicy policy;
    policy.plane[0] = policy.plane[1] = .002;
    const auto state                  = onDevice(policy);
    ASSERT_TRUE(state.result.success);
    EXPECT_EQ(state.policy.order, 2);
    EXPECT_EQ(state.policy.acceptedCount, 1);
    EXPECT_EQ(state.policy.mask, 2);
    EXPECT_DOUBLE_EQ(state.p(0), 5);
    EXPECT_NEAR(state.policy.arrival[1], startTime + 2e-10, 3e-21);
    expectSameExecution(onHost(policy), state);
}

TEST_F(RFEventStepTest, AStepStartingOnTheGapPlaneDoesNotApplyADuplicateImpulse) {
    DevicePolicy policy;
    policy.number    = 1;
    policy.plane[0]  = 0;
    const auto state = onDevice(policy);
    ASSERT_TRUE(state.result.success);
    EXPECT_EQ(state.policy.acceptedCount, 0);
    EXPECT_EQ(state.policy.mask, 0);
    EXPECT_DOUBLE_EQ(state.p(0), 1);
    EXPECT_NEAR(state.r(2), .01, 1e-14);
    expectSameExecution(onHost(policy), state);
}

TEST_F(RFEventStepTest, EndpointCrossingAppliesACompleteKickAtTheNominalEnd) {
    DevicePolicy policy;
    policy.number    = 1;
    policy.plane[0]  = .01;
    const auto state = onDevice(policy);
    ASSERT_TRUE(state.result.success);
    EXPECT_DOUBLE_EQ(state.result.duration, nominalStep);
    EXPECT_NEAR(state.r(2), .01, 1e-14);
    EXPECT_DOUBLE_EQ(state.p(0), 3);
    EXPECT_EQ(state.policy.acceptedCount, 1);
    EXPECT_EQ(state.policy.mask, 1);
    EXPECT_NEAR(state.policy.arrival[0], startTime + nominalStep, 3e-21);
    expectSameExecution(onHost(policy), state);
}

TEST_F(RFEventStepTest, FiniteProfilesClipAnEarlierPlaneWithoutDiscardingALaterGap) {
    DevicePolicy policy;
    policy.minimum[1] = .1;
    policy.maximum[1] = .2;
    // Inclusive finite-profile edge: x=0 belongs to gap zero.
    policy.minimum[0] = 0;
    policy.maximum[0] = .1;
    const auto state  = onDevice(policy);
    ASSERT_TRUE(state.result.success);
    EXPECT_EQ(state.policy.order, 1);
    EXPECT_EQ(state.policy.acceptedCount, 1);
    EXPECT_EQ(state.policy.mask, 1);
    EXPECT_DOUBLE_EQ(state.p(0), 3);
    EXPECT_DOUBLE_EQ(state.policy.arrival[1], 0);
    EXPECT_NEAR(state.policy.arrival[0], startTime + 3e-10, 3e-21);
    expectSameExecution(onHost(policy), state);
}

TEST_F(RFEventStepTest, StopsAfterACompleteKickWithoutAdvancingTheRemainder) {
    DevicePolicy policy;
    policy.stopAfter = 1;
    const auto state = onDevice(policy);
    ASSERT_TRUE(state.result.success);
    EXPECT_NEAR(state.result.duration, 1e-10, 3e-21);
    EXPECT_NEAR(state.r(2), .001, 1e-14);
    EXPECT_DOUBLE_EQ(state.p(0), 5);
    EXPECT_EQ(state.policy.order, 2);
    EXPECT_EQ(state.policy.acceptedCount, 1);
    EXPECT_EQ(state.policy.mask, 2);
    EXPECT_DOUBLE_EQ(state.policy.arrival[0], 0);
    expectSameExecution(onHost(policy), state);
}

TEST_F(RFEventStepTest, FailedBisectionReportsItsGapWithoutCommittingTrialMotion) {
    DevicePolicy policy;
    policy.failBisection = true;
    const auto state     = onDevice(policy);
    EXPECT_FALSE(state.result.success);
    EXPECT_EQ(state.diagnostic.failure, rf_event::Failure::FieldSupport);
    EXPECT_EQ(state.diagnostic.gap, 0);
    EXPECT_DOUBLE_EQ(state.diagnostic.time, startTime);
    EXPECT_DOUBLE_EQ(state.r(2), 0);
    EXPECT_DOUBLE_EQ(state.p(0), 1);
    EXPECT_EQ(state.policy.acceptedCount, 0);
    EXPECT_EQ(state.policy.mask, 0);
    expectSameExecution(onHost(policy), state);
}

TEST_F(RFEventStepTest, ARejectedLaterKickKeepsEarlierAcceptedEventsForTheCallerToRollback) {
    DevicePolicy policy;
    policy.failedKick = 0;
    const auto state  = onDevice(policy);
    EXPECT_FALSE(state.result.success);
    EXPECT_EQ(state.diagnostic.failure, rf_event::Failure::GapKick);
    EXPECT_EQ(state.diagnostic.gap, 0);
    EXPECT_NEAR(state.diagnostic.time, startTime + 3e-10, 3e-21);
    EXPECT_NEAR(state.diagnostic.position(2), .003, 1e-14);
    EXPECT_DOUBLE_EQ(state.diagnostic.momentum(0), 5);
    EXPECT_NEAR(state.r(2), .003, 1e-14);
    EXPECT_DOUBLE_EQ(state.p(0), 5);
    EXPECT_EQ(state.policy.order, 2);
    EXPECT_EQ(state.policy.acceptedCount, 1);
    EXPECT_EQ(state.policy.mask, 2);
    expectSameExecution(onHost(policy), state);
}
