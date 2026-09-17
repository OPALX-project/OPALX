#include "Algorithms/PassiveProbe.h"
#include "Ippl.h"
#include "Physics/Physics.h"
#include "gtest/gtest.h"

#include <array>
#include <cmath>
#include <limits>

namespace {
    using namespace passive_probe;

    void expectEndpoint(const Endpoint& actual, const Endpoint& expected) {
        EXPECT_DOUBLE_EQ(actual.time, expected.time);
        for (unsigned d = 0; d < 3; ++d) {
            EXPECT_DOUBLE_EQ(actual.position(d), expected.position(d));
            EXPECT_DOUBLE_EQ(actual.momentum(d), expected.momentum(d));
        }
    }

    void expectState(const State& actual, const State& expected) {
        expectEndpoint(actual.previous, expected.previous);
        EXPECT_EQ(actual.turns, expected.turns);
        EXPECT_EQ(actual.initialized, expected.initialized);
        EXPECT_EQ(actual.armed, expected.armed);
    }
}  // namespace

class PassiveProbeTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
    }
    static void TearDownTestSuite() { ippl::finalize(); }
};

TEST_F(PassiveProbeTest, DriftCrossingIsExactAndDoesNotMutateAcceptedEndpoints) {
    const Plane plane{Vector(0.3, -0.4, 0.7), Vector(0, 0, 3), 1e-9};
    const Vector p(0.02, -0.01, 0.3);
    const Vector velocity     = Physics::c * p / std::sqrt(1 + dot(p, p));
    const double crossingTime = 2e-8, dt = 5e-9;
    const Endpoint before{plane.origin - 0.4 * dt * velocity, p, crossingTime - 0.4 * dt};
    const Endpoint after{plane.origin + 0.6 * dt * velocity, p, crossingTime + 0.6 * dt};
    const Endpoint savedBefore = before, savedAfter = after;
    State state;
    Sample sample;
    EXPECT_EQ(update(plane, state, before, sample), Status::Initialized);
    EXPECT_TRUE(state.armed);
    EXPECT_EQ(update(plane, state, after, sample), Status::Crossing);
    EXPECT_EQ(sample.turn, 1u);
    EXPECT_EQ(state.turns, 1u);
    EXPECT_NEAR(sample.fraction, 0.4, 3e-16);
    EXPECT_NEAR(sample.crossing.time, crossingTime, 4e-24);
    for (unsigned d = 0; d < 3; ++d) {
        EXPECT_NEAR(sample.crossing.position(d), plane.origin(d), 1e-16);
        EXPECT_DOUBLE_EQ(sample.crossing.momentum(d), p(d));
    }
    expectEndpoint(before, savedBefore);
    expectEndpoint(after, savedAfter);
}

TEST_F(PassiveProbeTest, InitialPlaneAndOppositeCrossingDoNotCount) {
    const Plane plane;
    State state;
    Sample sample;
    EXPECT_EQ(update(plane, state, {Vector(0), Vector(0, 0, 1), 0}, sample), Status::Initialized);
    EXPECT_EQ(
            update(plane, state, {Vector(0, 0, 1), Vector(0, 0, 1), 1}, sample), Status::Advanced);
    EXPECT_EQ(
            update(plane, state, {Vector(0, 0, -1), Vector(0, 0, -1), 2}, sample),
            Status::Advanced);
    EXPECT_EQ(state.turns, 0u);
    EXPECT_TRUE(state.armed);
    EXPECT_EQ(
            update(plane, state, {Vector(0, 0, -0.5), Vector(0, 0, 1), 3}, sample),
            Status::Advanced);
    EXPECT_EQ(update(plane, state, {Vector(0), Vector(0, 0, 1), 4}, sample), Status::Crossing);
    EXPECT_DOUBLE_EQ(sample.fraction, 1.0);
    EXPECT_DOUBLE_EQ(sample.crossing.time, 4.0);
    EXPECT_EQ(update(plane, state, {Vector(0), Vector(0, 0, 1), 5}, sample), Status::Advanced);
    EXPECT_EQ(state.turns, 1u);
    EXPECT_EQ(
            update(plane, state, {Vector(0, 0, -1), Vector(0, 0, 1), 6}, sample), Status::Advanced);
    EXPECT_EQ(
            update(plane, state, {Vector(0, 0, 1), Vector(0, 0, 1), 7}, sample), Status::Crossing);
    EXPECT_EQ(state.turns, 2u);
}

TEST_F(PassiveProbeTest, ArmingTolerancePreventsNearPlaneJitterAndChecksMomentumDirection) {
    const Plane plane{Vector(0), Vector(0, 0, 1), 1e-3};
    State state;
    Sample sample;
    EXPECT_EQ(update(plane, state, {Vector(0), Vector(0, 0, 1), 0}, sample), Status::Initialized);
    for (unsigned i = 1; i <= 10; ++i) {
        const Endpoint endpoint{Vector(0, 0, (i % 2 ? -1 : 1) * 5e-4), Vector(0, 0, 1), double(i)};
        EXPECT_EQ(update(plane, state, endpoint, sample), Status::Advanced);
    }
    EXPECT_EQ(state.turns, 0u);
    EXPECT_EQ(
            update(plane, state, {Vector(0, 0, -0.1), Vector(0, 0, -1), 11}, sample),
            Status::Advanced);
    EXPECT_EQ(
            update(plane, state, {Vector(0, 0, 0.1), Vector(0, 0, -1), 12}, sample),
            Status::Advanced);
    EXPECT_FALSE(state.armed);
    EXPECT_EQ(state.turns, 0u);
}

TEST_F(PassiveProbeTest, ObliqueTranslatedPlaneUsesNormalizedSignedDistance) {
    const Plane plane{Vector(2, -3, 1), Vector(2, -1, 2), 0.1};
    const Vector unitNormal = plane.normal / 3.0;
    State state;
    Sample sample;
    const Endpoint before{plane.origin - 0.2 * unitNormal, 0.7 * unitNormal, 1};
    const Endpoint after{plane.origin + 0.3 * unitNormal, 0.7 * unitNormal, 2};
    EXPECT_EQ(update(plane, state, before, sample), Status::Initialized);
    EXPECT_TRUE(state.armed);
    EXPECT_EQ(update(plane, state, after, sample), Status::Crossing);
    EXPECT_NEAR(sample.crossing.time, 1.4, 1e-15);
    for (unsigned d = 0; d < 3; ++d)
        EXPECT_NEAR(sample.crossing.position(d), plane.origin(d), 5e-16);
}

TEST_F(PassiveProbeTest, ConstantElectricAccelerationHasSecondOrderObservationError) {
    // Exact relativistic longitudinal trajectory for dp/(mc)/dt = acceleration.
    // This tests the observation interpolation independently of any integrator.
    constexpr double p0 = 0.2, acceleration = 2e7, crossingTime = 1e-7;
    const auto endpoint = [&](double t) {
        const double p = p0 + acceleration * t;
        const double z =
                Physics::c / acceleration * (std::sqrt(1 + p * p) - std::sqrt(1 + p0 * p0));
        return Endpoint{Vector(0.2, -0.1, z), Vector(0, 0, p), crossingTime + t};
    };
    double previousTimeError = 0, previousMomentumError = 0;
    for (const double dt : {4e-9, 2e-9, 1e-9, 5e-10}) {
        State state;
        Sample sample;
        ASSERT_EQ(update(Plane{}, state, endpoint(-0.4 * dt), sample), Status::Initialized);
        ASSERT_EQ(update(Plane{}, state, endpoint(0.6 * dt), sample), Status::Crossing);
        const double timeError     = std::abs(sample.crossing.time - crossingTime);
        const double momentumError = std::abs(sample.crossing.momentum(2) - p0);
        EXPECT_GT(timeError, 0);
        EXPECT_GT(momentumError, 0);
        EXPECT_NEAR(sample.crossing.position(2), 0, 1e-16);
        EXPECT_NEAR(momentumError, acceleration * timeError, 2e-15);
        if (previousTimeError != 0) {
            EXPECT_GT(previousTimeError / timeError, 3.8);
            EXPECT_LT(previousTimeError / timeError, 4.1);
            EXPECT_GT(previousMomentumError / momentumError, 3.8);
            EXPECT_LT(previousMomentumError / momentumError, 4.1);
        }
        previousTimeError     = timeError;
        previousMomentumError = momentumError;
    }
}

TEST_F(PassiveProbeTest, DisabledObservationLeavesHistoryAndSampleUntouched) {
    State state{{Vector(1), Vector(2), 3}, 7, true, true};
    Sample sample{{Vector(4), Vector(5), 6}, 7, 0.3};
    const State savedState   = state;
    const Sample savedSample = sample;
    const double nan         = std::numeric_limits<double>::quiet_NaN();
    EXPECT_EQ(
            update(Plane{}, state, {Vector(nan), Vector(nan), nan}, sample, false),
            Status::Disabled);
    expectState(state, savedState);
    expectEndpoint(sample.crossing, savedSample.crossing);
    EXPECT_EQ(sample.turn, savedSample.turn);
    EXPECT_DOUBLE_EQ(sample.fraction, savedSample.fraction);
}

TEST_F(PassiveProbeTest, RejectsInvalidPlanesAndEndpointsWithoutPartialWrites) {
    const State initial{{Vector(0, 0, -1), Vector(0, 0, 1), 1}, 3, true, true};
    const Sample initialSample{{Vector(4), Vector(5), 6}, 7, 0.3};
    const Endpoint valid{Vector(0, 0, 1), Vector(0, 0, 1), 2};
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    const auto check = [&](const Plane& plane, const Endpoint& endpoint, State state) {
        const State savedState = state;
        Sample sample          = initialSample;
        EXPECT_EQ(update(plane, state, endpoint, sample), Status::InvalidInput);
        // NaN inputs are only in plane/accepted; history in this test is finite.
        expectState(state, savedState);
        expectEndpoint(sample.crossing, initialSample.crossing);
        EXPECT_EQ(sample.turn, initialSample.turn);
        EXPECT_DOUBLE_EQ(sample.fraction, initialSample.fraction);
    };
    check({Vector(0), Vector(0), 1e-9}, valid, initial);
    check({Vector(nan, 0, 0), Vector(0, 0, 1), 1e-9}, valid, initial);
    check({Vector(0), Vector(0, inf, 1), 1e-9}, valid, initial);
    check({Vector(0), Vector(1e308, 0, 1), 1e-9}, valid, initial);
    for (double tolerance : {0., -1., nan, inf})
        check({Vector(0), Vector(0, 0, 1), tolerance}, valid, initial);
    check(Plane{}, {Vector(0, nan, 1), valid.momentum, 2}, initial);
    check(Plane{}, {valid.position, Vector(0, inf, 1), 2}, initial);
    check(Plane{}, {valid.position, valid.momentum, nan}, initial);
    // Finite endpoints can still overflow the time difference or interpolation.
    State largeTime         = initial;
    largeTime.previous.time = -1e308;
    check(Plane{}, {valid.position, valid.momentum, 1e308}, largeTime);
    State largeCoordinate                = initial;
    largeCoordinate.previous.position(0) = -1e308;
    check(Plane{}, {Vector(1e308, 0, 1), valid.momentum, 2}, largeCoordinate);
}

TEST_F(PassiveProbeTest, RoundedEqualTimesAdvanceHistoryAndAllowSpatiallyBracketedCrossings) {
    constexpr double time = 4.57e-6, substep = 1e-22;
    ASSERT_GT(substep, 0);
    ASSERT_EQ(time + substep, time);
    ASSERT_EQ(time + 2 * substep, time);
    State state;
    Sample sample;
    const Endpoint initial{Vector(0, 0, -2e-9), Vector(0, 0, 0.2), time - 1e-15};
    const Endpoint near{Vector(0, 0, -3e-14), Vector(0, 0, 0.7), time};
    const Endpoint coalesced{Vector(0, 0, -1e-14), Vector(0, 0, 0.3), time + substep};
    const Endpoint after{Vector(0, 0, 1e-14), Vector(0, 0, 0.5), time + 2 * substep};
    ASSERT_EQ(update(Plane{}, state, initial, sample), Status::Initialized);
    ASSERT_TRUE(state.armed);
    ASSERT_EQ(update(Plane{}, state, near, sample), Status::Advanced);
    EXPECT_EQ(update(Plane{}, state, coalesced, sample), Status::Advanced);
    expectEndpoint(state.previous, coalesced);
    EXPECT_EQ(state.turns, 0u);
    EXPECT_EQ(sample.turn, 0u);
    ASSERT_EQ(update(Plane{}, state, after, sample), Status::Crossing);
    EXPECT_DOUBLE_EQ(sample.fraction, 0.5);
    EXPECT_DOUBLE_EQ(sample.crossing.time, time);
    EXPECT_DOUBLE_EQ(sample.crossing.position(2), 0);
    EXPECT_DOUBLE_EQ(sample.crossing.momentum(2), 0.4);
    EXPECT_EQ(sample.turn, 1u);
    expectEndpoint(state.previous, after);
    // Re-observing the same represented endpoint cannot duplicate the crossing.
    EXPECT_EQ(update(Plane{}, state, after, sample), Status::Advanced);
    EXPECT_EQ(state.turns, 1u);
    const Endpoint later{
            Vector(0, 0, 2e-14), Vector(0, 0, 0.5),
            std::nextafter(time, std::numeric_limits<double>::infinity())};
    EXPECT_EQ(update(Plane{}, state, later, sample), Status::Advanced);
    expectEndpoint(state.previous, later);
}

TEST_F(PassiveProbeTest, DecreasingTimeAndCounterOverflowLeaveHistoryUntouched) {
    State state{{Vector(0, 0, -1), Vector(0, 0, 1), 1}, 3, true, true};
    Sample sample;
    const State initial = state;
    for (double time : {std::nextafter(1., 0.), 0., -1.}) {
        EXPECT_EQ(
                update(Plane{}, state, {Vector(0, 0, 1), Vector(0, 0, 1), time}, sample),
                Status::DecreasingTime);
        expectState(state, initial);
    }
    state.turns      = std::numeric_limits<std::uint64_t>::max();
    const State full = state;
    EXPECT_EQ(
            update(Plane{}, state, {Vector(0, 0, 1), Vector(0, 0, 1), 2}, sample),
            Status::TurnOverflow);
    expectState(state, full);
    EXPECT_EQ(sample.turn, 0u);
}

TEST_F(PassiveProbeTest, DeviceAndHostObservationsAgree) {
    constexpr unsigned count = 5;
    Kokkos::View<State*> states("passive histories", count);
    Kokkos::View<Sample*> samples("passive crossings", count);
    Kokkos::View<Status*> statuses("passive statuses", count);
    const Plane plane{Vector(0.2, -0.4, 0.5), Vector(0, 0, 2), 1e-9};
    Kokkos::parallel_for(
            "passive accepted endpoints", count, KOKKOS_LAMBDA(unsigned i) {
                State state;
                Sample sample;
                const Endpoint before{plane.origin - Vector(0, 0, 0.1), Vector(0, 0, 1), 4.57e-6};
                const Endpoint after{
                        plane.origin + Vector(0, 0, 0.1 * (i + 1)), Vector(0, 0, 1),
                        4.57e-6 + (i == 3 ? 3.25e-22 : 1e-9)};
                update(plane, state, before, sample);
                statuses(i) = update(plane, state, after, sample, i != 4);
                states(i)   = state;
                samples(i)  = sample;
            });
    const auto stateHost  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), states);
    const auto sampleHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), samples);
    const auto statusHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), statuses);
    for (unsigned i = 0; i < count; ++i) {
        State state;
        Sample sample;
        const Endpoint before{plane.origin - Vector(0, 0, 0.1), Vector(0, 0, 1), 4.57e-6};
        const Endpoint after{
                plane.origin + Vector(0, 0, 0.1 * (i + 1)), Vector(0, 0, 1),
                4.57e-6 + (i == 3 ? 3.25e-22 : 1e-9)};
        update(plane, state, before, sample);
        EXPECT_EQ(statusHost(i), update(plane, state, after, sample, i != 4));
        if (i == 3) {
            EXPECT_EQ(statusHost(i), Status::Crossing);
            EXPECT_DOUBLE_EQ(sampleHost(i).crossing.time, before.time);
        }
        expectState(stateHost(i), state);
        expectEndpoint(sampleHost(i).crossing, sample.crossing);
        EXPECT_EQ(sampleHost(i).turn, sample.turn);
        EXPECT_DOUBLE_EQ(sampleHost(i).fraction, sample.fraction);
    }
}
