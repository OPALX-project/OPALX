// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#include <limits>
#include "Algorithms/TrackReferenceStep.h"
#include "Fields/CyclotronRFProfile.h"
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
        double speed       = 1e7;

        KOKKOS_INLINE_FUNCTION size_t count() const { return number; }
        KOKKOS_INLINE_FUNCTION bool magnetic(Vector& r, Vector&, double h, double, double) {
            ++magneticCalls;
            if (failBisection && h > .4e-9 && h < .6e-9) return false;
            r(2) += speed * h;
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

    // In-memory flat voltage profile exercises the actual proton RF energy
    // impulse, without an input deck, external data or an OPALX process.
    template <class Profile>
    struct ProtonGapPolicy {
        Profile profile;
        bool stopAtKick = false, visitedGap = false;
        KOKKOS_INLINE_FUNCTION size_t count() const { return 1; }
        KOKKOS_INLINE_FUNCTION bool magnetic(Vector& r, Vector& p, double h, double, double) {
            r += (Physics::c * h / Kokkos::sqrt(1 + dot(p, p))) * p;
            return true;
        }
        KOKKOS_INLINE_FUNCTION double distance(size_t, const Vector& r) const {
            return r(2) - .003;
        }
        KOKKOS_INLINE_FUNCTION bool supports(size_t, const Vector& r) const {
            return Kokkos::abs(r(0)) <= 1;
        }
        KOKKOS_INLINE_FUNCTION bool kick(size_t, const Vector&, Vector& p, double t, double m) {
            return CyclotronRFKick{500000, 1, 0, 0, 1}.apply(profile, .5, t, m, p);
        }
        KOKKOS_INLINE_FUNCTION bool visited(size_t) const { return visitedGap; }
        KOKKOS_INLINE_FUNCTION void markVisited(size_t) { visitedGap = true; }
        KOKKOS_INLINE_FUNCTION bool accepted(size_t, const Vector&, const Vector&, double) const {
            return stopAtKick;
        }
    };

    struct ProtonResult {
        Vector position, momentum;
        rf_event::Result result;
        bool kicked;
    };

    // Namespace scope also keeps the result type and launch operator accessible
    // to NVCC: generated GTest TestBody methods are private.
    struct TerminalProtonKernel {
        Kokkos::View<double**> profile;
        Kokkos::View<ProtonResult*> output;
        Vector start, initialP;
        double duration;
        unsigned long long frameSteps = 1;
        KOKKOS_INLINE_FUNCTION void operator()(size_t i) const {
            Vector r = start, p = initialP;
            ProtonGapPolicy<decltype(profile)> policy{profile};
            const auto result = rf_event::advance(
                    r, p, startTime, duration, mass, policy, nullptr, true, frameSteps);
            output(i) = {r, p, result, policy.visitedGap};
        }
    };

    struct Snapshot {
        Vector r = Vector(0), p = Vector(1, 0, 1);
        DevicePolicy policy;
        rf_event::Result result;
        rf_event::Diagnostic diagnostic;
    };

    Snapshot onDevice(
            DevicePolicy policy, double dt = nominalStep, bool terminalEndpoint = false,
            Vector initialR = Vector(0), unsigned long long frameSteps = 1) {
        Kokkos::View<Snapshot*> result("RF event scheduling result", 1);
        Kokkos::parallel_for(
                "RF event scheduling", 1, KOKKOS_LAMBDA(size_t i) {
                    Snapshot state;
                    state.policy = policy;
                    state.r      = initialR;
                    state.result = rf_event::advance(
                            state.r, state.p, startTime, dt, mass, state.policy, &state.diagnostic,
                            terminalEndpoint, frameSteps);
                    result(i) = state;
                });
        return Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), result)(0);
    }

    Snapshot onHost(
            DevicePolicy policy, double dt = nominalStep, bool terminalEndpoint = false,
            Vector initialR = Vector(0), unsigned long long frameSteps = 1) {
        Snapshot state;
        state.r      = initialR;
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
                state.r, state.p, startTime, dt, mass, host, &state.diagnostic, terminalEndpoint,
                frameSteps);
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

TEST_F(RFEventStepTest, ReferenceTerminalDurationIncludesACoincidentRoundoffShortParticle) {
    DevicePolicy policy;
    policy.number    = 1;
    policy.plane[0]  = .003;
    policy.stopAfter = 1;
    // Replay the actual host event duration, including subtraction roundoff in
    // dt-remaining, rather than choosing a convenient hand-written final time.
    const auto reference = onHost(policy);
    ASSERT_TRUE(reference.result.success);
    ASSERT_EQ(reference.policy.acceptedCount, 1);
    const double duration = reference.result.duration;
    policy.stopAfter      = 0;
    const Vector particleStart(
            0, 0, policy.plane[0] - policy.speed * duration - 4 * Kokkos::epsilon_v<double>);
    const auto strict = onDevice(policy, duration, false, particleStart);
    ASSERT_TRUE(strict.result.success);
    ASSERT_LT(strict.r(2), policy.plane[0]);
    EXPECT_EQ(strict.policy.acceptedCount, 0);  // Exposes the original missed kick.
    const auto terminal = onDevice(policy, duration, true, particleStart);
    ASSERT_TRUE(terminal.result.success);
    EXPECT_DOUBLE_EQ(terminal.result.duration, duration);
    EXPECT_EQ(terminal.policy.acceptedCount, 1);
    EXPECT_DOUBLE_EQ(terminal.p(0), reference.p(0));  // Whole impulse, never a fraction.
    EXPECT_NEAR(terminal.policy.arrival[0], startTime + duration, 3e-21);
    EXPECT_DOUBLE_EQ(terminal.r(2), strict.r(2));  // No position projection or extra drift.
    expectSameExecution(onHost(policy, duration, true, particleStart), terminal);
}

TEST_F(RFEventStepTest, TerminalEndpointIncludesAccumulatedFrameArithmeticOnlyWithinItsBudget) {
    DevicePolicy policy;
    policy.number                           = 1;
    policy.plane[0]                         = .01;
    constexpr unsigned long long frameSteps = 1468;
    // The actual RF-stop fixture accumulates about 150 fm of signed plane
    // deficit over 1468 frame updates. A one-step arithmetic budget cannot
    // represent that history, even when the final step itself is very short.
    const Vector accumulated(0, 0, -150e-15);
    const auto singleStep = onDevice(policy, nominalStep, true, accumulated);
    EXPECT_TRUE(singleStep.result.success);
    EXPECT_EQ(singleStep.policy.acceptedCount, 0);
    const auto terminal = onDevice(policy, nominalStep, true, accumulated, frameSteps);
    EXPECT_TRUE(terminal.result.success);
    EXPECT_EQ(terminal.policy.acceptedCount, 1);
    EXPECT_DOUBLE_EQ(terminal.p(0), 3);
    EXPECT_DOUBLE_EQ(terminal.result.duration, nominalStep);
    expectSameExecution(onHost(policy, nominalStep, true, accumulated, frameSteps), terminal);
    const auto ordinary = onDevice(policy, nominalStep, false, accumulated, frameSteps);
    EXPECT_TRUE(ordinary.result.success);
    EXPECT_EQ(ordinary.policy.acceptedCount, 0);

    // Four times the accumulated arithmetic allowance is a resolvably later
    // crossing and must not acquire a terminal kick at the reference stop.
    const Vector later(0, 0, -256. * static_cast<double>(frameSteps) * Kokkos::epsilon_v<double>);
    const auto miss = onDevice(policy, nominalStep, true, later, frameSteps);
    EXPECT_TRUE(miss.result.success);
    EXPECT_LT(miss.r(2), policy.plane[0]);  // A later crossing is still upstream.
    EXPECT_EQ(miss.policy.acceptedCount, 0);
    EXPECT_DOUBLE_EQ(miss.p(0), 1);
    expectSameExecution(onHost(policy, nominalStep, true, later, frameSteps), miss);
}

TEST_F(RFEventStepTest, RejectsUnrepresentableTerminalAccumulationWithoutChangingOrdinarySteps) {
    DevicePolicy policy;
    policy.number            = 1;
    policy.plane[0]          = .01;
    constexpr auto excessive = std::numeric_limits<unsigned long long>::max();
    const auto rejected      = onDevice(policy, nominalStep, true, Vector(0), excessive);
    EXPECT_FALSE(rejected.result.success);
    EXPECT_DOUBLE_EQ(rejected.result.duration, 0);
    EXPECT_EQ(rejected.diagnostic.failure, rf_event::Failure::InvalidStep);
    EXPECT_EQ(rejected.policy.magneticCalls, 0);
    EXPECT_EQ(rejected.policy.acceptedCount, 0);
    EXPECT_DOUBLE_EQ(rejected.r(2), 0);
    EXPECT_DOUBLE_EQ(rejected.p(0), 1);
    expectSameExecution(onHost(policy, nominalStep, true, Vector(0), excessive), rejected);

    const Vector late(0, 0, -150e-15);
    const auto ordinary = onDevice(policy, nominalStep, false, late, excessive);
    EXPECT_TRUE(ordinary.result.success);
    EXPECT_EQ(ordinary.policy.acceptedCount, 0);
    EXPECT_DOUBLE_EQ(ordinary.p(0), 1);
    expectSameExecution(onHost(policy, nominalStep, false, late, excessive), ordinary);
}

TEST_F(RFEventStepTest, TerminalEndpointRoundoffBoundScalesWithLaboratoryCoordinates) {
    DevicePolicy policy;
    policy.number     = 1;
    policy.plane[0]   = .01;
    policy.minimum[0] = -2e6;
    policy.maximum[0] = 2e6;
    const Vector start(1e6, 0, -4e6 * Kokkos::epsilon_v<double>);
    const auto state = onDevice(policy, nominalStep, true, start);
    ASSERT_TRUE(state.result.success);
    ASSERT_LT(state.r(2), policy.plane[0]);
    EXPECT_EQ(state.policy.acceptedCount, 1);
    EXPECT_DOUBLE_EQ(state.result.duration, nominalStep);
    EXPECT_DOUBLE_EQ(state.p(0), 3);
    expectSameExecution(onHost(policy, nominalStep, true, start), state);
}

TEST_F(RFEventStepTest, TerminalEndpointDoesNotKickLaterRecedingOrUnsupportedParticles) {
    DevicePolicy policy;
    policy.number   = 1;
    policy.plane[0] = .01;
    const Vector late(0, 0, -128 * Kokkos::epsilon_v<double>);
    const auto later = onDevice(policy, nominalStep, true, late);
    ASSERT_TRUE(later.result.success);
    EXPECT_EQ(later.policy.acceptedCount, 0);
    EXPECT_DOUBLE_EQ(later.p(0), 1);
    expectSameExecution(onHost(policy, nominalStep, true, late), later);

    policy.minimum[0] = .1;
    policy.maximum[0] = .2;
    const Vector rounded(0, 0, -4 * Kokkos::epsilon_v<double>);
    const auto unsupported = onDevice(policy, nominalStep, true, rounded);
    EXPECT_TRUE(unsupported.result.success);
    EXPECT_EQ(unsupported.policy.acceptedCount, 0);
    expectSameExecution(onHost(policy, nominalStep, true, rounded), unsupported);

    policy.minimum[0]   = -1;
    policy.maximum[0]   = 1;
    policy.plane[0]     = 0;
    policy.speed        = -2 * Kokkos::epsilon_v<double> / nominalStep;
    const auto receding = onDevice(policy, nominalStep, true, rounded);
    EXPECT_TRUE(receding.result.success);
    EXPECT_EQ(receding.policy.acceptedCount, 0);
    expectSameExecution(onHost(policy, nominalStep, true, rounded), receding);

    policy.speed          = 0;
    const auto stationary = onDevice(policy, nominalStep, true, rounded);
    EXPECT_TRUE(stationary.result.success);
    EXPECT_EQ(stationary.policy.acceptedCount, 0);
    expectSameExecution(onHost(policy, nominalStep, true, rounded), stationary);
    const auto onPlane = onDevice(policy, nominalStep, true, Vector(0));
    EXPECT_TRUE(onPlane.result.success);
    EXPECT_EQ(onPlane.policy.acceptedCount, 0);
}

TEST_F(RFEventStepTest, TerminalKickMeetsReferenceEnergyForOneOrTwoDistributedProtons) {
    Kokkos::View<double**> profile("flat RF voltage profile", 2, 3);
    auto hostProfile = Kokkos::create_mirror_view(profile);
    for (unsigned i = 0; i < 2; ++i) {
        hostProfile(i, 0) = i;
        hostProfile(i, 1) = 1;
        hostProfile(i, 2) = 0;
    }
    Kokkos::deep_copy(profile, hostProfile);
    constexpr double initialEnergy = 72.8e6, targetEnergy = 73e6;
    const double gamma = 1 + initialEnergy / mass;
    const Vector initialP(0, 0, std::sqrt(gamma * gamma - 1));
    const double duration =
            track_reference::collectiveTerminalStep(ippl::Comm->getCommunicator(), [&] {
                Vector r(0), p = initialP;
                ProtonGapPolicy<decltype(hostProfile)> reference{hostProfile, true};
                const auto result =
                        rf_event::advance(r, p, startTime, nominalStep, mass, reference);
                if (!result.success || !reference.visitedGap)
                    throw std::runtime_error("Reference proton missed RF gap");
                return result.duration;
            });
    const double speed                      = Physics::c * initialP(2) / gamma;
    constexpr unsigned long long frameSteps = 1468;
    const Vector start(0, 0, .003 - speed * duration - 150e-15);
    for (int total : {1, 2}) {
        const int local = total / ippl::Comm->size()
                          + (ippl::Comm->rank() < total % ippl::Comm->size() ? 1 : 0);
        Kokkos::View<ProtonResult*> output("terminal proton states", local);
        Kokkos::parallel_for(
                "terminal proton RF kick", local,
                TerminalProtonKernel{profile, output, start, initialP, duration, frameSteps});
        const auto result = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), output);
        int localKicks    = 0;
        for (int i = 0; i < local; ++i) {
            EXPECT_TRUE(result(i).result.success);
            EXPECT_TRUE(result(i).kicked);
            localKicks += result(i).kicked ? 1 : 0;
            EXPECT_DOUBLE_EQ(result(i).result.duration, duration);
            const double energy =
                    (std::sqrt(1 + dot(result(i).momentum, result(i).momentum)) - 1) * mass;
            EXPECT_GE(energy, targetEnergy);
            // Eight gamma ulps expressed in eV, below 2 micro-eV for a proton;
            // this checks the full analytic gain, not just crossing the target.
            EXPECT_NEAR(
                    energy, initialEnergy + 500000 * std::cos(startTime + duration),
                    8 * Kokkos::epsilon_v<double> * mass);
            EXPECT_DOUBLE_EQ(result(i).momentum(0), 0);
            EXPECT_DOUBLE_EQ(result(i).momentum(1), 0);
        }
        int globalKicks = 0;
        MPI_Allreduce(
                &localKicks, &globalKicks, 1, MPI_INT, MPI_SUM, ippl::Comm->getCommunicator());
        EXPECT_EQ(globalKicks, total);
    }

    // A physical particle can exceed the target before the reference's common
    // final time. Its policy must finish the supplied duration after that kick,
    // rather than applying an individual energy stop and leaving its position
    // behind the bunch clock.
    Kokkos::View<ProtonResult*> continuing("off-reference proton state", 1);
    Kokkos::parallel_for(
            "off-reference proton RF kick", 1,
            TerminalProtonKernel{profile, continuing, Vector(0), initialP, nominalStep});
    const auto end = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), continuing);
    EXPECT_TRUE(end(0).result.success);
    EXPECT_TRUE(end(0).kicked);
    EXPECT_DOUBLE_EQ(end(0).result.duration, nominalStep);
    EXPECT_GT(end(0).position(2), .003);
    EXPECT_GE((std::sqrt(1 + dot(end(0).momentum, end(0).momentum)) - 1) * mass, targetEnergy);
}
