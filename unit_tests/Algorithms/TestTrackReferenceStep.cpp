#include "Algorithms/DefaultVisitor.h"
#include "Algorithms/DirectedTurnCounter.h"
#include "Algorithms/PartData.h"
#include "Algorithms/TrackReferenceStep.h"
#include "BeamlineCore/MonitorRep.h"
#include "Beamlines/FlaggedBeamline.h"
#include "Elements/OpalBeamline.h"
#include "PartBunch/PartBunch.h"
#include "Structure/Beam.h"
#include "Utilities/OpalException.h"
#include "gtest/gtest.h"

using track_reference::State;
using Vector = Vector_t<double, 3>;

namespace {
    class CountingMonitor : public MonitorRep {
    public:
        CountingMonitor() : MonitorRep("M") { getGeometry().setElementLength(1); }
        ElementBase* clone() const override { return new CountingMonitor(*this); }
        void initialise(PartBunch_t*) override {}
        void finalise() override {}
        bool applyToReferenceParticle(
                const Vector&, const Vector&, const double&, Vector&, Vector&) override {
            ++calls;
            return false;
        }
        int calls = 0;
    };
}  // namespace

class TrackReferenceStepTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
        gmsg = new Inform(nullptr, -1);
    }
    static void TearDownTestSuite() {
        delete gmsg;
        gmsg = nullptr;
        ippl::finalize();
    }
};

TEST_F(TrackReferenceStepTest, DriftAndMidpointTime) {
    const State start{Vector(0.2, -0.1, 0.3), Vector(0.3, 0.4, 1.2)};
    const double dt = 1e-10, time = 2e-9;
    int calls  = 0;
    auto field = [&](const Vector& r, const Vector& p, double t, Vector&, Vector&) {
        ++calls;
        EXPECT_DOUBLE_EQ(t, time - 0.5 * dt);
        for (unsigned d = 0; d < 3; ++d)
            EXPECT_NEAR(
                    r(d),
                    start.position(d) + 0.5 * Physics::c * dt * p(d) / std::sqrt(1 + dot(p, p)),
                    2e-16);
        return true;
    };
    const auto end = track_reference::advance(start, dt, time, 9.382720813e8, 1, field);
    EXPECT_EQ(calls, 1);
    EXPECT_TRUE(end.hitMaterial);
    for (unsigned d = 0; d < 3; ++d) {
        EXPECT_NEAR(
                end.position(d),
                start.position(d)
                        + Physics::c * dt * start.momentum(d)
                                  / std::sqrt(1 + dot(start.momentum, start.momentum)),
                2e-16);
        EXPECT_DOUBLE_EQ(end.momentum(d), start.momentum(d));
    }
    EXPECT_THROW(track_reference::advance(start, 0, time, 1, 1, field), std::invalid_argument);
}

TEST_F(TrackReferenceStepTest, MagneticBorisAngleAndEnergyForBothCharges) {
    const State start{Vector(0), Vector(0, 0, 1)};
    const double mass = 9.382720813e8, dt = 1e-10, fieldStrength = 0.8;
    for (double charge : {-1., 1.}) {
        const auto end = track_reference::advance(
                start, dt, dt, mass, charge,
                [&](const Vector&, const Vector&, double, Vector&, Vector& b) {
                    b(1) = fieldStrength;
                    return false;
                });
        const double angle = 2
                             * std::atan(
                                     0.5 * dt * charge * Physics::c * Physics::c * fieldStrength
                                     / (std::sqrt(2.) * mass));
        EXPECT_NEAR(end.momentum(0), -std::sin(angle), 2e-16);
        EXPECT_NEAR(end.momentum(2), std::cos(angle), 2e-16);
        EXPECT_NEAR(dot(end.momentum, end.momentum), 1., 3e-16);
    }
}

TEST_F(TrackReferenceStepTest, BeamlineTrialsDoNotWriteMonitors) {
    Beam beam;
    PartBunch_t bunch(
            {1.}, {1.}, {&beam}, {0}, 1., "LF2", opalx::spacecharge::CartesianDomainConfig3D{});
    FlaggedBeamline line;
    DefaultVisitor visitor(line, false, false);
    OpalBeamline lattice;
    CountingMonitor monitor;
    monitor.setCSTrafoGlobal2Local(CoordinateSystemTrafo(Vector(0), Quaternion()));
    monitor.fixPosition();
    lattice.visit(monitor, visitor, bunch);
    lattice.prepareSections();
    auto* stored = dynamic_cast<CountingMonitor*>(lattice.getElements().begin()->get());
    ASSERT_NE(stored, nullptr);
    const PartData reference(1., 9.382720813e8, 1e6);
    const State start{Vector(0, 0, 0.2), Vector(0, 0, 1)};
    const auto trial =
            track_reference::advanceInBeamline(lattice, reference, start, 1e-10, 1e-10, false);
    EXPECT_EQ(stored->calls, 0);
    const auto committed =
            track_reference::advanceInBeamline(lattice, reference, start, 1e-10, 1e-10, true);
    EXPECT_EQ(stored->calls, 1);
    for (unsigned d = 0; d < 3; ++d) {
        EXPECT_DOUBLE_EQ(trial.position(d), committed.position(d));
        EXPECT_DOUBLE_EQ(trial.momentum(d), committed.momentum(d));
    }
}

TEST_F(TrackReferenceStepTest, TerminalDriftInRotatedTranslatedPlane) {
    const Vector origin(1, 2, 3), normal(0.6, 0, 0.8);
    DirectedTurnCounter counter(origin, normal);
    const Vector start = origin - 0.003 * normal;
    counter.update(start, normal);
    const double dt = 1e-10, speed = Physics::c / std::sqrt(2.);
    int calls    = 0;
    auto advance = [&](double h) {
        ++calls;
        return State{start + h * speed * normal, normal};
    };
    const auto duration = counter.terminalStep(1, start, dt, advance);
    EXPECT_NEAR(duration, 0.003 / speed, 1e-12 * dt);
    EXPECT_EQ(counter.count(), 0);  // Trials cannot consume a turn.
    const auto end            = advance(duration);
    const Vector displacement = end.position - origin;
    EXPECT_GE(dot(displacement, normal), 0);
    EXPECT_LT(dot(displacement, normal), 1e-9);
    EXPECT_TRUE(counter.update(end.position, end.momentum));
    EXPECT_GT(calls, 1);
}

TEST_F(TrackReferenceStepTest, SkipsLaunchReverseEarlierAndDistantCrossings) {
    const Vector zero(0), forward(0, 0, 1), start(0, 0, -0.003);
    DirectedTurnCounter counter(zero, forward);
    const double dt = 1e-10;
    int calls       = 0;
    auto advance    = [&](double h) {
        ++calls;
        return State{start + Vector(0, 0, h * 1e8), forward};
    };
    EXPECT_EQ(counter.terminalStep(1, start, dt, advance), dt);  // Not armed.
    EXPECT_EQ(calls, 0);
    counter.update(start, forward);
    EXPECT_EQ(counter.terminalStep(2, start, dt, advance), dt);  // Not terminal turn.
    EXPECT_EQ(counter.terminalStep(1, Vector(0, 0, -1), dt, advance), dt);
    EXPECT_EQ(calls, 0);
    EXPECT_EQ(
            counter.terminalStep(
                    1, start, dt,
                    [&](double) {
                        return State{Vector(0, 0, 0.001), -forward};
                    }),
            dt);
    EXPECT_TRUE(counter.update(zero, forward));
    counter.update(start, forward);
    EXPECT_LT(counter.terminalStep(2, start, dt, advance), dt);
}

TEST_F(TrackReferenceStepTest, RejectsInvalidOrUnresolvedTerminalEvent) {
    const Vector start(0, 0, -0.001), forward(0, 0, 1);
    DirectedTurnCounter counter(Vector(0), forward);
    counter.update(start, forward);
    auto discontinuous = [&](double h) {
        return State{Vector(0, 0, h < 5e-11 ? -0.001 : 0.001), forward};
    };
    EXPECT_THROW(counter.terminalStep(1, start, 0, discontinuous), std::invalid_argument);
    EXPECT_THROW(counter.terminalStep(1, start, 1e-10, discontinuous), std::runtime_error);
    EXPECT_THROW(
            counter.terminalStep(
                    1, start, 1e-10,
                    [&](double h) {
                        return State{start + Vector(0, 0, h * 1e8), forward, true};
                    }),
            std::runtime_error);
    EXPECT_EQ(counter.count(), 0);
}

TEST_F(TrackReferenceStepTest, ExactEndpointAndNoCrossing) {
    const Vector start(0, 0, -0.001), forward(0, 0, 1);
    DirectedTurnCounter counter(Vector(0), forward);
    counter.update(start, forward);
    constexpr double dt = 1e-10;
    EXPECT_DOUBLE_EQ(
            counter.terminalStep(
                    1, start, dt,
                    [&](double h) {
                        return State{Vector(0, 0, 0.001 * (h / dt - 1)), forward};
                    }),
            dt);
    EXPECT_DOUBLE_EQ(
            counter.terminalStep(
                    1, start, dt,
                    [&](double) {
                        return State{start, forward};
                    }),
            dt);
}

TEST_F(TrackReferenceStepTest, ThreeTurnsInUniformMagnetLocalizesOnlyTheLastStep) {
    const double mass = 9.382720813e8, omega = Physics::c / std::sqrt(2.);
    const double dt = 2 * std::acos(-1.) / omega / 127.3;
    State state{Vector(1, 0, 0), Vector(0, 0, 1)};
    DirectedTurnCounter counter(state.position, state.momentum);
    double accumulatedAngle = 0;
    unsigned clipped        = 0;
    for (unsigned step = 0; step < 400 && counter.count() < 3; ++step) {
        const auto advance = [&](double h) {
            return track_reference::advance(
                    state, h, h, mass, 1,
                    [&](const Vector&, const Vector&, double, Vector&, Vector& b) {
                        b(1) = mass / Physics::c;
                        return false;
                    });
        };
        const double accepted = counter.terminalStep(3, state.position, dt, advance);
        if (accepted < dt) {
            ++clipped;
            EXPECT_EQ(counter.count(), 2);
        }
        accumulatedAngle += 2 * std::atan(0.5 * omega * accepted);
        state = advance(accepted);
        counter.update(state.position, state.momentum);
    }
    EXPECT_EQ(counter.count(), 3);
    EXPECT_EQ(clipped, 1);
    EXPECT_NEAR(accumulatedAngle, 6 * std::acos(-1.), 2e-12);
    EXPECT_NEAR(state.position(0), 1., 2e-12);
    EXPECT_GE(state.position(2), 0);
    EXPECT_LT(state.position(2), 1e-9);
    EXPECT_NEAR(dot(state.momentum, state.momentum), 1., 1e-13);
}

TEST_F(TrackReferenceStepTest, CollectiveSearchOnlyRunsOnRootAndPropagatesFailures) {
    const auto communicator = ippl::Comm->getCommunicator();
    int calls               = 0;
    EXPECT_DOUBLE_EQ(
            track_reference::collectiveTerminalStep(
                    communicator,
                    [&] {
                        ++calls;
                        return 3e-11;
                    }),
            3e-11);
    EXPECT_EQ(calls, ippl::Comm->rank() == 0 ? 1 : 0);
    EXPECT_THROW(
            track_reference::collectiveTerminalStep(
                    communicator,
                    []() -> double {
                        throw std::runtime_error("trial failure");
                    }),
            OpalException);
    EXPECT_THROW(
            track_reference::collectiveTerminalStep(
                    communicator,
                    []() -> double {
                        throw OpalException("test field", "outside support");
                    }),
            OpalException);
    EXPECT_THROW(
            track_reference::collectiveTerminalStep(
                    communicator,
                    [] {
                        return 0.;
                    }),
            OpalException);
    EXPECT_THROW(
            track_reference::collectiveTerminalStep(
                    communicator,
                    []() -> double {
                        throw 1;
                    }),
            OpalException);
}
