#include "Algorithms/DefaultVisitor.h"
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
