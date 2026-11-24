#include "../KeywordObjects.h"
#include "Algorithms/PartData.h"
#include "Algorithms/TrackReferenceStep.h"
#include "BeamlineCore/MonitorRep.h"
#include "BeamlineCore/MultipoleRep.h"
#include "Elements/OpalBeamline.h"
#include "PartBunch/PartBunch.h"
#include "Structure/Beam.h"
#include "Track/FreezeRun.h"
#include "Utilities/OpalException.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <array>
#include <memory>
#include <string>

using track_reference::State;
using Vector = Vector_t<double, 3>;

namespace {
    class CountingMonitor : public MonitorRep {
    public:
        CountingMonitor() : MonitorRep("M") { getGeometry().setElementLength(1); }
        void initialise(PartBunch_t*) override {}
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
    auto beamObject = opalx::test::makeKeywordObject<Beam>("BEAM", "BEAM");
    Beam& beam      = *beamObject;
    PartBunch_t bunch(
            {1.}, {1.}, {opalx::run::describeBeamForBunch(beam)}, {0}, 1.,
            opalx::spacecharge::CartesianDomainConfig3D{});
    OpalBeamline lattice;
    auto monitor = std::make_shared<CountingMonitor>();
    monitor->setCSTrafoGlobal2Local(CoordinateSystemTrafo(Vector(0), Quaternion()));
    monitor->fixPosition();
    lattice.addElement(monitor, bunch);
    lattice.prepareSections();
    lattice.buildElementMapping(defaultElementMappingGrid, defaultElementReach);
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

TEST_F(TrackReferenceStepTest, LineOrderMakesReferenceFieldSumsIndependentOfAddresses) {
    auto beamObject = opalx::test::makeKeywordObject<Beam>("BEAM", "BEAM");
    Beam& beam      = *beamObject;
    PartBunch_t bunch(
            {1.}, {1.}, {opalx::run::describeBeamForBunch(beam)}, {0}, 1.,
            opalx::spacecharge::CartesianDomainConfig3D{});
    // Line order sums first, second, third; the second lattice holds them in reverse address
    // order, so a sum in address order would round differently.
    const std::array<double, 3> dipoles{0.1, 0.2, 0.3};
    ASSERT_NE((dipoles[0] + dipoles[1]) + dipoles[2], (dipoles[2] + dipoles[1]) + dipoles[0]);
    const std::array<std::string, 3> names{"first", "second", "third"};
    const PartData reference(1., 9.382720813e8, 1e6);
    const State start{Vector(0, 0, 0.5), Vector(0, 0, 0.1)};
    std::array<std::array<double, 12>, 2> snapshots{};
    for (unsigned arrangement = 0; arrangement < 2; ++arrangement) {
        std::array<std::shared_ptr<MultipoleRep>, 3> magnets;
        for (auto& magnet : magnets)
            magnet = std::make_shared<MultipoleRep>("unassigned");
        std::sort(magnets.begin(), magnets.end());
        if (arrangement == 1) std::reverse(magnets.begin(), magnets.end());
        OpalBeamline lattice;
        for (unsigned role = 0; role < 3; ++role) {
            auto& magnet = *magnets[role];
            magnet.getGeometry().setElementLength(1.);
            lattice.addElement(magnets[role], bunch);
            magnet.setName(names[role]);
            magnet.setNormalComponent(0, dipoles[role]);
            magnet.setAperture(ApertureType::RECTANGULAR, {2., 2.});
            magnet.setCSTrafoGlobal2Local(CoordinateSystemTrafo());
            magnet.fixPosition();
        }
        lattice.prepareSections();
        lattice.buildElementMapping(defaultElementMappingGrid, defaultElementReach);
        Vector electric(0.), magnetic(0.);
        lattice.getFieldAt(start.position, start.momentum, 0., electric, magnetic);
        EXPECT_EQ(magnetic(1), (dipoles[0] + dipoles[1]) + dipoles[2]);
        const auto end =
                track_reference::advanceInBeamline(lattice, reference, start, 1e-9, 1e-9, false);
        for (unsigned d = 0; d < 3; ++d) {
            snapshots[arrangement][d]     = electric(d);
            snapshots[arrangement][3 + d] = magnetic(d);
            snapshots[arrangement][6 + d] = end.position(d);
            snapshots[arrangement][9 + d] = end.momentum(d);
        }
    }
    for (unsigned i = 0; i < snapshots[0].size(); ++i)
        EXPECT_EQ(snapshots[0][i], snapshots[1][i]) << "field/end-state component " << i;
}

TEST_F(TrackReferenceStepTest, MisalignmentMovesTheElementForTheBunchOnly) {
    auto beamObject = opalx::test::makeKeywordObject<Beam>("BEAM", "BEAM");
    PartBunch_t bunch(
            {1.}, {1.}, {opalx::run::describeBeamForBunch(*beamObject)}, {2}, 1.,
            opalx::spacecharge::CartesianDomainConfig3D{});
    // A quadrupole over lab z in [0, 1) m with a 4 mm bore, shifted by DX = 3 mm and DZ = 0.3 m.
    // In the misaligned frame a lab point has local x - 3 mm and local z - 0.3 m.
    auto quad = std::make_shared<MultipoleRep>("Q");
    quad->getGeometry().setElementLength(1.);
    quad->setNormalComponent(1, 2.);
    quad->setAperture(ApertureType::RECTANGULAR, {4e-3, 4e-3});
    quad->setCSTrafoGlobal2Local(CoordinateSystemTrafo());
    quad->setMisalignment(CoordinateSystemTrafo(Vector(3e-3, 0., 0.3), Quaternion()));
    quad->fixPosition();
    OpalBeamline lattice;
    lattice.addElement(quad, bunch);
    lattice.prepareSections();
    lattice.buildElementMapping(defaultElementMappingGrid, defaultElementReach);

    // Reference side: the element at its nominal place. At z = 0.15 m only the nominal body is
    // there; at x = -2 mm the nominal bore holds the point and the misaligned one does not.
    const Vector forward(0., 0., 1.);
    for (const Vector& point : {Vector(1e-3, 0., 0.15), Vector(-2e-3, 0., 0.5)}) {
        EXPECT_EQ(lattice.getElements(point).size(), 1u);
        Vector electric(0.), magnetic(0.), expectedE(0.), expectedB(0.);
        lattice.getFieldAt(point, forward, 0., electric, magnetic);
        quad->applyToReferenceParticle(point, forward, 0., expectedE, expectedB);
        EXPECT_NE(expectedB(1), 0.);
        for (unsigned d = 0; d < 3; ++d)
            EXPECT_NEAR(magnetic(d), expectedB(d), 1e-15) << d;
    }
    EXPECT_NO_THROW(lattice.getElementMapping().checkReferencePoint(Vector(-2e-3, 0., 0.5), 0.5));
    EXPECT_THROW(
            lattice.getElementMapping().checkReferencePoint(Vector(5e-3, 0., 0.5), 0.5),
            OpalException);
    // On the nominal axis the reference feels no quadrupole field.
    const PartData reference(1., 9.382720813e8, 1e6);
    const auto end = track_reference::advanceInBeamline(
            lattice, reference, State{Vector(0., 0., 0.5), forward}, 1e-10, 1e-10, false);
    EXPECT_FALSE(end.hitMaterial);
    EXPECT_EQ(end.momentum(0), 0.);
    EXPECT_EQ(end.momentum(1), 0.);

    // Bunch side: the misaligned element. The particle on the nominal axis gets the field at
    // misaligned x = -3 mm; the one at x = -2 mm is outside the misaligned bore.
    const std::array<Vector, 2> positions{Vector(0., 0., 0.5), Vector(-2e-3, 0., 0.5)};
    auto pc = bunch.getParticleContainer(0);
    pc->createParticles(positions.size());
    auto R = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->R.getView());
    auto P = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->P.getView());
    for (size_t i = 0; i < positions.size(); ++i) {
        R(i) = positions[i];
        P(i) = forward;
    }
    Kokkos::deep_copy(pc->R.getView(), R);
    Kokkos::deep_copy(pc->P.getView(), P);
    bunch.resetPcActive();
    bunch.setTimeStep(1e-12);

    lattice.applyExternalFields(bunch);
    const auto B = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->B.getView());
    const CoordinateSystemTrafo misaligned =
            quad->getMisalignment() * quad->getCSTrafoGlobal2Local();
    Vector expectedE(0.), expectedB(0.);
    quad->apply(misaligned.transformTo(positions[0]), forward, 0., expectedE, expectedB);
    EXPECT_NE(expectedB(1), 0.);
    for (unsigned d = 0; d < 3; ++d)
        EXPECT_NEAR(B(0)(d), expectedB(d), 1e-15) << d;

    EXPECT_EQ(lattice.applyElementApertures(bunch), 1u);
    const auto invalid =
            Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->InvalidMask.getView());
    EXPECT_FALSE(invalid(0));
    EXPECT_TRUE(invalid(1));
}
