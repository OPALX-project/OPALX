/**
 * \file TestMonitor.cpp
 * \brief Unit tests for Monitor behavior and temporal sampling.
 *
 * This test suite validates basic Monitor behavior using:
 *  - a minimal concrete Monitor implementation
 *  - direct access to internal state where needed
 *  - a single-rank initialized IPPL/MPI environment
 *
 * ---------------------------------------------------------------------------
 * Coverage
 * ---------------------------------------------------------------------------
 *
 * 1. Basic API
 *    - getType()
 *    - getRequiredNumberOfTimeSteps()
 *    - getFieldExtent()
 *
 * 2. Collection type bookkeeping
 *    - default collection type is SPATIAL
 *    - setCollectionType() updates the stored collection mode
 *
 * 3. Geometry / containment
 *    - isInside() accepts points inside the monitor length
 *    - isInside() rejects points outside the monitor length
 *
 * 4. Safe early-return behavior
 *    - apply(nullptr) returns false
 *    - applyToReferenceParticle() returns false without a reference bunch
 *    - goOffline() is safe before a LossDataSink is created
 *
 * 5. Temporal sampling
 *    - particle timestamps and mean time agree with the reference crossing time
 *    - drifted positions match an analytic field-free snapshot
 *
 * ---------------------------------------------------------------------------
 * Notes
 * ---------------------------------------------------------------------------
 *
 * - Monitor is field-free. The dummy field exists only to satisfy the abstract
 *   Component interface in the concrete test class.
 *
 * - Sampling uses a small PartBunch and an in-memory sink, without file output
 *   or a full tracker.
 */

#include <gtest/gtest.h>

#include "Ippl.h"

#define private public
#include "AbsBeamline/Monitor.h"
#include "BeamlineGeometry/Geometry.h"
#undef private

#include "Structure/Beam.h"

#include <array>
#include <cmath>
#include <memory>

// ---------------------------------------------------------------------------
// Minimal concrete Monitor
// ---------------------------------------------------------------------------
class TestMonitor : public Monitor {
public:
    TestMonitor() : Monitor("test_monitor") { geom_.setElementLength(0.01); }

    explicit TestMonitor(const std::string& name) : Monitor(name) { geom_.setElementLength(0.01); }

    ElementBase* clone() const override { return new TestMonitor(*this); }

    Geometry& getGeometry() override { return geom_; }

    const Geometry& getGeometry() const override { return geom_; }

    Plane getPlane() const override { return plane_; }

    void setPlane(Plane plane) { plane_ = plane; }

    void setLength(double length) { geom_.setElementLength(length); }

    void attachBunch(PartBunch_t* bunch) { RefPartBunch_m = bunch; }

private:
    Geometry geom_;
    Plane plane_ = OFF;
};

// ---------------------------------------------------------------------------
// Fixture
// ---------------------------------------------------------------------------
class MonitorTest : public ::testing::Test {
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

    void TearDown() override { Monitor::statFileEntries_sm.clear(); }
};

// ---------------------------------------------------------------------------
// Basic API
// ---------------------------------------------------------------------------
TEST_F(MonitorTest, GetType) {
    TestMonitor monitor;

    EXPECT_EQ(monitor.getType(), ElementType::MONITOR);
}

TEST_F(MonitorTest, RequiredNumberOfTimeStepsIsOne) {
    TestMonitor monitor;

    EXPECT_EQ(monitor.getRequiredNumberOfTimeSteps(), 1);
}

TEST_F(MonitorTest, GetFieldExtendUsesMonitorHalfLength) {
    TestMonitor monitor;

    double zBegin = 0.0;
    double zEnd   = 0.0;

    monitor.getFieldExtent(zBegin, zEnd);

    EXPECT_DOUBLE_EQ(zBegin, -0.005);
    EXPECT_DOUBLE_EQ(zEnd, 0.005);
}

// ---------------------------------------------------------------------------
// Collection type bookkeeping
// ---------------------------------------------------------------------------
TEST_F(MonitorTest, DefaultCollectionTypeIsSpatial) {
    TestMonitor monitor;

    EXPECT_EQ(monitor.type_m, CollectionType::SPATIAL);
}

TEST_F(MonitorTest, SetCollectionTypeUpdatesCollectionMode) {
    TestMonitor monitor;

    monitor.setCollectionType(CollectionType::TEMPORAL);

    EXPECT_EQ(monitor.type_m, CollectionType::TEMPORAL);

    monitor.setCollectionType(CollectionType::SPATIAL);

    EXPECT_EQ(monitor.type_m, CollectionType::SPATIAL);
}

// ---------------------------------------------------------------------------
// Plane bookkeeping
// ---------------------------------------------------------------------------
TEST_F(MonitorTest, PlaneDefaultsToOff) {
    TestMonitor monitor;

    EXPECT_EQ(monitor.getPlane(), Monitor::OFF);
}

TEST_F(MonitorTest, PlaneCanBeSetInConcreteTestMonitor) {
    TestMonitor monitor;

    monitor.setPlane(Monitor::XY);

    EXPECT_EQ(monitor.getPlane(), Monitor::XY);
}

// ---------------------------------------------------------------------------
// Geometry / containment
// ---------------------------------------------------------------------------
TEST_F(MonitorTest, IsInsideAcceptsPointsWithinElementLength) {
    TestMonitor monitor;
    monitor.setLength(0.02);

    EXPECT_TRUE(monitor.isInside(Vector_t<double, 3>(0.0, 0.0, 0.0)));
    EXPECT_TRUE(monitor.isInside(Vector_t<double, 3>(0.0, 0.0, 0.009)));
    EXPECT_TRUE(monitor.isInside(Vector_t<double, 3>(0.0, 0.0, -0.009)));
}

TEST_F(MonitorTest, IsInsideRejectsPointsOutsideElementLength) {
    TestMonitor monitor;
    monitor.setLength(0.02);

    EXPECT_FALSE(monitor.isInside(Vector_t<double, 3>(0.0, 0.0, 0.011)));
    EXPECT_FALSE(monitor.isInside(Vector_t<double, 3>(0.0, 0.0, -0.011)));
}

// ---------------------------------------------------------------------------
// Safe early returns
// ---------------------------------------------------------------------------
TEST_F(MonitorTest, ApplyNullParticleContainerIsSafe) {
    TestMonitor monitor;

    EXPECT_NO_THROW(monitor.apply(std::shared_ptr<ParticleContainer_t>()));
}

TEST_F(MonitorTest, ApplyToReferenceParticleWithoutReferenceBunchReturnsFalse) {
    TestMonitor monitor;

    Vector_t<double, 3> R(0.0);
    Vector_t<double, 3> P(0.0);
    Vector_t<double, 3> E(1.0);
    Vector_t<double, 3> B(2.0);

    monitor.goOnline(0.0);

    EXPECT_FALSE(monitor.applyToReferenceParticle(R, P, 0.0, E, B));
}

TEST_F(MonitorTest, FieldOnlyApplyOverloadLeavesFieldsUnchanged) {
    TestMonitor monitor;

    Vector_t<double, 3> R(0.0);
    Vector_t<double, 3> P(0.0);

    Vector_t<double, 3> E(0.0);
    E(0) = 1.0;
    E(1) = 2.0;
    E(2) = 3.0;

    Vector_t<double, 3> B(0.0);
    B(0) = 4.0;
    B(1) = 5.0;
    B(2) = 6.0;

    monitor.apply(R, P, 0.0, E, B);

    EXPECT_DOUBLE_EQ(E(0), 1.0);
    EXPECT_DOUBLE_EQ(E(1), 2.0);
    EXPECT_DOUBLE_EQ(E(2), 3.0);

    EXPECT_DOUBLE_EQ(B(0), 4.0);
    EXPECT_DOUBLE_EQ(B(1), 5.0);
    EXPECT_DOUBLE_EQ(B(2), 6.0);
}

TEST_F(MonitorTest, GoOfflineWithoutLossDataSinkIsSafe) {
    TestMonitor monitor;

    EXPECT_NO_THROW(monitor.goOffline());
}

TEST_F(MonitorTest, TemporalSamplesUseReferenceCrossingTime) {
    Beam exemplar;
    Beam* beam                = Beam::find("UNNAMED_BEAM");
    const double dt           = 1.0e-10;
    const double endTime      = 10.0 * dt;
    const double midpointTime = endTime - 0.5 * dt;
    const Vector_t<double, 3> momentum(0.0, 0.0, 2.0);
    const double speed                  = Physics::c * 2.0 / std::sqrt(5.0);
    const std::array<double, 2> offsets = {-2.0e-4, 2.0e-4};

    PartBunch_t bunch(
            {Physics::q_e}, {0.00051099895}, {beam}, {offsets.size()}, 1.0, "LF2",
            opalx::spacecharge::CartesianDomainConfig3D{.periodicParticleBoundary = true});
    bunch.setdT(dt);
    bunch.setT(endTime);
    bunch.setGlobalTrackStep(1);
    auto pc = bunch.getParticleContainer();
    pc->set_sPos(1.0);  // Avoid the pre-tracking crossing guard.
    pc->createParticles(offsets.size());

    // Cover crossings before, at, and after the particle timestep endpoint.
    for (const double fraction : {0.25, 0.5, 0.75}) {
        SCOPED_TRACE(fraction);
        const double crossingTime       = midpointTime + fraction * dt;
        const double referenceMidpointZ = -fraction * speed * dt;
        auto positions                  = pc->R.getHostMirror();
        auto momenta                    = pc->P.getHostMirror();
        auto ids                        = pc->ID.getHostMirror();
        for (size_t i = 0; i < offsets.size(); ++i) {
            // Analytic endpoint state for a particle at offsets[i] at crossingTime.
            positions(i) = Vector_t<double, 3>(
                    1.0e-3, -2.0e-3, offsets[i] + speed * (endTime - crossingTime));
            momenta(i) = momentum;
            ids(i)     = i;
        }
        Kokkos::deep_copy(pc->R.getView(), positions);
        Kokkos::deep_copy(pc->P.getView(), momenta);
        Kokkos::deep_copy(pc->ID.getView(), ids);

        TestMonitor monitor;
        monitor.attachBunch(&bunch);
        monitor.setCollectionType(CollectionType::TEMPORAL);
        // Empty output name suppresses save(), retaining records for inspection.
        monitor.lossDs_m = std::make_unique<LossDataSink>();
        monitor.goOnline(0.0);
        Vector_t<double, 3> electric(0.0), magnetic(0.0);
        EXPECT_FALSE(monitor.applyToReferenceParticle(
                Vector_t<double, 3>(0.0, 0.0, referenceMidpointZ), momentum, midpointTime, electric,
                magnetic));

        ASSERT_EQ(monitor.lossDs_m->particles_m.size(), offsets.size());
        for (size_t i = 0; i < offsets.size(); ++i) {
            const auto& sample = monitor.lossDs_m->particles_m[i];
            EXPECT_NEAR(sample.getTime(), crossingTime, 1.0e-22);
            EXPECT_NEAR(sample.getX(), 1.0e-3, 1.0e-14);
            EXPECT_NEAR(sample.getY(), -2.0e-3, 1.0e-14);
            EXPECT_NEAR(sample.getZ(), offsets[i], 1.0e-14);
        }
        const auto stats = monitor.lossDs_m->computeStatistics(1);
        ASSERT_EQ(stats.size(), 1u);
        EXPECT_NEAR(stats.begin()->refTime_m, crossingTime, 1.0e-22);
        EXPECT_NEAR(stats.begin()->tmean_m, crossingTime, 1.0e-22);
        EXPECT_NEAR(stats.begin()->trms_m, 0.0, 1.0e-22);
        EXPECT_DOUBLE_EQ(bunch.getT(), endTime);
    }
}
