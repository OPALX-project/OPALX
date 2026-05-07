/**
 * \file TestMonitor.cpp
 * \brief Unit tests for lightweight Monitor behavior.
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
 *    - bends()
 *    - getRequiredNumberOfTimeSteps()
 *    - getFieldExtend()
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
 * ---------------------------------------------------------------------------
 * Notes
 * ---------------------------------------------------------------------------
 *
 * - Monitor is field-free. The dummy field exists only to satisfy the abstract
 *   Component interface in the concrete test class.
 *
 * - These tests intentionally avoid file-output checks. They focus on behavior
 *   that can be verified without constructing a full PartBunch/ParticleContainer
 *   tracking setup.
 */

#include <gtest/gtest.h>

#include "Ippl.h"

#define private public
#include "AbsBeamline/Monitor.h"
#undef private

#include "Fields/EMField.h"

#include <memory>

// ---------------------------------------------------------------------------
// Dummy field
// ---------------------------------------------------------------------------
class DummyMonitorField : public EMField {
public:
    void scale(double) override {}
};

// ---------------------------------------------------------------------------
// Minimal concrete Monitor
// ---------------------------------------------------------------------------
class TestMonitor : public Monitor {
public:
    TestMonitor() : Monitor("test_monitor") {
        geom_.setElementLength(0.01);
    }

    explicit TestMonitor(const std::string& name) : Monitor(name) {
        geom_.setElementLength(0.01);
    }

    ElementBase* clone() const override {
        return new TestMonitor(*this);
    }

    StraightGeometry& getGeometry() override {
        return geom_;
    }

    const StraightGeometry& getGeometry() const override {
        return geom_;
    }

    EMField& getField() override {
        return field_;
    }

    const EMField& getField() const override {
        return field_;
    }

    Plane getPlane() const override {
        return plane_;
    }

    void setPlane(Plane plane) {
        plane_ = plane;
    }

    void setLength(double length) {
        geom_.setElementLength(length);
    }

private:
    StraightGeometry geom_;
    DummyMonitorField field_;
    Plane plane_ = OFF;
};


// ---------------------------------------------------------------------------
// Fixture
// ---------------------------------------------------------------------------
class MonitorTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
    }

    static void TearDownTestSuite() {
        ippl::finalize();
    }
};

// ---------------------------------------------------------------------------
// Basic API
// ---------------------------------------------------------------------------
TEST_F(MonitorTest, GetType) {
    TestMonitor monitor;

    EXPECT_EQ(monitor.getType(), ElementType::MONITOR);
}

TEST_F(MonitorTest, Bends) {
    TestMonitor monitor;

    EXPECT_FALSE(monitor.bends());
}

TEST_F(MonitorTest, RequiredNumberOfTimeStepsIsOne) {
    TestMonitor monitor;

    EXPECT_EQ(monitor.getRequiredNumberOfTimeSteps(), 1);
}

TEST_F(MonitorTest, GetFieldExtendUsesMonitorHalfLength) {
    TestMonitor monitor;

    double zBegin = 0.0;
    double zEnd = 0.0;

    monitor.getFieldExtend(zBegin, zEnd);

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
TEST_F(MonitorTest, ApplyNullParticleContainerReturnsFalse) {
    TestMonitor monitor;

    EXPECT_FALSE(monitor.apply(std::shared_ptr<ParticleContainer_t>()));
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

TEST_F(MonitorTest, FieldOnlyApplyOverloadReturnsFalseAndLeavesFieldsUnchanged) {
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

    EXPECT_FALSE(monitor.apply(R, P, 0.0, E, B));

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

