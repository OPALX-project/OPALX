#include "gtest/gtest.h"

#include "BeamlineCore/DriftRep.h"
#include "BeamlineGeometry/ElementGeometry.h"
#include "BeamlineGeometry/Misalignment.h"
#include "BeamlineGeometry/PlacedElement.h"
#include "BeamlineGeometry/PlacementPose.h"
#include "BeamlineGeometry/Port.h"
#include "BeamlineGeometry/SupportPlacement.h"

#include <cmath>

namespace {
    constexpr double tol = 1e-12;

    using Vector3 = ippl::Vector<double, 3>;

    Quaternion rotationAroundZ(double angle) {
        return Quaternion(std::cos(0.5 * angle), 0.0, 0.0, std::sin(0.5 * angle));
    }

    void expectVectorNear(const Vector3& actual, const Vector3& expected) {
        EXPECT_NEAR(actual(0), expected(0), tol);
        EXPECT_NEAR(actual(1), expected(1), tol);
        EXPECT_NEAR(actual(2), expected(2), tol);
    }
}  // namespace

class PlacedElementTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;

        ippl::initialize(argc, argv);
    }

    static void TearDownTestSuite() { ippl::finalize(); }
};

TEST_F(PlacedElementTest, ActualBodyTransformAppliesCorrectionAfterNominalPlacement) {
    PlacementPose nominal(
            CoordinateSystemTrafo(Vector3(1.0, 2.0, 3.0), rotationAroundZ(M_PI / 2.0)));
    Misalignment correction(
            CoordinateSystemTrafo(Vector3(0.5, -0.25, 0.0), CoordinateSystemTrafo().getRotation()));
    PlacedElement placed(nullptr, nominal, correction);
    Vector3 point(2.0, 0.0, -1.0);

    const Vector3 actual   = placed.getActualBodyTransform().transformTo(point);
    const Vector3 expected = correction.getNominalToActual().transformTo(
            nominal.getParentToNominal().transformTo(point));

    expectVectorNear(actual, expected);
}

TEST_F(PlacedElementTest, EntryExitAndSupportTransformsAreDerivedFromBodyFrame) {
    PlacementPose nominal(
            CoordinateSystemTrafo(Vector3(0.0, 0.0, 1.0), rotationAroundZ(M_PI / 4.0)));
    Misalignment correction(CoordinateSystemTrafo(Vector3(0.0, 0.0, 0.5), Quaternion()));
    ElementGeometry geometry(
            Port("entry", CoordinateSystemTrafo(Vector3(0.0, 0.0, -0.2), Quaternion())),
            Port("body", CoordinateSystemTrafo(Vector3(0.0, 0.0, 0.0), Quaternion())),
            Port("exit", CoordinateSystemTrafo(Vector3(0.0, 0.0, 0.8), Quaternion())));
    SupportPlacement support(CoordinateSystemTrafo(Vector3(0.1, 0.0, 0.0), Quaternion()));
    PlacedElement placed(nullptr, nominal, correction, geometry, support);
    Vector3 point(0.5, -0.5, 2.0);

    const CoordinateSystemTrafo actualBody = placed.getActualBodyTransform();
    expectVectorNear(
            placed.getEntryTransform().transformTo(point),
            geometry.getEntry().getBodyToPort().transformTo(actualBody.transformTo(point)));
    expectVectorNear(
            placed.getBodyTransform().transformTo(point),
            geometry.getBody().getBodyToPort().transformTo(actualBody.transformTo(point)));
    expectVectorNear(
            placed.getExitTransform().transformTo(point),
            geometry.getExit().getBodyToPort().transformTo(actualBody.transformTo(point)));
    expectVectorNear(
            placed.getSupportTransform().transformTo(point),
            support.getBodyToSupport().transformTo(actualBody.transformTo(point)));
}

TEST_F(PlacedElementTest, DefaultGeometryExposesNamedPorts) {
    ElementGeometry geometry;

    EXPECT_EQ(geometry.getEntry().getName(), "entry");
    EXPECT_EQ(geometry.getBody().getName(), "body");
    EXPECT_EQ(geometry.getExit().getName(), "exit");
}

TEST_F(PlacedElementTest, ElementBaseBridgeUsesExistingPlacementAndEdgeState) {
    DriftRep drift("D1");
    drift.setElementLength(1.5);
    drift.setCSTrafoGlobal2Local(
            CoordinateSystemTrafo(Vector3(1.0, 2.0, 3.0), rotationAroundZ(M_PI / 6.0)));
    drift.setMisalignment(CoordinateSystemTrafo(Vector3(0.1, 0.2, 0.3), Quaternion()));

    const PlacedElement placed = drift.getPlacedElement();

    expectVectorNear(placed.getNominalBodyTransform().getOrigin(), Vector3(1.0, 2.0, 3.0));
    expectVectorNear(placed.getNominalEntryTransform().getOrigin(), Vector3(1.0, 2.0, 3.0));
    EXPECT_NEAR(placed.getNominalExitTransform().getOrigin()(2), 4.5, tol);
    EXPECT_EQ(placed.getGeometry().getEntry().getName(), "entry");
    EXPECT_EQ(placed.getGeometry().getExit().getName(), "exit");
}

TEST_F(PlacedElementTest, ElementBaseExposesExplicitStraightElementPortContract) {
    DriftRep drift("D2");
    drift.setElementLength(1.5);

    const Port entry = drift.getEntryPort();
    const Port body  = drift.getBodyPort();
    const Port exit  = drift.getExitPort();

    EXPECT_EQ(entry.getName(), "entry");
    EXPECT_EQ(body.getName(), "body");
    EXPECT_EQ(exit.getName(), "exit");
    expectVectorNear(entry.getBodyToPort().getOrigin(), Vector3(0.0, 0.0, 0.0));
    expectVectorNear(body.getBodyToPort().getOrigin(), Vector3(0.0, 0.0, 0.0));
    expectVectorNear(exit.getBodyToPort().getOrigin(), Vector3(0.0, 0.0, 1.5));
}
