#include "gtest/gtest.h"

#include "BeamlineCore/DriftRep.h"
#include "Elements/OpalBeamline.h"

#include <cmath>

namespace {
    constexpr double tol = 1e-12;

    using Vector3 = ippl::Vector<double, 3>;

    Quaternion rotationAroundY(double angle) {
        return Quaternion(std::cos(0.5 * angle), 0.0, std::sin(0.5 * angle), 0.0);
    }

    void expectVectorNear(const Vector3& actual, const Vector3& expected) {
        EXPECT_NEAR(actual(0), expected(0), tol);
        EXPECT_NEAR(actual(1), expected(1), tol);
        EXPECT_NEAR(actual(2), expected(2), tol);
    }
}  // namespace

class OpalBeamlinePlacementTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;

        ippl::initialize(argc, argv);
    }

    static void TearDownTestSuite() { ippl::finalize(); }
};

TEST_F(OpalBeamlinePlacementTest, BridgeReturnsPlacedElementViewAndPreservesNominalQueries) {
    OpalBeamline beamline;
    auto drift = std::make_shared<DriftRep>("D2");
    drift->setElementLength(2.0);
    drift->setCSTrafoGlobal2Local(
            CoordinateSystemTrafo(Vector3(0.5, -1.0, 4.0), rotationAroundY(M_PI / 8.0)));
    drift->setMisalignment(CoordinateSystemTrafo(Vector3(0.25, 0.0, 0.0), Quaternion()));

    const PlacedElement placed = beamline.getPlacedElement(drift);
    const Vector3 point(1.0, 2.0, 3.0);

    expectVectorNear(
            beamline.transformToLocalCS(drift, point),
            placed.getNominalBodyTransform().transformTo(point));
    expectVectorNear(
            beamline.transformFromLocalCS(drift, point),
            placed.getNominalBodyTransform().transformFrom(point));
    expectVectorNear(
            beamline.rotateToLocalCS(drift, point),
            placed.getNominalBodyTransform().rotateTo(point));
    expectVectorNear(
            beamline.rotateFromLocalCS(drift, point),
            placed.getNominalBodyTransform().rotateFrom(point));
    expectVectorNear(
            beamline.getCSTrafoLab2Local(drift).getOrigin(),
            placed.getNominalBodyTransform().getOrigin());
    expectVectorNear(
            beamline.getNominalEntryTransform(drift).getOrigin(),
            placed.getNominalEntryTransform().getOrigin());
    expectVectorNear(
            beamline.getNominalExitTransform(drift).getOrigin(),
            placed.getNominalExitTransform().getOrigin());
    expectVectorNear(
            beamline.getMisalignment(drift).getOrigin(),
            placed.getMisalignment().getNominalToActual().getOrigin());
}

TEST_F(OpalBeamlinePlacementTest, PositionElementRelativeUsesPlacementPoseBridge) {
    OpalBeamline beamline(Vector3(0.0, 0.0, 5.0), Quaternion());
    auto drift = std::make_shared<DriftRep>("D3");
    drift->setPlacementPose(PlacementPose(
            CoordinateSystemTrafo(Vector3(1.0, 2.0, 3.0), rotationAroundY(M_PI / 10.0))));
    drift->fixPosition();

    CoordinateSystemTrafo expected = drift->getPlacementPose().getParentToNominal();
    expected *= beamline.getCSTrafoLab2Local();

    beamline.positionElementRelative(drift);

    expectVectorNear(
            beamline.getPlacedElement(drift).getNominalBodyTransform().getOrigin(),
            expected.getOrigin());
}
