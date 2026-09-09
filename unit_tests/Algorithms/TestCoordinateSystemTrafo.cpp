#include "gtest/gtest.h"

#include "Algorithms/CoordinateSystemTrafo.h"
#include "Ippl.h"

#include <cmath>

namespace {
    constexpr double tol = 1e-12;

    using Vector3 = ippl::Vector<double, 3>;

    Quaternion rotationAroundX(double angle) {
        return Quaternion(std::cos(0.5 * angle), std::sin(0.5 * angle), 0.0, 0.0);
    }

    Quaternion rotationAroundY(double angle) {
        return Quaternion(std::cos(0.5 * angle), 0.0, std::sin(0.5 * angle), 0.0);
    }

    Quaternion rotationAroundZ(double angle) {
        return Quaternion(std::cos(0.5 * angle), 0.0, 0.0, std::sin(0.5 * angle));
    }

    void expectVectorNear(const Vector3& actual, const Vector3& expected, double tolerance = tol) {
        EXPECT_NEAR(actual(0), expected(0), tolerance);
        EXPECT_NEAR(actual(1), expected(1), tolerance);
        EXPECT_NEAR(actual(2), expected(2), tolerance);
    }
}  // namespace

class CoordinateSystemTrafoTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc    = 0;
        char** argv = nullptr;

        ippl::initialize(argc, argv);
    }

    static void TearDownTestSuite() { ippl::finalize(); }
};

TEST_F(CoordinateSystemTrafoTest, DefaultConstructorIsIdentity) {
    CoordinateSystemTrafo trafo;
    Vector3 point(1.0, -2.0, 3.0);
    Vector3 vector(-4.0, 5.0, -6.0);

    expectVectorNear(trafo.getOrigin(), Vector3(0.0));
    EXPECT_DOUBLE_EQ(trafo.getRotation().real(), 1.0);
    EXPECT_DOUBLE_EQ(trafo.getRotation()(1), 0.0);
    EXPECT_DOUBLE_EQ(trafo.getRotation()(2), 0.0);
    EXPECT_DOUBLE_EQ(trafo.getRotation()(3), 0.0);

    expectVectorNear(trafo.transformTo(point), point);
    expectVectorNear(trafo.transformFrom(point), point);
    expectVectorNear(trafo.rotateTo(vector), vector);
    expectVectorNear(trafo.rotateFrom(vector), vector);
}

TEST_F(CoordinateSystemTrafoTest, TransformRoundTrip) {
    CoordinateSystemTrafo trafo(Vector3(1.5, -2.0, 0.75), rotationAroundZ(M_PI / 2.0));
    Vector3 point(2.5, 3.0, -1.0);

    const Vector3 local = trafo.transformTo(point);

    expectVectorNear(trafo.transformFrom(local), point);
}

TEST_F(CoordinateSystemTrafoTest, RotateRoundTrip) {
    CoordinateSystemTrafo trafo(Vector3(3.0, -1.0, 2.0), rotationAroundY(M_PI / 3.0));
    Vector3 vector(1.0, -4.0, 2.0);

    const Vector3 local = trafo.rotateTo(vector);

    expectVectorNear(trafo.rotateFrom(local), vector);
}

TEST_F(CoordinateSystemTrafoTest, SpaceChargeMomentumRotationPreservesBoostAndRoundTrip) {
    // Exercise the device-view rotation used around the binned SC solve, including
    // a straight beam, a large bend, and a zero-momentum particle (gamma=1).
    using View = Kokkos::View<Vector3*>;
    View momenta("momenta", 3);
    auto host = Kokkos::create_mirror_view(momenta);
    for (double angle : {0.0, 0.7, M_PI / 2}) {
        CoordinateSystemTrafo beamToRef(Vector3(7.0, -2.0, 3.0), rotationAroundY(angle));
        const auto refToBeam = beamToRef.inverted();
        const Vector3 original[] = {
                beamToRef.rotateTo(Vector3(0.0, 0.0, 0.4)),
                beamToRef.rotateTo(Vector3(0.01, -0.02, 0.41)), Vector3(0.0)};
        for (int i = 0; i < 3; ++i) {
            host(i) = original[i];
        }
        Kokkos::deep_copy(momenta, host);
        refToBeam.rotateBunchTo(momenta, 3);
        Kokkos::deep_copy(host, momenta);
        expectVectorNear(host(0), Vector3(0.0, 0.0, 0.4));
        expectVectorNear(host(1), Vector3(0.01, -0.02, 0.41));
        expectVectorNear(host(2), Vector3(0.0));
        for (int i = 0; i < 3; ++i) {
            EXPECT_NEAR(host(i).dot(host(i)), original[i].dot(original[i]), 1.0e-14);
        }
        beamToRef.rotateBunchTo(momenta, 3);
        Kokkos::deep_copy(host, momenta);
        for (int i = 0; i < 3; ++i) {
            expectVectorNear(host(i), original[i]);
        }
    }
}

TEST_F(CoordinateSystemTrafoTest, InverseMatchesTransformFromAndRotateFrom) {
    CoordinateSystemTrafo trafo(Vector3(-1.0, 2.0, 4.0), rotationAroundX(M_PI / 4.0));
    CoordinateSystemTrafo inverse = trafo.inverted();
    Vector3 point(3.0, -2.0, 0.5);
    Vector3 vector(-1.0, 5.0, 2.0);

    expectVectorNear(inverse.transformTo(trafo.transformTo(point)), point);
    expectVectorNear(inverse.transformTo(point), trafo.transformFrom(point));
    expectVectorNear(inverse.rotateTo(vector), trafo.rotateFrom(vector));
}

TEST_F(CoordinateSystemTrafoTest, InvertInPlaceMatchesInverted) {
    CoordinateSystemTrafo trafo(Vector3(0.5, -1.5, 2.5), rotationAroundY(M_PI / 5.0));
    CoordinateSystemTrafo inverse = trafo.inverted();

    trafo.invert();

    expectVectorNear(trafo.getOrigin(), inverse.getOrigin());
    EXPECT_NEAR(trafo.getRotation().real(), inverse.getRotation().real(), tol);
    EXPECT_NEAR(trafo.getRotation()(1), inverse.getRotation()(1), tol);
    EXPECT_NEAR(trafo.getRotation()(2), inverse.getRotation()(2), tol);
    EXPECT_NEAR(trafo.getRotation()(3), inverse.getRotation()(3), tol);
}

TEST_F(CoordinateSystemTrafoTest, CompositionAppliesRightThenLeftToPoints) {
    CoordinateSystemTrafo left(Vector3(2.0, 0.0, 1.0), rotationAroundZ(M_PI / 2.0));
    CoordinateSystemTrafo right(Vector3(-1.0, 3.0, 0.5), rotationAroundX(M_PI / 6.0));
    CoordinateSystemTrafo composed = left * right;
    Vector3 point(4.0, -2.0, 1.0);

    expectVectorNear(composed.transformTo(point), left.transformTo(right.transformTo(point)));
}

TEST_F(CoordinateSystemTrafoTest, CompositionAppliesRightThenLeftToVectors) {
    CoordinateSystemTrafo left(Vector3(1.0, -2.0, 0.0), rotationAroundY(M_PI / 7.0));
    CoordinateSystemTrafo right(Vector3(0.0, 1.0, 2.0), rotationAroundZ(M_PI / 3.0));
    CoordinateSystemTrafo composed = left * right;
    Vector3 vector(3.0, -1.0, 2.0);

    expectVectorNear(composed.rotateTo(vector), left.rotateTo(right.rotateTo(vector)));
}
