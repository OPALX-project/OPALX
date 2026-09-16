#include "Algorithms/BorisMidpoint.h"
#include "Algorithms/CoordinateSystemTrafo.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "gtest/gtest.h"

using Vector = boris_midpoint::Vector;
using State = boris_midpoint::State;

class BorisMidpointTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
    }
    static void TearDownTestSuite() { ippl::finalize(); }
};

TEST_F(BorisMidpointTest, PositionDependentElectricFieldUsesTheActualHalfDrift) {
    const Vector start(2e-5, -1e-5, 0.003), initialP(0.02, -0.01, 0.1);
    const double dt = 2e-10, mass = 9.382720813e8, gradient = 1e7;
    Kokkos::View<State*> results("midpoint electric trial", 1);
    Kokkos::View<Vector*> samples("sample positions", 1), liveP("unchanged momentum", 1);
    Kokkos::deep_copy(liveP, initialP);
    Kokkos::parallel_for("electric midpoint test", 1, KOKKOS_LAMBDA(int i) {
        const Vector midpoint = boris_midpoint::halfDrift(start, liveP(i), dt);
        samples(i) = midpoint;
        const Vector electric(gradient * midpoint(0), -0.5 * gradient * midpoint(1), 0);
        results(i) = boris_midpoint::endpoint(midpoint, liveP(i), electric, Vector(0), dt, mass, 1);
    });
    const auto result = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), results);
    const auto sampled = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), samples);
    const auto preserved = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), liveP);
    const double gamma = std::sqrt(1 + dot(initialP, initialP));
    const Vector midpoint = start + (0.5 * Physics::c * dt / gamma) * initialP;
    const Vector electric(gradient * midpoint(0), -0.5 * gradient * midpoint(1), 0);
    const Vector momentum = initialP + (dt * Physics::c / mass) * electric;
    const Vector position = midpoint + (0.5 * Physics::c * dt
            / std::sqrt(1 + dot(momentum, momentum))) * momentum;
    for (unsigned d = 0; d < 3; ++d) {
        EXPECT_NEAR(sampled(0)(d), midpoint(d), 2e-18);
        EXPECT_NEAR(result(0).momentum(d), momentum(d), 2e-16);
        EXPECT_NEAR(result(0).position(d), position(d), 2e-18);
        EXPECT_DOUBLE_EQ(preserved(0)(d), initialP(d));
    }
    const double prestepPx = initialP(0) + dt * Physics::c * gradient * start(0) / mass;
    EXPECT_GT(std::abs(result(0).momentum(0) - prestepPx), 1e-8);
}

TEST_F(BorisMidpointTest, CombinedSelfAndExternalFieldsMatchIndependentRotation) {
    const Vector midpoint(0.01, -0.02, 0.03), initialP(0.2, -0.1, 0.7);
    const Vector selfE(1e5, -2e5, 3e4), externalE(4e4, 1e4, -2e4);
    const Vector selfB(0.01, 0.02, -0.03), externalB(0.1, -0.5, 0.25);
    const Vector electric = selfE + externalE, magnetic = selfB + externalB;
    const double dt = 5e-10, mass = 9.382720813e8;
    Kokkos::View<State*> results("combined field trials", 2);
    Kokkos::parallel_for("combined field test", 2, KOKKOS_LAMBDA(int i) {
        results(i) = boris_midpoint::endpoint(
                midpoint, initialP, electric, magnetic, dt, mass, i == 0 ? 1. : -1.);
    });
    const auto result = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), results);
    for (unsigned i = 0; i < 2; ++i) {
        const double charge = i == 0 ? 1. : -1.;
        const Vector impulse = (0.5 * dt * charge * Physics::c / mass) * electric;
        const Vector beforeRotation = initialP + impulse;
        const double fieldNorm = std::sqrt(dot(magnetic, magnetic));
        const Vector axis = magnetic / fieldNorm;
        const double angle = 2 * std::atan(0.5 * dt * charge * Physics::c * Physics::c
                * fieldNorm / (mass * std::sqrt(1 + dot(beforeRotation, beforeRotation))));
        const Vector momentum = std::cos(angle) * beforeRotation
                + std::sin(angle) * cross(beforeRotation, axis)
                + (1 - std::cos(angle)) * dot(beforeRotation, axis) * axis + impulse;
        const Vector position = midpoint + (0.5 * Physics::c * dt
                / std::sqrt(1 + dot(momentum, momentum))) * momentum;
        for (unsigned d = 0; d < 3; ++d) {
            EXPECT_NEAR(result(i).momentum(d), momentum(d), 4e-16);
            EXPECT_NEAR(result(i).position(d), position(d), 2e-17);
        }
    }
}

TEST_F(BorisMidpointTest, RejectedTrialUndoesOnlyUncommittedHalfDrift) {
    const Vector initialP(0.2, 0.3, 1.);
    const double initialDt = 1e-9;
    Kokkos::View<Vector*> restored("rejected half drifts", 2);
    Kokkos::parallel_for("rollback midpoint test", 2, KOKKOS_LAMBDA(int i) {
        Vector position = i == 0 ? Vector(1e-6, -2e-6, 3e-6) : Vector(1, -2, 3);
        for (unsigned trial = 0; trial < 32; ++trial) {
            const double dt = Kokkos::ldexp(initialDt, -static_cast<int>(trial));
            // Use the existing tracker's scaled-position first push. The trial
            // endpoint must not overwrite this live midpoint or its momentum.
            Vector midpoint = position / (Physics::c * dt);
            BorisPusher().push(midpoint, initialP, dt);
            midpoint *= Physics::c * dt;
            const State rejected = boris_midpoint::endpoint(
                    midpoint, initialP, Vector(1e5, -2e5, 0), Vector(0, 0.8, 0),
                    dt, 9.382720813e8, 1.);
            (void)rejected;
            position = boris_midpoint::halfDrift(midpoint, initialP, -dt);
        }
        restored(i) = position;
    });
    const auto result = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), restored);
    for (unsigned i = 0; i < 2; ++i) {
        const Vector initial = i == 0 ? Vector(1e-6, -2e-6, 3e-6) : Vector(1, -2, 3);
        const double scale = std::max(Physics::c * initialDt,
                std::sqrt(dot(initial, initial)));
        // Covers 32 scaled forward/unscaled reverse arithmetic pairs. This is
        // an absolute floating-point roundoff bound, not an integration tolerance.
        const double bound = 64 * std::numeric_limits<double>::epsilon() * scale;
        for (unsigned d = 0; d < 3; ++d)
            EXPECT_NEAR(result(i)(d), initial(d), bound);
    }
}

TEST_F(BorisMidpointTest, RotatingTheFrameRotatesBothFieldContributions) {
    const CoordinateSystemTrafo frame(
            Vector(1, -2, 0.5), Quaternion(std::cos(0.3), 0., std::sin(0.3), 0.));
    const Vector start(0.01, -0.02, 0.03), initialP(0.2, -0.1, 0.7);
    const Vector selfE(1e5, -2e5, 3e4), externalE(4e4, 1e4, -2e4);
    const Vector selfB(0.01, 0.02, -0.03), externalB(0.1, -0.5, 0.25);
    const Vector transformedR = frame.transformFrom(start), transformedP = frame.rotateFrom(initialP);
    const Vector transformedE = frame.rotateFrom(selfE) + frame.rotateFrom(externalE);
    const Vector transformedB = frame.rotateFrom(selfB) + frame.rotateFrom(externalB);
    Kokkos::View<State*> results("rotated trials", 2);
    Kokkos::parallel_for("rotated midpoint test", 2, KOKKOS_LAMBDA(int i) {
        const Vector r = i == 0 ? start : transformedR;
        const Vector p = i == 0 ? initialP : transformedP;
        const Vector e = i == 0 ? Vector(selfE + externalE) : transformedE;
        const Vector b = i == 0 ? Vector(selfB + externalB) : transformedB;
        results(i) = boris_midpoint::endpoint(
                boris_midpoint::halfDrift(r, p, 5e-10), p, e, b, 5e-10, 9.382720813e8, 1.);
    });
    const auto result = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), results);
    const Vector position = frame.transformFrom(result(0).position);
    const Vector momentum = frame.rotateFrom(result(0).momentum);
    for (unsigned d = 0; d < 3; ++d) {
        EXPECT_NEAR(result(1).position(d), position(d), 6e-16);
        EXPECT_NEAR(result(1).momentum(d), momentum(d), 4e-16);
    }
}

TEST_F(BorisMidpointTest, MagneticEnergyAndZeroChargeDriftMatchAnalyticResults) {
    const Vector midpoint(0.001, -0.002, 0.003), p(0, 0, 1.);
    const double dt = 1e-10, mass = 9.382720813e8, magnetic = 0.8;
    Kokkos::View<State*> results("magnetic and neutral trials", 3);
    Kokkos::parallel_for("magnetic and neutral test", 3, KOKKOS_LAMBDA(int i) {
        const double charge = i == 0 ? 1. : (i == 1 ? -1. : 0.);
        const Vector electric = i == 2 ? Vector(1e5, -2e5, 3e4) : Vector(0);
        results(i) = boris_midpoint::endpoint(
                midpoint, p, electric, Vector(0, magnetic, 0), dt, mass, charge);
    });
    const auto result = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), results);
    for (unsigned i = 0; i < 3; ++i) {
        const double charge = i == 0 ? 1. : (i == 1 ? -1. : 0.);
        const double angle = 2 * std::atan(0.5 * dt * charge * Physics::c * Physics::c
                * magnetic / (std::sqrt(2.) * mass));
        const Vector momentum(-std::sin(angle), 0, std::cos(angle));
        const Vector position = midpoint + (0.5 * Physics::c * dt / std::sqrt(2.)) * momentum;
        EXPECT_NEAR(dot(result(i).momentum, result(i).momentum), 1., 4e-16);
        for (unsigned d = 0; d < 3; ++d) {
            EXPECT_NEAR(result(i).momentum(d), momentum(d), 2e-16);
            EXPECT_NEAR(result(i).position(d), position(d), 4e-18);
        }
    }
}
