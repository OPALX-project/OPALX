#include "Algorithms/OneTurnMap.h"
#include "Algorithms/PartData.h"
#include "Elements/OpalBeamline.h"
#include "Utilities/OpalException.h"
#include "gtest/gtest.h"

namespace {
    OneTurnMap::Settings settings(unsigned points = 512) {
        OneTurnMap::Settings s;
        s.momentum = 1;
        s.dt       = 2 * std::acos(-1.0) * std::sqrt(2.0) / Physics::c / points;
        s.maxPath  = 9;
        s.maxSteps = 2 * points;
        return s;
    }
    OneTurnMap::Advance uniformField(
            const ExternalFieldRayTracker& tracker, const Vector_t<double, 3>& field) {
        return [&tracker, field](
                       const OneTurnMap::State& ray, double dt,
                       std::vector<OneTurnMap::Step>* accepted) {
            auto step = tracker.step(ray, dt, [field](const auto&, auto& e, auto& b) {
                e = Vector_t<double, 3>(0.0);
                b = field;
                return false;
            });
            if (accepted) accepted->push_back(step);
            return step.end;
        };
    }
}  // namespace

TEST(OneTurnMap, ReturnsAfterFullTurnForEitherCharge) {
    OpalBeamline beamline;
    for (double charge : {-1.0, 1.0}) {
        PartData reference(charge, 938272081.3, 938272081.3);
        ExternalFieldRayTracker tracker(
                beamline, reference, ExternalFieldRayTracker::IntegrationMethod::RK4);
        OneTurnMap map(
                uniformField(
                        tracker,
                        Vector_t<double, 3>(0, -reference.getM() / Physics::c / charge, 0)),
                CoordinateSystemTrafo(), settings());
        const auto result = map({0, 0, 0, 0});
        EXPECT_NEAR(result.ray.pathLength, 2 * std::acos(-1.0), 3e-8);
        EXPECT_GT(result.steps, 500u);  // Neither launch nor half-turn crossing.
        EXPECT_NEAR(result.coordinates[0], 0, 3e-8);
        EXPECT_NEAR(result.coordinates[1], 0, 3e-8);
        EXPECT_NEAR(result.sectionResidual, 0, 1e-12);
        EXPECT_GT(result.ray.momentum(2), 0);
    }
}

TEST(OneTurnMap, FixedSectionJacobianUsesMomentumNotSlope) {
    OpalBeamline beamline;
    PartData reference(1, 938272081.3, 938272081.3);
    ExternalFieldRayTracker tracker(
            beamline, reference, ExternalFieldRayTracker::IntegrationMethod::DOP853);
    OneTurnMap map(
            uniformField(tracker, Vector_t<double, 3>(0, -reference.getM() / Physics::c, 0)),
            CoordinateSystemTrafo(), settings(128));
    const auto matrix = map.jacobian({0.01, 0.02, 0.03, 0}, {1e-4, 1e-4, 1e-4, 1e-4});
    // Uniform B: a full transverse circle, free vertical drift.
    for (unsigned i = 0; i < 4; ++i)
        for (unsigned j = 0; j < 4; ++j) {
            const double expected = i == j ? 1.0 : (i == 2 && j == 3 ? 2 * std::acos(-1.0) : 0.0);
            EXPECT_NEAR(matrix[i][j], expected, 2e-8) << i << ',' << j;
        }
}

TEST(OneTurnMap, RejectsInvalidLaunchAndBudgets) {
    OpalBeamline beamline;
    PartData reference(1, 938272081.3, 938272081.3);
    ExternalFieldRayTracker tracker(beamline, reference);
    auto s     = settings();
    s.maxSteps = 2;
    OneTurnMap map(tracker, CoordinateSystemTrafo(), s);
    EXPECT_THROW(map({0, 1, 0, 0}), OpalException);
    EXPECT_THROW(map({0, 0, 0, 0}), OpalException);  // Drift never returns.
    EXPECT_THROW(map.jacobian({0, 0, 0, 0}, {0, 1e-4, 1e-4, 1e-4}), OpalException);
    s.dt = -1;
    EXPECT_THROW(OneTurnMap(tracker, CoordinateSystemTrafo(), s), OpalException);
}

TEST(OneTurnMap, LostTrialsPropagate) {
    OneTurnMap map(
            [](const auto&, double, auto*) -> OneTurnMap::State {
                throw OpalException("test", "material interception");
            },
            CoordinateSystemTrafo(), settings());
    EXPECT_THROW(map({0, 0, 0, 0}), OpalException);
}

TEST(OneTurnMap, TranslatedRotatedSectionAndPathLimit) {
    OpalBeamline beamline;
    PartData reference(1, 938272081.3, 938272081.3);
    ExternalFieldRayTracker tracker(
            beamline, reference, ExternalFieldRayTracker::IntegrationMethod::DOP853);
    const double angle = 0.7;
    CoordinateSystemTrafo section(
            Vector_t<double, 3>(3, -2, 1),
            Quaternion(std::cos(angle / 2), 0, std::sin(angle / 2), 0));
    const auto fields = uniformField(
            tracker, section.rotateFrom(Vector_t<double, 3>(0, -reference.getM() / Physics::c, 0)));
    auto s = settings(128);
    const OneTurnMap map(fields, section, s);
    const auto result = map({0.01, 0.02, 0.03, 0});
    EXPECT_NEAR(result.coordinates[0], 0.01, 1e-11);
    EXPECT_NEAR(result.coordinates[1], 0.02, 1e-11);
    EXPECT_NEAR(result.coordinates[2], 0.03, 1e-11);
    EXPECT_NEAR(result.sectionResidual, 0, 1e-12);
    s.maxPath = 3;
    EXPECT_THROW(OneTurnMap(fields, section, s)({0, 0, 0, 0}), OpalException);
}

TEST(OneTurnMap, PerturbedRaysHaveIndependentArrivalTimes) {
    OpalBeamline beamline;
    PartData reference(1, 938272081.3, 938272081.3);
    ExternalFieldRayTracker tracker(
            beamline, reference, ExternalFieldRayTracker::IntegrationMethod::DOP853);
    OneTurnMap map(
            [&](const auto& ray, double dt, auto* accepted) {
                auto step = tracker.step(ray, dt, [&](const auto& state, auto& e, auto& b) {
                    e = Vector_t<double, 3>(0.0);
                    b = Vector_t<double, 3>(
                            0, -reference.getM() / Physics::c * (1 + 0.1 * state.position(0)), 0);
                    return false;
                });
                if (accepted) accepted->push_back(step);
                return step.end;
            },
            CoordinateSystemTrafo(), settings(128));
    const auto a = map({-0.01, 0, 0, 0}), b = map({0.01, 0, 0, 0});
    EXPECT_GT(std::abs(a.ray.time - b.ray.time), 1e-12);
    EXPECT_NEAR(a.sectionResidual, 0, 1e-12);
    EXPECT_NEAR(b.sectionResidual, 0, 1e-12);
}

TEST(OneTurnMap, ReturnTimeConvergesWithRK4Refinement) {
    OpalBeamline beamline;
    PartData reference(1, 938272081.3, 938272081.3);
    ExternalFieldRayTracker tracker(
            beamline, reference, ExternalFieldRayTracker::IntegrationMethod::RK4);
    const auto fields =
            uniformField(tracker, Vector_t<double, 3>(0, -reference.getM() / Physics::c, 0));
    const double exact  = 2 * std::acos(-1.0) * std::sqrt(2.0) / Physics::c;
    const double coarse = std::abs(
            OneTurnMap(fields, CoordinateSystemTrafo(), settings(64))({0, 0, 0, 0}).ray.time
            - exact);
    const double fine = std::abs(
            OneTurnMap(fields, CoordinateSystemTrafo(), settings(128))({0, 0, 0, 0}).ray.time
            - exact);
    EXPECT_GT(coarse / fine, 10);
    EXPECT_LT(coarse / fine, 22);
}

namespace {
    // Exact circular reference with independently specified harmonic transverse
    // motion. The integer winding is absent from the endpoint eigenphases.
    OneTurnMap::Advance harmonicRing(double qx, double qy) {
        return [qx, qy](const auto& ray, double dt, std::vector<OneTurnMap::Step>* accepted) {
            auto state = [&](double h) {
                const double a = ray.time, b = a + h;
                const Vector_t<double, 3> normal(std::cos(a), 0, std::sin(a));
                const Vector_t<double, 3> base(std::cos(a) - 1, 0, std::sin(a));
                const Vector_t<double, 3> offset = ray.position - base;
                const double x = dot(offset, normal), px = dot(ray.momentum, normal);
                const double y = ray.position(1), py = ray.momentum(1);
                const double nx  = x * std::cos(qx * h) + px * std::sin(qx * h) / qx;
                const double npx = -qx * x * std::sin(qx * h) + px * std::cos(qx * h);
                const double ny  = y * std::cos(qy * h) + py * std::sin(qy * h) / qy;
                const double npy = -qy * y * std::sin(qy * h) + py * std::cos(qy * h);
                OneTurnMap::State end;
                end.time       = b;
                end.pathLength = ray.pathLength + h;
                end.position =
                        Vector_t<double, 3>((1 + nx) * std::cos(b) - 1, ny, (1 + nx) * std::sin(b));
                const double tangent = std::sqrt(1 - npx * npx - npy * npy);
                end.momentum         = Vector_t<double, 3>(
                        npx * std::cos(b) - tangent * std::sin(b), npy,
                        npx * std::sin(b) + tangent * std::cos(b));
                return end;
            };
            const auto end = state(dt);
            if (accepted) {
                OneTurnMap::Step interval;
                interval.midpoint = state(dt / 2);
                interval.end      = end;
                interval.duration = dt;
                accepted->push_back(interval);
            }
            return end;
        };
    }
}  // namespace

TEST(OneTurnMap, FullTunesResolveIntegerAndConjugateBranch) {
    auto s     = settings();
    s.dt       = 2 * std::acos(-1.0) / 512;
    s.maxSteps = 1024;
    const OneTurnMap::Coordinates zero{}, h{1e-5, 1e-5, 1e-5, 1e-5};
    for (const auto tunes : {std::array<double, 2>{3.41, 2.72}, {1.41, .72}}) {
        OneTurnMap map(harmonicRing(tunes[0], tunes[1]), CoordinateSystemTrafo(), s);
        const auto matrix = map.jacobian(zero, h);
        const auto full   = map.fullTunes(zero, h, matrix);
        ASSERT_TRUE(full.available) << full.diagnostic;
        EXPECT_NEAR(full.tunes[0], tunes[0], 1e-8);
        EXPECT_NEAR(full.tunes[1], tunes[1], 1e-8);
        EXPECT_GT(full.samples, 500u);
        auto coupled  = matrix;
        coupled[0][2] = .01;
        EXPECT_FALSE(map.fullTunes(zero, h, coupled).available);
        auto unstable  = matrix;
        unstable[0][0] = 10;
        EXPECT_FALSE(map.fullTunes(zero, h, unstable).available);
    }
}

TEST(OneTurnMap, FullTunesRejectCoarseSamplingAndNonPlanarReference) {
    auto s     = settings();
    s.dt       = 2 * std::acos(-1.0) / 16;
    s.maxSteps = 64;
    const OneTurnMap::Coordinates zero{}, h{1e-5, 1e-5, 1e-5, 1e-5};
    OneTurnMap map(harmonicRing(5.41, 3.72), CoordinateSystemTrafo(), s);
    const auto matrix = map.jacobian(zero, h);
    EXPECT_FALSE(map.fullTunes(zero, h, matrix).available);
    auto offplane = zero;
    offplane[2]   = 1e-3;
    EXPECT_FALSE(map.fullTunes(offplane, h, matrix).available);
}

TEST(OneTurnMap, FullTunesRejectIntermediateCouplingEvenWhenReturnMapDecouples) {
    auto s          = settings();
    s.dt            = 2 * std::acos(-1.0) / 512;
    s.maxSteps      = 1024;
    const auto base = harmonicRing(3.41, 2.72);
    auto shear      = [](OneTurnMap::State ray, double sign) {
        const double t = ray.time;
        const Vector_t<double, 3> normal(std::cos(t), 0, std::sin(t));
        const Vector_t<double, 3> origin(std::cos(t) - 1, 0, std::sin(t));
        const Vector_t<double, 3> displacement = ray.position - origin;
        ray.position(1) += sign * std::sin(t) * dot(displacement, normal);
        return ray;
    };
    OneTurnMap map(
            [base, shear](const auto& ray, double dt, std::vector<OneTurnMap::Step>* accepted) {
                const auto start = shear(ray, -1);
                const auto end   = shear(base(start, dt, nullptr), 1);
                if (accepted) {
                    OneTurnMap::Step interval;
                    interval.midpoint = shear(base(start, dt / 2, nullptr), 1);
                    interval.end      = end;
                    interval.duration = dt;
                    accepted->push_back(interval);
                }
                return end;
            },
            CoordinateSystemTrafo(), s);
    const OneTurnMap::Coordinates zero{}, h{1e-5, 1e-5, 1e-5, 1e-5};
    const auto matrix = map.jacobian(zero, h);
    EXPECT_NEAR(matrix[2][0], 0, 1e-8);
    const auto full = map.fullTunes(zero, h, matrix);
    EXPECT_FALSE(full.available);
    EXPECT_NE(full.diagnostic.find("coupled transport"), std::string::npos) << full.diagnostic;
}
