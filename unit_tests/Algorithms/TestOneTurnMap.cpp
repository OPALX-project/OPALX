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
