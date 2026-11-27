/**
 * @file TestSpaceChargeMap.cpp
 * @brief Independent 6D event/force, canonical, covariance and TPS composition regressions.
 * All fixtures are constructed in C++; no OPAL input files or particles are run.
 */
#include <gtest/gtest.h>
#include <array>
#include <cmath>
#include <limits>

#include "Algorithms/MatchedCovariance.h"
#include "Algorithms/SpaceChargeMap.h"
#include "Ippl.h"
#include "Physics/Physics.h"
#include "SpaceCharge/CartesianPIC3D/RelativisticFieldComposer.h"
#include "Utilities/OpalException.h"

namespace {
    using namespace opalx::maps;
    using Chart              = SpaceChargeChart;
    using Matrix             = SpaceChargeMatrix;
    using Point              = std::array<double, 6>;
    using Vector             = opalx::spacecharge::analytic3d::Vector;
    using SpatialMatrix      = opalx::spacecharge::analytic3d::Matrix;
    using Profile            = opalx::spacecharge::analytic3d::Profile;
    constexpr Chart charts[] = {
            Chart::SectionSlopes, Chart::CommonTimeMechanical, Chart::CommonTimeCanonical};
    constexpr Profile profiles[] = {Profile::UniformEllipsoid, Profile::GaussianRms};

    Matrix identity() {
        Matrix result{};
        for (unsigned i = 0; i < 6; ++i)
            result[i][i] = 1;
        return result;
    }
    Matrix transpose(const Matrix& matrix) {
        Matrix result{};
        for (unsigned i = 0; i < 6; ++i)
            for (unsigned j = 0; j < 6; ++j)
                result[i][j] = matrix[j][i];
        return result;
    }
    // Independent small dense arithmetic, not the production TPS/Kokkos products.
    Matrix multiply(const Matrix& a, const Matrix& b) {
        Matrix result{};
        for (unsigned i = 0; i < 6; ++i)
            for (unsigned j = 0; j < 6; ++j)
                for (unsigned k = 0; k < 6; ++k)
                    result[i][j] += a[i][k] * b[k][j];
        return result;
    }
    Point apply(const Matrix& a, const Point& x) {
        Point result{};
        for (unsigned i = 0; i < 6; ++i)
            for (unsigned j = 0; j < 6; ++j)
                result[i] += a[i][j] * x[j];
        return result;
    }
    Matrix jacobian(const tps::Map<double>& map) {
        const auto j = map.jacobian(Point{});
        Matrix result{};
        for (unsigned i = 0; i < 6; ++i)
            for (unsigned k = 0; k < 6; ++k)
                result[i][k] = j(i, k);
        return result;
    }
    void expectMatrix(const Matrix& actual, const Matrix& expected, double tolerance = 2e-12) {
        for (unsigned i = 0; i < 6; ++i)
            for (unsigned j = 0; j < 6; ++j)
                EXPECT_NEAR(
                        actual[i][j], expected[i][j], tolerance * (1 + std::abs(expected[i][j])))
                        << i << ',' << j;
    }
    Matrix canonicalForm() {
        Matrix result{};
        for (unsigned i = 0; i < 6; i += 2) {
            result[i][i + 1] = 1;
            result[i + 1][i] = -1;
        }
        return result;
    }
    SpatialMatrix rotation() { return {{{0.36, 0.48, -0.8}, {-0.8, 0.6, 0}, {0.48, 0.64, 0.6}}}; }
    SpaceChargeReference reference(bool rotated = true, bool fields = true) {
        SpaceChargeReference ref;
        ref.massEV               = 8e8;
        ref.chargeE              = 1;
        auto& orbit              = ref.orbit;
        orbit.position           = {0.13, -0.22, 0.41};
        orbit.positionCorrection = {1e-17, -2e-17, 3e-17};
        orbit.time               = 2.3e-9;
        orbit.timeCorrection     = 1e-25;
        orbit.pathLength         = 1.25;
        const SpatialMatrix axes =
                rotated ? rotation() : SpatialMatrix{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
        const Vector u{0.4, -0.3, 2}, e{2e5, -3e5, 1e5}, b{0.3, -0.4, 0.7};
        for (unsigned j = 0; j < 3; ++j) {
            orbit.xAxis(j) = axes[0][j];
            orbit.yAxis(j) = axes[1][j];
            orbit.sAxis(j) = axes[2][j];
            for (unsigned i = 0; i < 3; ++i) {
                orbit.momentum(j) += axes[i][j] * u[i];
                if (fields) {
                    ref.electricFieldLab[j] += axes[i][j] * e[i];
                    ref.magneticFieldLab[j] += axes[i][j] * b[i];
                }
            }
        }
        return ref;
    }
    Matrix covariance() {
        const Point width{0.002, 0.0002, 0.003, 0.0003, 0.004, 0.0004};
        Matrix l{};
        for (unsigned i = 0; i < 6; ++i)
            for (unsigned j = 0; j <= i; ++j)
                l[i][j] = width[i] * (i == j ? 1 : 0.07 * (1 + (i + j) % 3));
        return multiply(l, transpose(l));
    }
    Vector crossProduct(const Vector& a, const Vector& b) {
        return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
    }
    Vector sectionVector(const SpaceChargeReference& ref, const Vector& lab) {
        Vector result{};
        for (unsigned j = 0; j < 3; ++j) {
            result[0] += ref.orbit.xAxis(j) * lab[j];
            result[1] += ref.orbit.yAxis(j) * lab[j];
            result[2] += ref.orbit.sAxis(j) * lab[j];
        }
        return result;
    }
    Vector sectionMomentum(const SpaceChargeReference& ref) {
        return sectionVector(
                ref, {ref.orbit.momentum(0), ref.orbit.momentum(1), ref.orbit.momentum(2)});
    }

    // Independent nonlinear event construction and RK4 Lorentz flow back to time t0.
    // This deliberately does not call the production section-intersection Jacobian.
    Point commonTimeEvent(const SpaceChargeReference& ref, const Point& x, bool canonical) {
        const auto u0      = sectionMomentum(ref);
        const double total = std::hypot(u0[0], u0[1], u0[2]);
        const double beta  = total / std::hypot(1.0, total);
        const double a = u0[0] / u0[2] + x[1], b = u0[1] / u0[2] + x[3];
        const double us = total * (1 + x[5]) / std::sqrt(1 + a * a + b * b);
        Point state{x[0], a * us, x[2], b * us, 0, us};
        const auto e        = sectionVector(ref, ref.electricFieldLab);
        const auto magnetic = sectionVector(ref, ref.magneticFieldLab);
        const double kappa  = ref.chargeE * Physics::c / ref.massEV;
        auto derivative     = [&](const Point& p) {
            const double gamma = std::sqrt(1 + p[1] * p[1] + p[3] * p[3] + p[5] * p[5]);
            const Vector velocity{
                    Physics::c * p[1] / gamma, Physics::c * p[3] / gamma,
                    Physics::c * p[5] / gamma};
            const auto magneticForce = crossProduct(velocity, magnetic);
            Point result{};
            for (unsigned i = 0; i < 3; ++i) {
                result[2 * i]     = velocity[i];
                result[2 * i + 1] = kappa * (e[i] + magneticForce[i]);
            }
            return result;
        };
        const double dt = x[4] / (beta * Physics::c);
        const auto k1   = derivative(state);
        auto displaced  = [&](const Point& k, double fraction) {
            Point result = state;
            for (unsigned i = 0; i < 6; ++i)
                result[i] += fraction * dt * k[i];
            return result;
        };
        const auto k2 = derivative(displaced(k1, 0.5));
        const auto k3 = derivative(displaced(k2, 0.5));
        const auto k4 = derivative(displaced(k3, 1));
        for (unsigned i = 0; i < 6; ++i)
            state[i] += dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
        for (unsigned i = 0; i < 3; ++i)
            state[2 * i + 1] -= u0[i];
        if (canonical) {
            const auto gauge = crossProduct(magnetic, {state[0], state[2], state[4]});
            for (unsigned i = 0; i < 3; ++i)
                state[2 * i + 1] += 0.5 * kappa * gauge[i];
        }
        return state;
    }
    Matrix eventJacobian(const SpaceChargeReference& ref, bool canonical, double step) {
        Matrix result{};
        for (unsigned j = 0; j < 6; ++j) {
            Point plus{}, minus{};
            plus[j]      = step;
            minus[j]     = -step;
            const auto p = commonTimeEvent(ref, plus, canonical);
            const auto m = commonTimeEvent(ref, minus, canonical);
            for (unsigned i = 0; i < 6; ++i)
                result[i][j] = (p[i] - m[i]) / (2 * step);
        }
        return result;
    }
    double difference(const Matrix& a, const Matrix& b) {
        double result = 0;
        for (unsigned i = 0; i < 6; ++i)
            for (unsigned j = 0; j < 6; ++j)
                result = std::hypot(result, a[i][j] - b[i][j]);
        return result;
    }
    double trace(const Matrix& matrix) {
        double result = 0;
        for (unsigned i = 0; i < 6; ++i)
            result += matrix[i][i];
        return result;
    }
    class SpaceChargeMapTest : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
        }
        static void TearDownTestSuite() { ippl::finalize(); }
        void SetUp() override { order = tps::Tps<double>::getGlobalTruncOrder(); }
        void TearDown() override { tps::Tps<double>::setGlobalTruncOrder(order); }
        int order;
    };
}  // namespace

TEST_F(SpaceChargeMapTest, AlignedFieldFreeScalingAndArrivalTimeSign) {
    auto ref           = reference(false, false);
    ref.orbit.momentum = {0, 0, 2};
    const SpaceChargeCoordinates coordinates(ref);
    auto expected = identity();
    for (unsigned i : {1u, 3u, 5u})
        expected[i][i] = 2;
    expectMatrix(coordinates.jacobian(Chart::SectionSlopes, Chart::CommonTimeMechanical), expected);
    expectMatrix(coordinates.jacobian(Chart::SectionSlopes, Chart::CommonTimeCanonical), expected);
    const auto map = coordinates.map(Chart::SectionSlopes, Chart::CommonTimeMechanical, "arrival");
    const auto output = map.evaluate(Point{0, 0, 0, 0, 0.001, 0});
    EXPECT_DOUBLE_EQ(output[4], 0.001);  // Early section arrival is ahead at the common time.
    EXPECT_EQ(coordinates.chart(Chart::SectionSlopes), linearTransferMapChart(ref.orbit));
    EXPECT_EQ(
            coordinates.chart(Chart::CommonTimeMechanical).referenceTime,
            ref.orbit.time - ref.orbit.timeCorrection);
}

TEST_F(SpaceChargeMapTest, IndependentNonlinearEventsAndLorentzFlowMatchFullJacobian) {
    for (bool rotated : {false, true}) {
        const auto ref = reference(rotated, true);
        const SpaceChargeCoordinates coordinates(ref);
        for (bool canonical : {false, true}) {
            const auto kind  = canonical ? Chart::CommonTimeCanonical : Chart::CommonTimeMechanical;
            const auto exact = coordinates.jacobian(Chart::SectionSlopes, kind);
            expectMatrix(eventJacobian(ref, canonical, 1e-5), exact, 2e-9);
            const double coarse = difference(eventJacobian(ref, canonical, 1e-3), exact);
            const double fine   = difference(eventJacobian(ref, canonical, 5e-4), exact);
            EXPECT_GT(coarse / fine, 3.9);
            EXPECT_LT(coarse / fine, 4.1);
        }
    }
}

TEST_F(SpaceChargeMapTest, EveryChartPairHasInverseAndPreservesCovariance) {
    const SpaceChargeCoordinates coordinates(reference());
    const auto initial = covariance();
    for (auto from : charts)
        for (auto to : charts) {
            const auto forward = coordinates.jacobian(from, to);
            const auto inverse = coordinates.jacobian(to, from);
            expectMatrix(multiply(inverse, forward), identity());
            const auto converted = coordinates.convertCovariance(initial, from, to);
            const auto expected  = multiply(multiply(forward, initial), transpose(forward));
            const auto restored  = coordinates.convertCovariance(converted, to, from);
            for (unsigned i = 0; i < 6; ++i)
                for (unsigned j = 0; j < 6; ++j) {
                    const double scale = std::sqrt(initial[i][i] * initial[j][j]);
                    EXPECT_NEAR(restored[i][j], initial[i][j], 2e-12 * scale);
                    EXPECT_NEAR(
                            converted[i][j], expected[i][j],
                            2e-12 * std::sqrt(converted[i][i] * converted[j][j]));
                }
            const auto roundTrip = tps::compose(
                    coordinates.map(to, from, "inverse"), coordinates.map(from, to, "forward"), 1);
            expectMatrix(jacobian(roundTrip), identity());
        }
}

TEST_F(SpaceChargeMapTest, FullMagneticCanonicalGaugeAndExactMetadataCompatibility) {
    auto ref = reference(false, true);
    const SpaceChargeCoordinates coordinates(ref);
    const auto c = coordinates.jacobian(Chart::CommonTimeMechanical, Chart::CommonTimeCanonical);
    const double kappa = ref.chargeE * Physics::c / ref.massEV;
    const Point p{0.001, 0.002, -0.003, 0.004, 0.005, -0.006};
    const auto converted = apply(c, p);
    const auto gauge     = crossProduct(ref.magneticFieldLab, {p[0], p[2], p[4]});
    for (unsigned i = 0; i < 3; ++i)
        EXPECT_NEAR(converted[2 * i + 1], p[2 * i + 1] + 0.5 * kappa * gauge[i], 1e-15);
    const auto omega = multiply(multiply(transpose(c), canonicalForm()), c);
    const auto transverse =
            MatchedCovariance::canonicalTransform(ref.magneticFieldLab[2], ref.chargeE, ref.massEV);
    for (unsigned i = 0; i < 4; ++i)
        for (unsigned j = 0; j < 4; ++j)
            EXPECT_DOUBLE_EQ(c[i][j], transverse[i][j]);
    EXPECT_GT(difference(omega, canonicalForm()), 0.1);
    EXPECT_EQ(coordinates.chart(Chart::CommonTimeCanonical).canonicalPairs.size(), 3u);
    EXPECT_EQ(coordinates.chart(Chart::CommonTimeCanonical).parameters.size(), 4u);
    ref.magneticFieldLab[0] += 0.1;
    const SpaceChargeCoordinates changedGauge(ref);
    EXPECT_EQ(
            coordinates.chart(Chart::CommonTimeMechanical),
            changedGauge.chart(Chart::CommonTimeMechanical));
    EXPECT_NE(
            coordinates.chart(Chart::CommonTimeCanonical),
            changedGauge.chart(Chart::CommonTimeCanonical));
    auto fieldFree = reference(false, false);
    const SpaceChargeCoordinates firstSpecies(fieldFree);
    fieldFree.massEV *= 2;
    const SpaceChargeCoordinates secondSpecies(fieldFree);
    for (auto kind : {Chart::CommonTimeMechanical, Chart::CommonTimeCanonical}) {
        EXPECT_NE(firstSpecies.chart(kind), secondSpecies.chart(kind));
        EXPECT_THROW(
                (void)tps::compose(
                        secondSpecies.map(kind, Chart::SectionSlopes, "other-species"),
                        firstSpecies.map(Chart::SectionSlopes, kind, "species")),
                std::invalid_argument);
    }
    EXPECT_THROW(
            (void)tps::compose(
                    changedGauge.map(
                            Chart::CommonTimeCanonical, Chart::CommonTimeMechanical, "out"),
                    coordinates.map(Chart::CommonTimeMechanical, Chart::CommonTimeCanonical, "in")),
            std::invalid_argument);
}

TEST_F(SpaceChargeMapTest, KickDensityUsesFullSectionCovarianceAndReferenceDirection) {
    const SpaceChargeCoordinates coordinates(reference());
    const auto sigma    = covariance();
    const auto t        = coordinates.jacobian(Chart::SectionSlopes, Chart::CommonTimeMechanical);
    const auto sameTime = multiply(multiply(t, sigma), transpose(t));
    SpatialMatrix d{};
    for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j)
            d[i][j] = (i == j ? 1 : 0)
                      + (coordinates.gamma() - 1) * coordinates.direction()[i]
                                * coordinates.direction()[j];
    const auto kick = makeSpaceChargeKick(
            coordinates, sigma, Chart::SectionSlopes, Profile::GaussianRms, 1e-10, 0.03, "sc");
    for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j) {
            double expected = 0;
            for (unsigned a = 0; a < 3; ++a)
                for (unsigned b = 0; b < 3; ++b)
                    expected += d[i][a] * sameTime[2 * a][2 * b] * d[j][b];
            EXPECT_NEAR(
                    kick.restSpatialCovariance[i][j], expected,
                    2e-12
                            * std::sqrt(
                                    kick.restSpatialCovariance[i][i]
                                    * kick.restSpatialCovariance[j][j]));
        }
    // The section x-z correlations and transverse reference slopes contribute to density.
    EXPECT_GT(std::abs(sameTime[0][4]), 1e-7);
    const auto mechanicalSigma =
            coordinates.convertCovariance(sigma, Chart::SectionSlopes, Chart::CommonTimeMechanical);
    const auto canonicalSigma =
            coordinates.convertCovariance(sigma, Chart::SectionSlopes, Chart::CommonTimeCanonical);
    const auto canonicalKick = makeSpaceChargeKick(
            coordinates, canonicalSigma, Chart::CommonTimeCanonical, Profile::GaussianRms, 1e-10,
            0.03, "can");
    const auto mechanicalKick = makeSpaceChargeKick(
            coordinates, mechanicalSigma, Chart::CommonTimeMechanical, Profile::GaussianRms, 1e-10,
            0.03, "mech");
    expectMatrix(jacobian(canonicalKick.polynomial), jacobian(mechanicalKick.polynomial));
    const auto equivalent = tps::compose(
            coordinates.map(Chart::CommonTimeMechanical, Chart::SectionSlopes, "back"),
            tps::compose(
                    mechanicalKick.polynomial,
                    coordinates.map(Chart::SectionSlopes, Chart::CommonTimeMechanical, "toTime")));
    expectMatrix(jacobian(kick.polynomial), jacobian(equivalent));
}

TEST_F(SpaceChargeMapTest, CanonicalKickIsSymplecticAndPreservesAllThreeMomentInvariants) {
    const SpaceChargeCoordinates coordinates(reference());
    const auto sigma = covariance();
    for (auto profile : profiles) {
        const auto kick = makeSpaceChargeKick(
                coordinates, sigma, Chart::CommonTimeCanonical, profile, 1e-10, 0.03, "sc");
        const auto k = jacobian(kick.polynomial);
        expectMatrix(multiply(multiply(transpose(k), canonicalForm()), k), canonicalForm());
        const auto transported = multiply(multiply(k, sigma), transpose(k));
        const auto originalSquare =
                multiply(multiply(canonicalForm(), sigma), multiply(canonicalForm(), sigma));
        const auto transportedSquare = multiply(
                multiply(canonicalForm(), transported), multiply(canonicalForm(), transported));
        auto before = identity(), after = identity();
        for (unsigned power = 1; power <= 3; ++power) {
            before = multiply(before, originalSquare);
            after  = multiply(after, transportedSquare);
            EXPECT_NEAR(trace(after), trace(before), 2e-11 * std::abs(trace(before)));
        }
        for (unsigned i = 0; i < 3; ++i)
            for (unsigned j = 0; j < 3; ++j) {
                EXPECT_DOUBLE_EQ(kick.forceGradient[i][j], kick.forceGradient[j][i]);
                EXPECT_GE(kick.forceAbsoluteError[i][j], 0);
            }
    }
}

TEST_F(SpaceChargeMapTest, KicksAgreeWithExistingLorentzFieldsAndIndependentSiTimeImpulse) {
    using opalx::spacecharge::FieldCompositionPolicy;
    using opalx::spacecharge::LorentzBoost;
    for (double totalMomentum : {0.01, 2.0, 50.0}) {
        auto ref = reference(false, false);
        const Vector n{0.36, 0.48, 0.8};
        ref.orbit.momentum = {totalMomentum * n[0], totalMomentum * n[1], totalMomentum * n[2]};
        const SpaceChargeCoordinates coordinates(ref);
        for (auto profile : profiles) {
            const double h  = 0.03;
            const auto kick = makeSpaceChargeKick(
                    coordinates, covariance(), Chart::CommonTimeMechanical, profile, 1e-10, h,
                    "sc");
            const Point point{2e-4, 1e-5, -3e-4, -2e-5, 1e-4, 3e-5};
            const auto output = kick.polynomial.evaluate(point);
            LorentzBoost::Vector restPosition(0), electricRest(0), electricLab, magneticLab;
            const double dot = n[0] * point[0] + n[1] * point[2] + n[2] * point[4];
            for (unsigned i = 0; i < 3; ++i)
                restPosition[i] = point[2 * i] + (coordinates.gamma() - 1) * n[i] * dot;
            for (unsigned i = 0; i < 3; ++i)
                for (unsigned j = 0; j < 3; ++j)
                    electricRest[i] += kick.restTensor.gradient[i][j] * restPosition[j];
            FieldCompositionPolicy policy;
            policy.gamma        = coordinates.gamma();
            policy.meanMomentum = {
                    ref.orbit.momentum(0), ref.orbit.momentum(1), ref.orbit.momentum(2)};
            const auto boost = LorentzBoost::from(policy);
            boost.toLab(electricRest, electricLab, magneticLab);
            const auto forcePerCharge = electricLab + cross(boost.velocity, magneticLab);
            const double massKg       = ref.massEV * Physics::q_e / (Physics::c * Physics::c);
            const double dt           = h / (coordinates.beta() * Physics::c);
            const double impulse      = dt * ref.chargeE * Physics::q_e / (massKg * Physics::c);
            for (unsigned i = 0; i < 3; ++i) {
                EXPECT_DOUBLE_EQ(output[2 * i], point[2 * i]);
                const double expected = point[2 * i + 1] + impulse * forcePerCharge[i];
                EXPECT_NEAR(
                        output[2 * i + 1], expected,
                        2e-10 * std::abs(impulse * forcePerCharge[i]) + 1e-16);
            }
        }
    }
}

TEST_F(SpaceChargeMapTest, SignedChargePathAndRestEnergyScaling) {
    auto ref = reference(false, false);
    const SpaceChargeCoordinates coordinates(ref);
    for (auto profile : profiles) {
        const auto original =
                jacobian(makeSpaceChargeKick(
                                 coordinates, covariance(), Chart::CommonTimeCanonical, profile,
                                 1e-10, 0.03, "a")
                                 .polynomial);
        const auto negative =
                jacobian(makeSpaceChargeKick(
                                 coordinates, covariance(), Chart::CommonTimeCanonical, profile,
                                 -2e-10, 0.03, "b")
                                 .polynomial);
        const auto longer = jacobian(makeSpaceChargeKick(
                                             coordinates, covariance(), Chart::CommonTimeCanonical,
                                             profile, 1e-10, 0.06, "c")
                                             .polynomial);
        ref.chargeE       = -1;
        const SpaceChargeCoordinates negativeSpecies(ref);
        const auto likeCharge =
                jacobian(makeSpaceChargeKick(
                                 negativeSpecies, covariance(), Chart::CommonTimeCanonical, profile,
                                 -1e-10, 0.03, "d")
                                 .polynomial);
        ref.chargeE = 1;
        ref.massEV *= 2;
        const SpaceChargeCoordinates heavierSpecies(ref);
        const auto heavier = jacobian(makeSpaceChargeKick(
                                              heavierSpecies, covariance(),
                                              Chart::CommonTimeCanonical, profile, 1e-10, 0.03, "e")
                                              .polynomial);
        ref.massEV /= 2;
        for (unsigned i = 0; i < 3; ++i)
            for (unsigned j = 0; j < 3; ++j) {
                const double coefficient = original[2 * i + 1][2 * j];
                EXPECT_NEAR(
                        negative[2 * i + 1][2 * j], -2 * coefficient,
                        2e-12 * std::abs(coefficient));
                EXPECT_NEAR(
                        longer[2 * i + 1][2 * j], 2 * coefficient, 2e-12 * std::abs(coefficient));
                EXPECT_NEAR(
                        heavier[2 * i + 1][2 * j], coefficient / 2, 2e-12 * std::abs(coefficient));
            }
        expectMatrix(likeCharge, original);
    }
}

TEST_F(SpaceChargeMapTest, ExtremePathAndChargeKeepRepresentableImpulse) {
    auto ref           = reference(false, false);
    ref.orbit.momentum = {0, 0, 0.01};
    const SpaceChargeCoordinates coordinates(ref);
    const auto huge = makeSpaceChargeKick(
            coordinates, covariance(), Chart::CommonTimeCanonical, Profile::GaussianRms, 1e-300,
            1e308, "huge");
    const auto equivalent = makeSpaceChargeKick(
            coordinates, covariance(), Chart::CommonTimeCanonical, Profile::GaussianRms, 1e8, 1,
            "equivalent");
    expectMatrix(jacobian(huge.polynomial), jacobian(equivalent.polynomial), 2e-11);
    EXPECT_THROW(
            (void)makeSpaceChargeKick(
                    coordinates, covariance(), Chart::CommonTimeCanonical, Profile::GaussianRms, 1,
                    1e308, "overflow"),
            OpalException);
    EXPECT_THROW(
            (void)makeSpaceChargeKick(
                    coordinates, covariance(), Chart::CommonTimeCanonical, Profile::GaussianRms,
                    1e-300, 1e-300, "underflow"),
            OpalException);
}

TEST_F(SpaceChargeMapTest, NearlyTangentSlopeDerivativesRetainSmallPositiveTerms) {
    auto ref = reference(false, false);
    for (bool alongX : {true, false}) {
        ref.orbit.momentum = {alongX ? 1.0 : 0.0, alongX ? 0.0 : 1.0, 1e-9};
        const SpaceChargeCoordinates coordinates(ref);
        const auto matrix = coordinates.jacobian(Chart::SectionSlopes, Chart::CommonTimeMechanical);
        const auto component = alongX ? 1u : 3u;
        EXPECT_NEAR(matrix[component][component], 1e-27, 1e-40);
        EXPECT_LT(matrix[5][component], 0);
    }
    ref.orbit.momentum = {1, 0, 0.5e-12};
    EXPECT_THROW((void)SpaceChargeCoordinates(ref), OpalException);
}

TEST_F(SpaceChargeMapTest, ReversibleKicksAndAmbientOrderRestorationInEveryChart) {
    const SpaceChargeCoordinates coordinates(reference());
    for (int ambient : {0, 1, 4, tps::Tps<double>::EXACT}) {
        tps::Tps<double>::setGlobalTruncOrder(ambient);
        for (auto kind : charts) {
            const auto forward = makeSpaceChargeKick(
                    coordinates, covariance(), kind, Profile::GaussianRms, 1e-10, 0.03, "forward");
            const auto backward = makeSpaceChargeKick(
                    coordinates, covariance(), kind, Profile::GaussianRms, 1e-10, -0.03,
                    "backward");
            expectMatrix(
                    jacobian(tps::compose(backward.polynomial, forward.polynomial, 1)), identity());
            EXPECT_EQ(tps::Tps<double>::getGlobalTruncOrder(), ambient);
            for (unsigned i = 0; i < 6; ++i)
                EXPECT_DOUBLE_EQ(forward.polynomial.polynomial()[i][0], 0);
        }
    }
}

TEST_F(SpaceChargeMapTest, ExternalSegmentAndKickComposeWithExactSharedCharts) {
    auto entrance = reference(false, false), exit = entrance;
    entrance.orbit.momentum = {0, 0, 2};
    exit.orbit.momentum     = entrance.orbit.momentum;
    const double length = 0.4, gamma = std::sqrt(5.0), beta = 2 / gamma;
    exit.orbit.position(2) += length;
    exit.orbit.time += length / (beta * Physics::c);
    exit.orbit.pathLength += length;
    const SpaceChargeCoordinates in(entrance), out(exit);
    LinearTransferMapDraft draft;
    draft.entrance = entrance.orbit;
    draft.exit     = exit.orbit;
    for (unsigned i = 0; i < 6; ++i)
        draft.matrix(i, i) = 1;
    draft.matrix(0, 1) = draft.matrix(2, 3) = length;
    draft.matrix(4, 5)                      = length / (gamma * gamma);
    const auto external                     = toTpsLinearTransferMap(draft, "drift");
    const auto canonical                    = tps::compose(
            out.map(Chart::SectionSlopes, Chart::CommonTimeCanonical, "exit-chart"),
            tps::compose(
                    external,
                    in.map(Chart::CommonTimeCanonical, Chart::SectionSlopes, "entrance-chart")));
    auto expected  = identity();
    expected[0][1] = expected[2][3] = length / 2;
    expected[4][5]                  = length / (2 * gamma * gamma);
    expectMatrix(jacobian(canonical), expected);
    const auto sigma = multiply(multiply(expected, covariance()), transpose(expected));
    const auto kick  = makeSpaceChargeKick(
            out, sigma, Chart::CommonTimeCanonical, Profile::GaussianRms, 1e-10, 0.03, "sc");
    const auto composed = tps::compose(kick.polynomial, canonical);
    expectMatrix(jacobian(composed), multiply(jacobian(kick.polynomial), expected));
    EXPECT_GT(
            difference(
                    multiply(jacobian(kick.polynomial), expected),
                    multiply(expected, jacobian(kick.polynomial))),
            1e-5);
    expectMatrix(
            multiply(multiply(transpose(jacobian(composed)), canonicalForm()), jacobian(composed)),
            canonicalForm());
}

TEST_F(SpaceChargeMapTest, ZeroChargePathAndNeutralSpeciesGiveExactIdentity) {
    auto ref = reference();
    for (unsigned mode = 0; mode < 3; ++mode) {
        if (mode == 2) ref.chargeE = 0;
        const SpaceChargeCoordinates coordinates(ref);
        for (auto kind : charts) {
            const auto kick = makeSpaceChargeKick(
                    coordinates, covariance(), kind, Profile::GaussianRms, mode == 0 ? 0 : 1e-10,
                    mode == 1 ? 0 : 0.03, "zero");
            const auto actual   = jacobian(kick.polynomial);
            const auto expected = identity();
            for (unsigned i = 0; i < 6; ++i)
                for (unsigned j = 0; j < 6; ++j)
                    EXPECT_DOUBLE_EQ(actual[i][j], expected[i][j]);
        }
    }
}

TEST_F(SpaceChargeMapTest, RejectsInvalidReferenceCovarianceAndKickInputs) {
    const double nan      = std::numeric_limits<double>::quiet_NaN();
    const double infinity = std::numeric_limits<double>::infinity();
    for (unsigned mode = 0; mode < 10; ++mode) {
        auto ref = reference();
        if (mode == 0) ref.massEV = 0;
        if (mode == 1) ref.massEV = nan;
        if (mode == 2) ref.chargeE = infinity;
        if (mode == 3) ref.electricFieldLab[0] = nan;
        if (mode == 4) ref.magneticFieldLab[2] = infinity;
        if (mode == 5) ref.orbit.momentum = {0, 0, 0};
        if (mode == 6) ref.orbit.momentum = -ref.orbit.sAxis;
        if (mode == 7) ref.orbit.xAxis *= 2;
        if (mode == 8) ref.orbit.timeCorrection = infinity;
        if (mode == 9) ref.orbit.momentum = ref.orbit.xAxis;
        EXPECT_THROW((void)SpaceChargeCoordinates(ref), OpalException) << mode;
    }
    const SpaceChargeCoordinates coordinates(reference());
    for (unsigned mode = 0; mode < 5; ++mode) {
        auto sigma = covariance();
        if (mode == 0) sigma[0][0] = 0;
        if (mode == 1) sigma[5][5] = nan;
        if (mode == 2) sigma[0][2] = infinity;
        if (mode == 3) sigma[2][0] += 1e-4;
        if (mode == 4) sigma[0][2] = sigma[2][0] = 1;
        EXPECT_THROW(
                (void)coordinates.convertCovariance(
                        sigma, Chart::SectionSlopes, Chart::CommonTimeCanonical),
                OpalException);
        EXPECT_THROW(
                (void)makeSpaceChargeKick(
                        coordinates, sigma, Chart::SectionSlopes, Profile::GaussianRms, 0, 0,
                        "invalid"),
                OpalException);
    }
    EXPECT_THROW((void)coordinates.chart(static_cast<Chart>(99)), OpalException);
    EXPECT_THROW(
            (void)coordinates.jacobian(static_cast<Chart>(99), Chart::SectionSlopes),
            OpalException);
    EXPECT_THROW(
            (void)coordinates.map(Chart::SectionSlopes, Chart::CommonTimeCanonical, ""),
            OpalException);
    for (double path : {nan, infinity})
        EXPECT_THROW(
                (void)makeSpaceChargeKick(
                        coordinates, covariance(), Chart::CommonTimeMechanical,
                        Profile::GaussianRms, 1e-10, path, "invalid"),
                OpalException);
    EXPECT_THROW(
            (void)makeSpaceChargeKick(
                    coordinates, covariance(), Chart::CommonTimeMechanical, Profile::GaussianRms,
                    nan, 0.03, "invalid"),
            OpalException);
    EXPECT_THROW(
            (void)makeSpaceChargeKick(
                    coordinates, covariance(), Chart::CommonTimeMechanical,
                    static_cast<Profile>(99), 1e-10, 0.03, "invalid"),
            OpalException);
    EXPECT_THROW(
            (void)makeSpaceChargeKick(
                    coordinates, covariance(), Chart::CommonTimeCanonical, Profile::GaussianRms,
                    1e-10, 0.03, ""),
            OpalException);
    opalx::spacecharge::analytic3d::IntegrationSettings exhausted;
    exhausted.maximumIntervals = 1;
    EXPECT_THROW(
            (void)makeSpaceChargeKick(
                    coordinates, covariance(), Chart::CommonTimeMechanical, Profile::GaussianRms,
                    1e-10, 0.03, "invalid", exhausted),
            OpalException);
    EXPECT_EQ(tps::Tps<double>::getGlobalTruncOrder(), order);
}

TEST_F(SpaceChargeMapTest, HostCoordinateAndKickResultsAgreeAcrossMpiRanks) {
    const SpaceChargeCoordinates coordinates(reference());
    const auto kick = makeSpaceChargeKick(
            coordinates, covariance(), Chart::CommonTimeCanonical, Profile::GaussianRms, 1e-10,
            0.03, "sc");
    const auto matrix = jacobian(kick.polynomial);
    std::array<double, 36> values{}, minimum{}, maximum{};
    for (unsigned i = 0; i < 6; ++i)
        for (unsigned j = 0; j < 6; ++j)
            values[6 * i + j] = matrix[i][j];
    MPI_Allreduce(values.data(), minimum.data(), 36, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(values.data(), maximum.data(), 36, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    EXPECT_EQ(minimum, maximum);
}
