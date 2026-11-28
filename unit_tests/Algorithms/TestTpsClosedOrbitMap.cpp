#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
#include "Algorithms/LinearMapEigenAnalysis.h"
#include "Algorithms/TpsClosedOrbitMap.h"
#include "Ippl.h"
#include "Utilities/OpalException.h"
#include "gtest/gtest.h"

namespace {
    using EV     = LinearMapEigenAnalysis;
    using Matrix = EV::Matrix;
    using Point  = OneTurnMap::Coordinates;
    using Series = tps::Tps<double>;

    Matrix rotations(double qx = 0.155, double qy = 0.89, double rx = 1, double ry = 1) {
        Matrix result{};
        for (unsigned plane = 0; plane < 2; ++plane) {
            const unsigned i    = 2 * plane;
            const double phase  = 2 * std::acos(-1.0) * (plane ? qy : qx);
            const double radius = plane ? ry : rx;
            result[i][i] = result[i + 1][i + 1] = radius * std::cos(phase);
            result[i][i + 1]                    = radius * std::sin(phase);
            result[i + 1][i]                    = -radius * std::sin(phase);
        }
        return result;
    }
    Matrix coupled(const Matrix& matrix) {
        const double c = std::cos(0.6), s = std::sin(0.6);
        const Matrix rotation{{{c, 0, s, 0}, {0, c, 0, s}, {-s, 0, c, 0}, {0, -s, 0, c}}};
        Matrix result{};
        for (unsigned i = 0; i < 4; ++i)
            for (unsigned j = 0; j < 4; ++j)
                for (unsigned k = 0; k < 4; ++k)
                    for (unsigned l = 0; l < 4; ++l)
                        result[i][j] += rotation[i][k] * matrix[k][l] * rotation[j][l];
        return result;
    }
    ClosedOrbitSolver::Result orbit(const Matrix& matrix = rotations()) {
        ClosedOrbitSolver::Result result;
        result.status      = ClosedOrbitSolver::Status::Converged;
        result.coordinates = {0.02, 0.03, -0.01, -0.04};
        result.residual    = {1e-11, 2e-12, -3e-11, -4e-12};
        result.matrix      = matrix;
        return result;
    }
    tps::Map<double> rechart(
            const tps::Map<double>& map, const tps::PhaseSpaceChart& input,
            const tps::PhaseSpaceChart& output, int order = 1) {
        return tps::Map<double>(tps::MapContext(input, output, order), map.polynomial());
    }
    void expectSpectrum(const EV::Result& actual, const EV::Result& expected) {
        EXPECT_EQ(actual.stability, expected.stability);
        EXPECT_EQ(actual.nearInteger, expected.nearInteger);
        EXPECT_NEAR(actual.maximumResidual, expected.maximumResidual, 1e-13);
        EXPECT_NEAR(actual.basisReciprocalCondition, expected.basisReciprocalCondition, 1e-13);
        std::array<bool, 4> used{};
        // Compare eigenvalue sets; LAPACK's ordering and vector phases are not a contract.
        for (unsigned expectedIndex = 0; expectedIndex < 4; ++expectedIndex) {
            const auto value = expected.eigenvalues[expectedIndex];
            unsigned best    = 4;
            double distance  = std::numeric_limits<double>::infinity();
            for (unsigned i = 0; i < 4; ++i) {
                const double candidate = std::abs(actual.eigenvalues[i] - value);
                if (!used[i] && candidate < distance) {
                    distance = candidate;
                    best     = i;
                }
            }
            ASSERT_LT(best, 4u);
            used[best] = true;
            EXPECT_LT(distance, 1e-13);
            // Stable fixtures have distinct modes. Check vector directions up to
            // arbitrary complex phase, catching a transposed Jacobian as well.
            if (expected.stability == EV::Stability::Stable) {
                EV::Complex overlap{};
                for (unsigned i = 0; i < 4; ++i)
                    overlap += std::conj(expected.eigenvectors[expectedIndex][i])
                               * actual.eigenvectors[best][i];
                EXPECT_NEAR(std::abs(overlap), 1, 1e-13);
            }
        }
        ASSERT_EQ(actual.modes.size(), expected.modes.size());
        for (unsigned i = 0; i < expected.modes.size(); ++i) {
            EXPECT_EQ(actual.modes[i].nearInteger, expected.modes[i].nearInteger);
            ASSERT_EQ(actual.modes[i].tune.has_value(), expected.modes[i].tune.has_value());
            ASSERT_EQ(
                    actual.modes[i].complementaryTune.has_value(),
                    expected.modes[i].complementaryTune.has_value());
            if (actual.modes[i].tune) {
                EXPECT_NEAR(*actual.modes[i].tune, *expected.modes[i].tune, 1e-13);
                EXPECT_NEAR(
                        *actual.modes[i].complementaryTune, *expected.modes[i].complementaryTune,
                        1e-13);
            }
        }
    }
    class TpsClosedOrbitMapTest : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
        }
        static void TearDownTestSuite() { ippl::finalize(); }
        void SetUp() override { previousOrder = Series::getGlobalTruncOrder(); }
        void TearDown() override { Series::setGlobalTruncOrder(previousOrder); }
        int previousOrder;
    };
}  // namespace

TEST_F(TpsClosedOrbitMapTest, KeepsSolvedReferenceResidualAndRotatedMechanicalChart) {
    const auto source = orbit(coupled(rotations()));
    const auto saved  = source;
    const CoordinateSystemTrafo section(
            Vector_t<double, 3>(3, -2, 1), Quaternion(std::cos(0.2), 0, std::sin(0.2), 0));
    for (int order : {0, 1, 4, Series::EXACT}) {
        Series::setGlobalTruncOrder(order);
        const auto map    = opalx::maps::toTpsClosedOrbitMap(source, section, 2.3);
        const auto& chart = map.context().input();
        EXPECT_EQ(chart, map.context().output());
        EXPECT_EQ(chart.convention, tps::CoordinateConvention::mechanical);
        EXPECT_TRUE(chart.canonicalPairs.empty());
        EXPECT_EQ(
                chart.reference,
                (std::vector<double>(source.coordinates.begin(), source.coordinates.end())));
        EXPECT_FALSE(chart.referenceTime);
        EXPECT_EQ(chart.referenceMomentum, 2.3);
        EXPECT_EQ(chart.momentumUnit, "p/(mc)");
        const std::array<std::string, 4> names{"x", "px", "y", "py"};
        for (unsigned i = 0; i < 4; ++i) {
            EXPECT_EQ(chart.coordinates[i].name, names[i]);
            EXPECT_EQ(chart.coordinates[i].unit, i % 2 ? "p/(mc)" : "m");
            EXPECT_EQ(chart.coordinates[i].scale, 1);
            EXPECT_EQ(map.polynomial()[i][0], source.residual[i]);
        }
        for (unsigned i = 0; i < 3; ++i) {
            EXPECT_EQ(chart.frame.origin[i], section.getOrigin()(i));
            for (unsigned j = 0; j < 3; ++j)
                EXPECT_EQ(chart.frame.axes[3 * i + j], section.getRotationMatrix()(i, j));
        }
        const auto matrix = map.jacobian(Point{});
        const Point delta{0.001, -0.002, 0.003, -0.004};
        const auto evaluated = map.evaluate(delta);
        const auto zero      = map.evaluate(Point{});
        for (unsigned i = 0; i < 4; ++i) {
            EXPECT_EQ(zero[i], source.residual[i]);
            double expected = source.residual[i];
            for (unsigned j = 0; j < 4; ++j) {
                EXPECT_EQ(matrix(i, j), source.matrix[i][j]);
                expected += source.matrix[i][j] * delta[j];
            }
            EXPECT_NEAR(evaluated[i], expected, 2e-13);
        }
        (void)EV::analyze(map, EV::Settings{});
        EXPECT_EQ(Series::getGlobalTruncOrder(), order);
    }
    EXPECT_EQ(source.coordinates, saved.coordinates);
    EXPECT_EQ(source.residual, saved.residual);
    EXPECT_EQ(source.matrix, saved.matrix);
}

TEST_F(TpsClosedOrbitMapTest, SpectrumMatchesExistingLapackAcrossStabilityCases) {
    Matrix identity{};
    for (unsigned i = 0; i < 4; ++i)
        identity[i][i] = 1;
    auto minusOne  = identity;
    minusOne[0][0] = minusOne[1][1] = -1;
    auto realUnstable               = rotations();
    realUnstable[0][0]              = 2;
    realUnstable[1][1]              = 0.5;
    realUnstable[0][1] = realUnstable[1][0] = 0;
    auto defective                          = rotations(0.2, 0.2);
    defective[0][2] = defective[1][3] = 1;
    for (const auto& matrix :
         {rotations(), coupled(rotations(0.155, 0.31)), rotations(1e-8, 0.3), identity, minusOne,
          realUnstable, coupled(rotations(0.2, 0.2, 1.1, 1 / 1.1)), rotations(0.2, 0.3, 0.9, 0.8),
          defective}) {
        const auto map =
                opalx::maps::toTpsClosedOrbitMap(orbit(matrix), CoordinateSystemTrafo(), 2.3);
        expectSpectrum(EV::analyze(map, EV::Settings{}), EV::analyze(matrix, EV::Settings{}));
    }
    const auto map    = opalx::maps::toTpsClosedOrbitMap(orbit(), CoordinateSystemTrafo(), 2.3);
    const auto result = EV::analyze(map, EV::Settings{});
    ASSERT_EQ(result.modes.size(), 2u);
    EXPECT_NEAR(*result.modes[0].tune, 0.11, 1e-13);
    EXPECT_NEAR(*result.modes[0].complementaryTune, 0.89, 1e-13);
    std::ostringstream oldReport, tpsReport;
    EV::writeReport(oldReport, EV::analyze(orbit().matrix, EV::Settings{}));
    EV::writeReport(tpsReport, result);
    EXPECT_EQ(tpsReport.str(), oldReport.str());
}

TEST_F(TpsClosedOrbitMapTest, ChartNormalizationDoesNotDoubleApplyEigenanalysisScales) {
    auto source = orbit(coupled(rotations(0.155, 0.31)));
    EV::Settings settings;
    settings.scales = {100, 0.001, 10, 0.01};
    for (unsigned i = 0; i < 4; ++i)
        for (unsigned j = 0; j < 4; ++j)
            source.matrix[i][j] *= settings.scales[i] / settings.scales[j];
    const auto unitMap = opalx::maps::toTpsClosedOrbitMap(source, CoordinateSystemTrafo(), 2.3);
    auto chart         = unitMap.context().input();
    // Powers of two make the normalization round trip exact for these bounded inputs.
    const Point scales{2, 0.5, 4, 0.25};
    tps::VpsMap<double> polynomial(4, 4);
    for (unsigned i = 0; i < 4; ++i) {
        chart.coordinates[i].scale = scales[i];
        polynomial[i]              = Series(1, 4);
        polynomial[i][0]           = source.residual[i] / scales[i];
        for (unsigned j = 0; j < 4; ++j)
            polynomial[i][j + 1] = source.matrix[i][j] * scales[j] / scales[i];
    }
    const tps::Map<double> normalized(tps::MapContext(chart, chart, 1), polynomial);
    expectSpectrum(EV::analyze(normalized, settings), EV::analyze(source.matrix, settings));
    expectSpectrum(EV::analyze(unitMap, settings), EV::analyze(source.matrix, settings));
}

TEST_F(TpsClosedOrbitMapTest, SolverToTpsToSpectrumUsesNoAdditionalMapEvaluations) {
    const auto matrix = coupled(rotations(0.155, 0.31));
    const Point target{0.02, 0.003, -0.01, 0.004};
    unsigned calls    = 0;
    const auto solved = ClosedOrbitSolver::solve(
            [&](const Point& point) {
                ++calls;
                auto result = target;
                for (unsigned i = 0; i < 4; ++i)
                    for (unsigned j = 0; j < 4; ++j)
                        result[i] += matrix[i][j] * (point[j] - target[j]);
                return result;
            },
            Point{}, ClosedOrbitSolver::Settings{});
    ASSERT_EQ(solved.status, ClosedOrbitSolver::Status::Converged);
    const auto count = calls;
    EXPECT_EQ(count, solved.evaluations);
    const auto map    = opalx::maps::toTpsClosedOrbitMap(solved, CoordinateSystemTrafo(), 2.3);
    const auto result = EV::analyze(map, EV::Settings{});
    EXPECT_EQ(calls, count);
    expectSpectrum(result, EV::analyze(solved.matrix, EV::Settings{}));
    ASSERT_EQ(result.modes.size(), 2u);
    EXPECT_NEAR(*result.modes[0].tune, 0.155, 1e-10);
    EXPECT_NEAR(*result.modes[1].tune, 0.31, 1e-10);
}

TEST_F(TpsClosedOrbitMapTest, RejectsUnsupportedChartsAndNonfiniteCoefficients) {
    const auto map   = opalx::maps::toTpsClosedOrbitMap(orbit(), CoordinateSystemTrafo(), 2.3);
    const auto chart = map.context().input();
    for (unsigned change = 0; change < 7; ++change) {
        auto altered = chart;
        if (change == 0) altered.coordinates[1] = {"xprime", "1", 1};
        if (change == 1) {
            altered.convention     = tps::CoordinateConvention::canonical;
            altered.canonicalPairs = {{0, 1}, {2, 3}};
        }
        if (change == 2) {
            altered.referenceMomentum.reset();
            altered.momentumUnit.clear();
        }
        if (change == 3) altered.referenceTime = 0;
        if (change == 4) altered.coordinates[0].unit = "mm";
        if (change == 5) altered.referenceMomentum = 0.01;  // No real longitudinal momentum.
        if (change == 6) altered.frame.origin[0] += 0.1;
        const auto invalid = rechart(map, chart, altered);
        EXPECT_THROW((void)EV::analyze(invalid, EV::Settings{}), OpalException) << change;
        if (change < 6) {
            const auto sameChart = rechart(map, altered, altered);
            EXPECT_THROW((void)EV::analyze(sameChart, EV::Settings{}), OpalException) << change;
        }
    }
    const auto higherOrder = rechart(map, chart, chart, 2);
    EXPECT_THROW((void)EV::analyze(higherOrder, EV::Settings{}), OpalException);
    auto six = chart;
    six.coordinates.resize(6, {"extra", "1", 1});
    six.coordinates[5].name = "other";
    six.reference.resize(6);
    const auto sixDimensional = tps::identityMap<double>(six, 1);
    EXPECT_THROW((void)EV::analyze(sixDimensional, EV::Settings{}), OpalException);
    for (unsigned index : {0u, 1u}) {
        auto polynomial      = map.polynomial();
        polynomial[0][index] = std::numeric_limits<double>::quiet_NaN();
        const tps::Map<double> invalid(tps::MapContext(chart, chart, 1), polynomial);
        EXPECT_THROW((void)EV::analyze(invalid, EV::Settings{}), OpalException);
    }
}

TEST_F(TpsClosedOrbitMapTest, RejectsInvalidClosedOrbitAndFixedMomentumLaunch) {
    const auto valid = orbit();
    for (unsigned change = 0; change < 8; ++change) {
        auto invalid = valid;
        if (change == 0) invalid.status = ClosedOrbitSolver::Status::MaxIterations;
        if (change == 1) invalid.coordinates[0] = std::numeric_limits<double>::quiet_NaN();
        if (change == 2) invalid.residual[2] = std::numeric_limits<double>::infinity();
        if (change == 3) invalid.matrix[3][2] = std::numeric_limits<double>::quiet_NaN();
        if (change == 4) invalid.coordinates[1] = 2.3;
        if (change == 5) invalid.coordinates[3] = 3;
        if (change == 6) invalid.coordinates[1] = std::numeric_limits<double>::max();
        if (change == 7) invalid.matrix[0][1] = std::numeric_limits<double>::infinity();
        EXPECT_THROW(
                (void)opalx::maps::toTpsClosedOrbitMap(invalid, CoordinateSystemTrafo(), 2.3),
                OpalException)
                << change;
    }
    for (double momentum :
         {0.0, -1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::max()})
        EXPECT_THROW(
                (void)opalx::maps::toTpsClosedOrbitMap(valid, CoordinateSystemTrafo(), momentum),
                OpalException);
    const CoordinateSystemTrafo invalidSection(
            Vector_t<double, 3>(std::numeric_limits<double>::infinity(), 0, 0), Quaternion());
    EXPECT_THROW((void)opalx::maps::toTpsClosedOrbitMap(valid, invalidSection, 2.3), OpalException);
}
