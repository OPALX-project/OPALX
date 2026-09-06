#include <limits>
#include <sstream>
#include "Algorithms/ClosedOrbitSolver.h"
#include "Algorithms/LinearMapEigenAnalysis.h"
#include "Utilities/OpalException.h"
#include "gtest/gtest.h"

namespace {
    using EV = LinearMapEigenAnalysis;
    using M  = EV::Matrix;
    M rotations(double q1, double q2, double r1 = 1, double r2 = 1) {
        M a{};
        for (unsigned k = 0; k < 2; ++k) {
            const double q = k ? q2 : q1, r = k ? r2 : r1;
            const double c  = r * std::cos(2 * std::acos(-1.0) * q),
                         s  = r * std::sin(2 * std::acos(-1.0) * q);
            a[2 * k][2 * k] = a[2 * k + 1][2 * k + 1] = c;
            a[2 * k][2 * k + 1]                       = s;
            a[2 * k + 1][2 * k]                       = -s;
        }
        return a;
    }
    M coupled(const M& a) {
        const double c = std::cos(0.6), s = std::sin(0.6);
        const M q{{{c, 0, s, 0}, {0, c, 0, s}, {-s, 0, c, 0}, {0, -s, 0, c}}};
        M b{};
        for (unsigned i = 0; i < 4; ++i)
            for (unsigned j = 0; j < 4; ++j)
                for (unsigned k = 0; k < 4; ++k)
                    for (unsigned l = 0; l < 4; ++l)
                        b[i][j] += q[i][k] * a[k][l] * q[j][l];
        return b;
    }
}  // namespace

TEST(LinearMapEigenAnalysis, StableRotationsAndBranchAmbiguity) {
    const M a = rotations(0.155, 0.89), copy = a;
    const auto result = EV::analyze(a, EV::Settings{});
    EXPECT_EQ(a, copy);
    ASSERT_EQ(result.stability, EV::Stability::Stable);
    ASSERT_EQ(result.modes.size(), 2u);
    EXPECT_NEAR(*result.modes[0].tune, 0.11, 1e-13);
    EXPECT_NEAR(*result.modes[0].complementaryTune, 0.89, 1e-13);
    EXPECT_NEAR(*result.modes[1].tune, 0.155, 1e-13);
    EXPECT_LT(result.maximumResidual, 1e-13);
    for (auto lambda : result.eigenvalues)
        EXPECT_NEAR(std::abs(lambda), 1, 1e-13);
}

TEST(LinearMapEigenAnalysis, CoupledModesAndCoordinateScaling) {
    auto matrix     = coupled(rotations(0.155, 0.31));
    auto settings   = EV::Settings{};
    settings.scales = {100, 0.001, 10, 0.01};
    for (unsigned i = 0; i < 4; ++i)
        for (unsigned j = 0; j < 4; ++j)
            matrix[i][j] *= settings.scales[i] / settings.scales[j];
    const auto result = EV::analyze(matrix, settings);
    ASSERT_EQ(result.stability, EV::Stability::Stable);
    ASSERT_EQ(result.modes.size(), 2u);
    EXPECT_NEAR(*result.modes[0].tune, 0.155, 1e-13);
    EXPECT_NEAR(*result.modes[1].tune, 0.31, 1e-13);
    EXPECT_LT(result.maximumResidual, 1e-13);
}

TEST(LinearMapEigenAnalysis, RealAndComplexInstabilities) {
    auto a  = rotations(0.15, 0.3);
    a[0][0] = 2;
    a[1][1] = 0.5;
    a[0][1] = a[1][0] = 0;
    EXPECT_EQ(EV::analyze(a, EV::Settings{}).stability, EV::Stability::Unstable);
    const auto complex = EV::analyze(coupled(rotations(0.2, 0.2, 1.1, 1 / 1.1)), EV::Settings{});
    EXPECT_EQ(complex.stability, EV::Stability::Unstable);
    ASSERT_EQ(complex.modes.size(), 2u);
    for (const auto& mode : complex.modes)
        EXPECT_FALSE(mode.tune.has_value());
    EXPECT_EQ(
            EV::analyze(rotations(0.2, 0.3, 0.9, 0.8), EV::Settings{}).stability,
            EV::Stability::NonUnitCircle);
}

TEST(LinearMapEigenAnalysis, NearIntegerAndNeutralRealModes) {
    const auto near = EV::analyze(rotations(1e-8, 0.3), EV::Settings{});
    EXPECT_EQ(near.stability, EV::Stability::Marginal);
    EXPECT_TRUE(near.nearInteger);
    ASSERT_EQ(near.modes.size(), 2u);
    EXPECT_TRUE(near.modes[0].nearInteger);
    EXPECT_NEAR(*near.modes[0].tune, 1e-8, 1e-14);
    M identity{};
    for (unsigned i = 0; i < 4; ++i)
        identity[i][i] = 1;
    const auto unit = EV::analyze(identity, EV::Settings{});
    EXPECT_EQ(unit.stability, EV::Stability::Marginal);
    EXPECT_TRUE(unit.nearInteger);
    EXPECT_TRUE(unit.modes.empty());
    identity[0][0] = identity[1][1] = -1;
    EXPECT_EQ(EV::analyze(identity, EV::Settings{}).stability, EV::Stability::Marginal);
}

TEST(LinearMapEigenAnalysis, DefectiveComplexUnitCircleMapIsNotStable) {
    auto a  = rotations(0.2, 0.2);
    a[0][2] = a[1][3] = 1;
    const auto result = EV::analyze(a, EV::Settings{});
    EXPECT_EQ(result.stability, EV::Stability::Marginal);
    EXPECT_LT(result.basisReciprocalCondition, 1e-10);
}

TEST(LinearMapEigenAnalysis, SolverToEigenanalysisEndToEnd) {
    const auto matrix = coupled(rotations(0.155, 0.31));
    const OneTurnMap::Coordinates target{0.02, 0.003, -0.01, 0.004};
    const auto orbit = ClosedOrbitSolver::solve(
            [&](const auto& u) {
                auto v = target;
                for (unsigned i = 0; i < 4; ++i)
                    for (unsigned j = 0; j < 4; ++j)
                        v[i] += matrix[i][j] * (u[j] - target[j]);
                return v;
            },
            {0, 0, 0, 0}, ClosedOrbitSolver::Settings{});
    ASSERT_EQ(orbit.status, ClosedOrbitSolver::Status::Converged);
    const auto result = EV::analyze(orbit.matrix, EV::Settings{});
    ASSERT_EQ(result.stability, EV::Stability::Stable);
    ASSERT_EQ(result.modes.size(), 2u);
    EXPECT_NEAR(*result.modes[0].tune, 0.155, 1e-10);
    EXPECT_NEAR(*result.modes[1].tune, 0.31, 1e-10);
}

TEST(LinearMapEigenAnalysis, InvalidInputIsRejected) {
    auto a  = rotations(0.2, 0.3);
    a[0][0] = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(EV::analyze(a, EV::Settings{}), OpalException);
    auto settings      = EV::Settings{};
    settings.scales[0] = 0;
    EXPECT_THROW(EV::analyze(rotations(0.2, 0.3), settings), OpalException);
}

TEST(LinearMapEigenAnalysis, ReportExplainsBranchesAndResonance) {
    std::ostringstream report;
    EV::writeReport(report, EV::analyze(rotations(1e-8, 0.3), EV::Settings{}));
    EXPECT_NE(report.str().find("MARGINAL"), std::string::npos);
    EXPECT_NE(report.str().find("Near-integer warning: yes"), std::string::npos);
    EXPECT_NE(report.str().find("conjugate branch"), std::string::npos);
    EXPECT_NE(report.str().find("no x/y assignment"), std::string::npos);
}
