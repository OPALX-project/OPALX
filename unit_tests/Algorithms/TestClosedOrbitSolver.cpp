#include <limits>
#include <random>
#include "Algorithms/ClosedOrbitSolver.h"
#include "Utilities/OpalException.h"
#include "gtest/gtest.h"

namespace {
    using Solver = ClosedOrbitSolver;
    using V      = Solver::Coordinates;
    using M      = Solver::Matrix;
    Solver::Map affine(const M& matrix, const V& fixed) {
        return [=](const V& u) {
            V v = fixed;
            for (unsigned i = 0; i < 4; ++i)
                for (unsigned j = 0; j < 4; ++j)
                    v[i] += matrix[i][j] * (u[j] - fixed[j]);
            return v;
        };
    }
}  // namespace

TEST(ClosedOrbitSolver, CoupledAffineMapAndFinalJacobian) {
    // Two stable rotations conjugated by an orthogonal cross-plane rotation.
    const double c = std::cos(0.4), s = std::sin(0.4), d = std::cos(1.2), t = std::sin(1.2);
    M matrix{
            {{(c + d) / 2, (s + t) / 2, (c - d) / 2, (s - t) / 2},
             {-(s + t) / 2, (c + d) / 2, (-s + t) / 2, (c - d) / 2},
             {(c - d) / 2, (s - t) / 2, (c + d) / 2, (s + t) / 2},
             {(-s + t) / 2, (c - d) / 2, -(s + t) / 2, (c + d) / 2}}};
    const V target{0.02, 0.003, -0.01, 0.004};
    const auto result = Solver::solve(affine(matrix, target), {0, 0, 0, 0}, Solver::Settings{});
    ASSERT_EQ(result.status, Solver::Status::Converged) << result.message;
    EXPECT_GT(result.evaluations, 9u);
    EXPECT_LE(result.iterations.size(), 2u);
    EXPECT_GT(result.pivotRatio, 0);
    for (unsigned i = 0; i < 4; ++i) {
        EXPECT_NEAR(result.coordinates[i], target[i], 1e-10);
        EXPECT_NEAR(result.residual[i], 0, 1e-10);
        for (unsigned j = 0; j < 4; ++j)
            EXPECT_NEAR(result.matrix[i][j], matrix[i][j], 1e-9);
    }
}

TEST(ClosedOrbitSolver, CoordinateScalingAndColumnPermutation) {
    Solver::Settings settings;
    settings.scales = {1000, 0.001, 10, 0.01};
    V target{2, 0.0002, -0.1, 0.003};
    // Unequal columns force QR pivoting. Similarity scaling changes units only.
    M matrix{{{0.2, 0.8, 0, 0}, {-0.3, 0.1, 0, 0}, {0, 0, 0.7, 0.2}, {0, 0, -0.1, 0.4}}};
    for (unsigned i = 0; i < 4; ++i) {
        settings.finiteDifferenceSteps[i] = 1e-5 * settings.scales[i];
        for (unsigned j = 0; j < 4; ++j)
            matrix[i][j] *= settings.scales[i] / settings.scales[j];
    }
    const auto result = Solver::solve(affine(matrix, target), {0, 0, 0, 0}, settings);
    ASSERT_EQ(result.status, Solver::Status::Converged) << result.message;
    for (unsigned i = 0; i < 4; ++i)
        EXPECT_NEAR(result.coordinates[i], target[i], 1e-10);
}

TEST(ClosedOrbitSolver, DampsNonlinearAndLostTrials) {
    unsigned losses = 0;
    Solver::Map map = [&](const V& u) {
        if (u[0] > 2) {
            ++losses;
            throw OpalException("test", "outside aperture");
        }
        V v{};
        v[0] = u[0] + u[0] * u[0] - 1;
        return v;
    };
    const auto result = Solver::solve(map, {0.1, 0, 0, 0}, Solver::Settings{});
    ASSERT_EQ(result.status, Solver::Status::Converged) << result.message;
    EXPECT_NEAR(result.coordinates[0], 1, 1e-10);
    EXPECT_GT(losses, 0u);
    ASSERT_FALSE(result.iterations.empty());
    EXPECT_LT(result.iterations.front().damping, 1);
    auto settings    = Solver::Settings{};
    settings.damping = false;
    EXPECT_EQ(
            Solver::solve(map, {0.1, 0, 0, 0}, settings).status, Solver::Status::LineSearchFailed);
}

TEST(ClosedOrbitSolver, RejectsSingularAndNearSingularJacobian) {
    auto identity = [](const V& u) {
        return u;
    };
    EXPECT_EQ(
            Solver::solve(identity, {0, 0, 0, 0}, Solver::Settings{}).status,
            Solver::Status::SingularJacobian);
    auto settings          = Solver::Settings{};
    settings.rankTolerance = 1e-6;
    M matrix{};
    matrix[0][0]      = 1 - 1e-8;
    const auto result = Solver::solve(affine(matrix, {0, 0, 0, 0}), {0, 0, 0, 0}, settings);
    EXPECT_EQ(result.status, Solver::Status::SingularJacobian);
    EXPECT_LT(result.pivotRatio, settings.rankTolerance);
}

TEST(ClosedOrbitSolver, IterationLimitAndEvaluationFailures) {
    auto settings          = Solver::Settings{};
    settings.maxIterations = 1;
    auto nonlinear         = [](const V& u) {
        V v{};
        v[0] = u[0] + u[0] * u[0] - 1;
        return v;
    };
    EXPECT_EQ(
            Solver::solve(nonlinear, {0.5, 0, 0, 0}, settings).status,
            Solver::Status::MaxIterations);
    auto loss = [](const V&) -> V {
        throw OpalException("test", "lost launch");
    };
    EXPECT_EQ(Solver::solve(loss, {0, 0, 0, 0}, settings).status, Solver::Status::EvaluationFailed);
    auto probeLoss = [](const V& u) -> V {
        if (u[1] != 0) throw OpalException("test", "lost derivative probe");
        return {};
    };
    EXPECT_EQ(
            Solver::solve(probeLoss, {0, 0, 0, 0}, settings).status,
            Solver::Status::EvaluationFailed);
    auto nan = [](const V&) {
        return V{std::numeric_limits<double>::quiet_NaN(), 0, 0, 0};
    };
    EXPECT_EQ(Solver::solve(nan, {0, 0, 0, 0}, settings).status, Solver::Status::EvaluationFailed);
}

TEST(ClosedOrbitSolver, FreshVerificationAndCorrectionCriterion) {
    unsigned calls = 0;
    auto changing  = [&](const V&) {
        ++calls;
        return calls == 10 ? V{0.01, 0, 0, 0} : V{};
    };
    const auto result = Solver::solve(changing, {0, 0, 0, 0}, Solver::Settings{});
    EXPECT_EQ(result.status, Solver::Status::VerificationFailed);
    EXPECT_EQ(result.evaluations, 10u);  // Initial + 8 probes + fresh turn.
    // Small residual alone is insufficient when M-I is weak in one direction.
    M matrix{};
    matrix[0][0]           = 1 - 1e-6;
    auto settings          = Solver::Settings{};
    settings.rankTolerance = 1e-9;
    const auto corrected   = Solver::solve(affine(matrix, {1e-5, 0, 0, 0}), {0, 0, 0, 0}, settings);
    ASSERT_EQ(corrected.status, Solver::Status::Converged) << corrected.message;
    EXPECT_FALSE(corrected.iterations.empty());
    EXPECT_NEAR(corrected.coordinates[0], 1e-5, 1e-10);
}

TEST(ClosedOrbitSolver, InvalidControlsThrow) {
    auto map = [](const V&) {
        return V{};
    };
    auto settings      = Solver::Settings{};
    settings.scales[2] = 0;
    EXPECT_THROW(Solver::solve(map, {0, 0, 0, 0}, settings), OpalException);
    settings                          = Solver::Settings{};
    settings.finiteDifferenceSteps[1] = -1;
    EXPECT_THROW(Solver::solve(map, {0, 0, 0, 0}, settings), OpalException);
}

TEST(ClosedOrbitSolver, DenseKnownSystemsExercisePivotedSolve) {
    std::mt19937 generator(544);
    std::uniform_real_distribution<double> random(-0.5, 0.5);
    for (unsigned sample = 0; sample < 40; ++sample) {
        M matrix{};
        V fixed{}, initial{};
        for (unsigned i = 0; i < 4; ++i) {
            fixed[i]   = random(generator);
            initial[i] = random(generator);
            for (unsigned j = 0; j < 4; ++j)
                matrix[i][j] = random(generator);
            matrix[i][i] -= 2;  // M-I strictly diagonally dominant: nonsingular oracle.
        }
        const auto result = Solver::solve(affine(matrix, fixed), initial, Solver::Settings{});
        ASSERT_EQ(result.status, Solver::Status::Converged) << sample << ' ' << result.message;
        for (unsigned i = 0; i < 4; ++i)
            EXPECT_NEAR(result.coordinates[i], fixed[i], 1e-10);
    }
}

TEST(ClosedOrbitSolver, FailedLineSearchPreservesLastAcceptedState) {
    auto settings          = Solver::Settings{};
    settings.maxBacktracks = 2;
    auto map               = [](const V& u) -> V {
        if (std::abs(u[0]) > 1e-4) throw OpalException("test", "outside aperture");
        return {1, 0, 0, 0};
    };
    const V initial{};
    const auto result = Solver::solve(map, initial, settings);
    EXPECT_EQ(result.status, Solver::Status::LineSearchFailed);
    EXPECT_EQ(result.coordinates, initial);
    EXPECT_TRUE(result.iterations.empty());
    EXPECT_NE(result.message.find("outside aperture"), std::string::npos);
}

TEST(ClosedOrbitSolver, CentralDifferenceJacobianHasSecondOrderError) {
    auto map = [](const V& u) {
        V v{};
        for (unsigned i = 0; i < 4; ++i)
            v[i] = 0.2 * u[i] + 0.1 * u[i] * u[i] * u[i];
        return v;
    };
    auto settings = Solver::Settings{};
    settings.finiteDifferenceSteps.fill(0.01);
    const auto coarse = Solver::solve(map, {0, 0, 0, 0}, settings);
    settings.finiteDifferenceSteps.fill(0.005);
    const auto fine = Solver::solve(map, {0, 0, 0, 0}, settings);
    ASSERT_EQ(coarse.status, Solver::Status::Converged);
    ASSERT_EQ(fine.status, Solver::Status::Converged);
    EXPECT_NEAR((coarse.matrix[0][0] - 0.2) / (fine.matrix[0][0] - 0.2), 4, 1e-8);
}

TEST(ClosedOrbitSolver, BacktrackingRejectsNonImprovingFiniteTrials) {
    auto map = [](const V& u) {
        V v{};
        v[0] = u[0] + u[0] * u[0] - 1;
        return v;
    };
    const auto result = Solver::solve(map, {0.1, 0, 0, 0}, Solver::Settings{});
    ASSERT_EQ(result.status, Solver::Status::Converged) << result.message;
    ASSERT_FALSE(result.iterations.empty());
    EXPECT_DOUBLE_EQ(result.iterations.front().damping, 0.25);
    EXPECT_NEAR(result.coordinates[0], 1, 1e-10);
}
