// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#include "Algorithms/ClosedOrbitSolver.h"
#include <algorithm>
#include <cmath>
#include "Utilities/OpalException.h"

namespace {
    using V = ClosedOrbitSolver::Coordinates;
    using M = ClosedOrbitSolver::Matrix;
    void require(bool valid, const char* message) {
        if (!valid) throw OpalException("ClosedOrbitSolver", message);
    }
    bool finite(const V& v) {
        return std::all_of(v.begin(), v.end(), [](double x) {
            return std::isfinite(x);
        });
    }
    double norm(const V& v) {
        double n = 0;
        for (double x : v)
            n = std::hypot(n, x);
        return n;
    }
    // Four-column Householder QR with exact recomputation of trailing column
    // norms (no downdate cancellation). Apply reflectors to RHS, back substitute,
    // and undo column permutation. Rank rejection is relative to the first pivot.
    bool correction(M a, V b, double cutoff, V& x, double& ratio) {
        std::array<unsigned, 4> permutation{0, 1, 2, 3};
        double largest = 0, smallest = 0;
        for (unsigned k = 0; k < 4; ++k) {
            unsigned pivot = k;
            double best    = -1;
            for (unsigned j = k; j < 4; ++j) {
                double n = 0;
                for (unsigned i = k; i < 4; ++i)
                    n = std::hypot(n, a[i][j]);
                if (n > best) {
                    best  = n;
                    pivot = j;
                }
            }
            if (k == 0) largest = best;
            smallest = k == 0 ? best : std::min(smallest, best);
            ratio    = largest > 0 ? smallest / largest : 0;
            if (!std::isfinite(best) || !(best > 0) || ratio <= cutoff) return false;
            for (unsigned i = 0; i < 4; ++i)
                std::swap(a[i][k], a[i][pivot]);
            std::swap(permutation[k], permutation[pivot]);
            // Normalize before forming v, avoiding squared matrix magnitudes.
            V v{};
            for (unsigned i = k; i < 4; ++i)
                v[i] = a[i][k] / best;
            v[k] += std::copysign(1.0, v[k]);
            const double vn = norm(v);
            for (double& value : v)
                value /= vn;
            for (unsigned j = k; j < 4; ++j) {
                double dot = 0;
                for (unsigned i = k; i < 4; ++i)
                    dot += v[i] * a[i][j];
                for (unsigned i = k; i < 4; ++i)
                    a[i][j] -= 2 * v[i] * dot;
            }
            double dot = 0;
            for (unsigned i = k; i < 4; ++i)
                dot += v[i] * b[i];
            for (unsigned i = k; i < 4; ++i)
                b[i] -= 2 * v[i] * dot;
        }
        V z{};
        for (int i = 3; i >= 0; --i) {
            double rhs = b[i];
            for (unsigned j = i + 1; j < 4; ++j)
                rhs -= a[i][j] * z[j];
            z[i] = rhs / a[i][i];
        }
        for (unsigned i = 0; i < 4; ++i)
            x[permutation[i]] = z[i];
        return finite(x);
    }
}  // namespace

ClosedOrbitSolver::Result ClosedOrbitSolver::solve(
        const OneTurnMap& map, const Coordinates& initial, const Settings& settings) {
    return solve(
            Map([&](const Coordinates& u) {
                return map(u).coordinates;
            }),
            initial, settings);
}

ClosedOrbitSolver::Result ClosedOrbitSolver::solve(
        const Map& map, const Coordinates& initial, const Settings& s) {
    require(bool(map) && finite(initial), "Missing map or nonfinite initial coordinates.");
    for (unsigned i = 0; i < 4; ++i)
        require(std::isfinite(s.scales[i]) && s.scales[i] > 0
                        && std::isfinite(2 * s.finiteDifferenceSteps[i])
                        && s.finiteDifferenceSteps[i] > 0,
                "Scales and finite-difference steps must be finite and positive.");
    require(std::isfinite(s.positionTolerance) && s.positionTolerance > 0
                    && std::isfinite(s.momentumTolerance) && s.momentumTolerance > 0
                    && std::isfinite(s.rankTolerance) && s.rankTolerance > 0 && s.rankTolerance < 1
                    && s.maxIterations > 0 && s.maxBacktracks <= 1000,
            "Invalid solver controls.");
    Result result;
    result.coordinates  = initial;
    const auto evaluate = [&](const V& u) {
        require(finite(u), "Nonfinite trial coordinates.");
        ++result.evaluations;
        const V value = map(u);
        require(finite(value), "Map returned nonfinite coordinates.");
        return value;
    };
    const auto residual = [&](const V& u) {
        V f = evaluate(u);
        for (unsigned i = 0; i < 4; ++i)
            f[i] -= u[i];
        require(finite(f), "Nonfinite closure residual.");
        return f;
    };
    const auto small = [&](const V& v) {
        for (unsigned i = 0; i < 4; ++i)
            if (std::abs(v[i]) > (i % 2 == 0 ? s.positionTolerance : s.momentumTolerance))
                return false;
        return true;
    };
    const auto scaledNorm = [&](V v) {
        for (unsigned i = 0; i < 4; ++i)
            v[i] /= s.scales[i];
        require(finite(v) && std::isfinite(norm(v)), "Nonfinite scaled residual.");
        return norm(v);
    };
    try {
        result.residual = residual(initial);
        for (unsigned iteration = 0;; ++iteration) {
            for (unsigned j = 0; j < 4; ++j) {
                V plus = result.coordinates, minus = plus;
                const double h = s.finiteDifferenceSteps[j];
                plus[j] += h;
                minus[j] -= h;
                require(plus[j] != result.coordinates[j] && minus[j] != result.coordinates[j],
                        "Finite-difference step below coordinate resolution.");
                const V a = evaluate(plus), b = evaluate(minus);
                for (unsigned i = 0; i < 4; ++i)
                    result.matrix[i][j] = (a[i] - b[i]) / (2 * h);
            }
            M a = result.matrix;
            V rhs{}, delta{};
            for (unsigned i = 0; i < 4; ++i) {
                rhs[i] = -result.residual[i] / s.scales[i];
                for (unsigned j = 0; j < 4; ++j)
                    a[i][j] = (a[i][j] - (i == j ? 1 : 0)) * s.scales[j] / s.scales[i];
                require(finite(a[i]) && std::isfinite(rhs[i]), "Nonfinite scaled Newton system.");
            }
            if (!correction(a, rhs, s.rankTolerance, delta, result.pivotRatio)) {
                result.status  = Status::SingularJacobian;
                result.message = "Rank-deficient or numerically unsolvable scaled (M-I).";
                return result;
            }
            for (unsigned i = 0; i < 4; ++i)
                delta[i] *= s.scales[i];
            require(finite(delta), "Nonfinite Newton correction.");
            if (small(result.residual) && small(delta)) {
                result.status   = Status::VerificationFailed;
                result.residual = residual(result.coordinates);  // Fresh, not cached.
                if (small(result.residual)) {
                    result.status  = Status::Converged;
                    result.message = "Closed orbit verified by a fresh map evaluation.";
                } else
                    result.message = "Fresh traversal did not verify closure.";
                return result;
            }
            if (iteration == s.maxIterations) {
                result.status  = Status::MaxIterations;
                result.message = "Maximum accepted Newton iterations reached.";
                return result;
            }
            const double merit = scaledNorm(result.residual);
            bool accepted      = false;
            std::string lastLoss;
            const unsigned backtracks = s.damping ? s.maxBacktracks : 0;
            for (unsigned backtrack = 0; backtrack <= backtracks; ++backtrack) {
                const double lambda = std::ldexp(1.0, -static_cast<int>(backtrack));
                V trial             = result.coordinates, step{};
                for (unsigned i = 0; i < 4; ++i) {
                    step[i] = lambda * delta[i];
                    trial[i] += step[i];
                }
                try {
                    const V f = residual(trial);
                    if (scaledNorm(f) <= (1 - 1e-4 * lambda) * merit) {
                        result.coordinates = trial;
                        result.residual    = f;
                        result.iterations.push_back({trial, f, step, lambda, result.pivotRatio});
                        accepted = true;
                        break;
                    }
                } catch (const OpalException& error) {
                    lastLoss = error.what();
                }
            }
            if (!accepted) {
                result.status  = Status::LineSearchFailed;
                result.message = "No admissible residual-reducing Newton trial.";
                if (!lastLoss.empty()) result.message += " Last trial failure: " + lastLoss;
                return result;
            }
        }
    } catch (const OpalException& error) {
        if (result.status != Status::VerificationFailed) result.status = Status::EvaluationFailed;
        result.message = error.what();
        return result;
    }
}
