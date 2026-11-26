#include "Algorithms/MatchedCovariance.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include "Algorithms/LinearMapEigenAnalysis.h"
#include "Algorithms/TpsClosedOrbitMap.h"
#include "Physics/Physics.h"
#include "tps/Algebra/Matrix.h"

// LAPACKE's unused complex declarations must also be portable in device builds.
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
#endif
#include <lapacke.h>
#if defined(__clang__)
#pragma clang diagnostic pop
#endif
#undef lapack_complex_float
#undef lapack_complex_double

namespace {
    using HostMatrix = tps::Matrix<double>;
    using Matrix     = MatchedCovariance::Matrix;
    void require(bool valid, const char* message) {
        if (!valid) throw OpalException("MATCH", message);
    }
    HostMatrix host(const Matrix& a) {
        HostMatrix result(4, 4);
        for (unsigned i = 0; i < 4; ++i)
            for (unsigned j = 0; j < 4; ++j)
                result(i, j) = a[i][j];
        return result;
    }
    Matrix stored(const HostMatrix& a) {
        Matrix result{};
        for (unsigned i = 0; i < 4; ++i)
            for (unsigned j = 0; j < 4; ++j) {
                require(std::isfinite(a(i, j)), "Nonfinite matching matrix.");
                result[i][j] = a(i, j);
            }
        return result;
    }
    double norm(const HostMatrix& a) {
        double result = 0;
        for (unsigned i = 0; i < 4; ++i)
            for (unsigned j = 0; j < 4; ++j)
                result = std::hypot(result, a(i, j));
        require(std::isfinite(result), "Matching matrix norm overflow.");
        return result;
    }
    double residual(const HostMatrix& a, const HostMatrix& reference) {
        const double denominator = norm(reference);
        require(denominator > 0, "Vanishing matching residual scale.");
        return norm(a - reference) / denominator;
    }
}  // namespace

MatchedCovariance::Matrix MatchedCovariance::canonicalTransform(
        double field, double charge, double mass) {
    require(std::isfinite(field) && std::isfinite(charge) && charge != 0 && std::isfinite(mass)
                    && mass > 0,
            "Invalid reference field, charge or rest energy.");
    const double b = charge * Physics::c * field / mass;
    require(std::isfinite(b), "Canonical conversion overflow.");
    Matrix result{};
    for (unsigned i = 0; i < 4; ++i)
        result[i][i] = 1;
    result[1][2] = -0.5 * b;
    result[3][0] = 0.5 * b;
    return result;
}

MatchedCovariance::Result MatchedCovariance::compute(
        const ClosedOrbitOptics& optics, const std::array<double, 2>& eps,
        const Settings& settings) {
    for (double e : eps)
        require(std::isfinite(e) && e > 0, "Mode emittances must be positive finite metres.");
    for (double tolerance : {settings.symplecticTolerance, settings.periodicTolerance})
        require(std::isfinite(tolerance) && tolerance > 0 && tolerance < 1,
                "Invalid matching residual tolerance.");
    require(std::isfinite(settings.minimumModeSeparation) && settings.minimumModeSeparation > 0
                    && settings.minimumModeSeparation < 0.25,
            "Invalid mode separation.");
    require(std::isfinite(optics.relativeEnergyDrift)
                    && std::abs(optics.relativeEnergyDrift) <= 1e-8,
            "COF energy drift exceeds 1e-8; refine the static-magnetic return map.");
    const auto& launch    = optics.launch;
    const double momentum = std::hypot(launch.momentum[0], launch.momentum[1], launch.momentum[2]);
    ClosedOrbitSolver::Result linearization;
    linearization.status      = ClosedOrbitSolver::Status::Converged;
    linearization.coordinates = optics.coordinates;
    linearization.residual    = optics.residual;
    linearization.matrix      = optics.matrix;
    const auto map = opalx::maps::toTpsClosedOrbitMap(linearization, launch.section, momentum);
    const auto m   = map.jacobian(OneTurnMap::Coordinates{});
    const auto c =
            host(canonicalTransform(optics.normalMagneticField, launch.chargeE, launch.massEV));
    auto inverse         = HostMatrix::identity(4);
    inverse(1, 2)        = -c(1, 2);
    inverse(3, 0)        = -c(3, 0);
    const auto canonical = c * m * inverse;
    LinearMapEigenAnalysis::Settings ev;
    ev.scales           = optics.scales;
    const auto spectrum = LinearMapEigenAnalysis::analyze(stored(canonical), ev);
    require(spectrum.stability == LinearMapEigenAnalysis::Stability::Stable
                    && spectrum.modes.size() == 2,
            "MATCH requires two stable, well-conditioned oscillatory modes.");
    Result result;
    result.emittances = eps;
    for (unsigned k = 0; k < 2; ++k) {
        require(spectrum.modes[k].tune.has_value(), "Missing oscillatory mode phase.");
        result.tunes[k] = *spectrum.modes[k].tune;
        require(result.tunes[k] > settings.minimumModeSeparation
                        && 0.5 - result.tunes[k] > settings.minimumModeSeparation,
                "Marginal modes near integer or half-integer phase are unsupported.");
    }
    require(result.tunes[1] - result.tunes[0] > settings.minimumModeSeparation,
            "Degenerate or insufficiently separated mode phases are unsupported.");
    auto j  = HostMatrix(4, 4, 0.0);
    j(0, 1) = j(2, 3) = 1;
    j(1, 0) = j(3, 2) = -1;
    auto scaled       = canonical;
    auto form         = j;
    for (unsigned i = 0; i < 4; ++i)
        for (unsigned k = 0; k < 4; ++k) {
            scaled(i, k) *= ev.scales[k] / ev.scales[i];
            form(i, k) *= ev.scales[i] * ev.scales[k];
        }
    result.symplecticResidual = residual(scaled.transpose() * form * scaled, form);
    require(result.symplecticResidual <= settings.symplecticTolerance,
            "One-turn symplectic residual exceeds tolerance; refine COF DT/FDSTEP.");

    HostMatrix basis(4, 4, 0.0), target(4, 4, 0.0);
    for (unsigned k = 0; k < 2; ++k) {
        unsigned column = 4;
        for (unsigned n = 0; n < 4; ++n)
            if (spectrum.eigenvalues[n] == spectrum.modes[k].eigenvalue) column = n;
        require(column < 4, "Mode eigenvector is missing.");
        for (unsigned i = 0; i < 4; ++i) {
            const auto value    = spectrum.eigenvectors[column][i] * ev.scales[i];
            basis(i, 2 * k)     = value.real();
            basis(i, 2 * k + 1) = value.imag();
        }
        double area = 0;
        for (unsigned i : {0u, 2u})
            area += basis(i, 2 * k) * basis(i + 1, 2 * k + 1)
                    - basis(i + 1, 2 * k) * basis(i, 2 * k + 1);
        require(std::isfinite(area) && std::abs(area) > std::numeric_limits<double>::min(),
                "Vanishing symplectic eigenvector norm.");
        const double normalization = std::sqrt(std::abs(area));
        for (unsigned i = 0; i < 4; ++i) {
            basis(i, 2 * k) /= normalization;
            basis(i, 2 * k + 1) /= normalization;
        }
        target(2 * k, 2 * k + 1) = std::copysign(1.0, area);
        target(2 * k + 1, 2 * k) = -target(2 * k, 2 * k + 1);
    }
    result.normalizationResidual = residual(basis.transpose() * j * basis, target);
    require(result.normalizationResidual <= settings.symplecticTolerance,
            "Mode basis is not symplectically orthogonal within tolerance.");
    auto weighted = basis;
    for (unsigned i = 0; i < 4; ++i)
        for (unsigned k = 0; k < 4; ++k)
            weighted(i, k) *= eps[k / 2];
    auto sigma = inverse * (weighted * basis.transpose()) * inverse.transpose();
    // Remove only multiplication roundoff asymmetry, before checking periodicity.
    sigma                 = (sigma + sigma.transpose()) * 0.5;
    auto mechanicalScaled = m;
    auto sigmaScaled      = sigma;
    for (unsigned i = 0; i < 4; ++i)
        for (unsigned k = 0; k < 4; ++k) {
            mechanicalScaled(i, k) *= ev.scales[k] / ev.scales[i];
            sigmaScaled(i, k) /= ev.scales[i] * ev.scales[k];
        }
    result.periodicResidual =
            residual(mechanicalScaled * sigmaScaled * mechanicalScaled.transpose(), sigmaScaled);
    require(result.periodicResidual <= settings.periodicTolerance,
            "Covariance periodicity residual exceeds tolerance; refine COF DT/FDSTEP.");
    // Test positive definiteness after diagonal equilibration, avoiding mixed units.
    std::array<double, 16> correlation{};
    for (unsigned i = 0; i < 4; ++i) {
        require(std::isfinite(sigma(i, i)) && sigma(i, i) > 0, "Nonpositive covariance diagonal.");
        for (unsigned k = 0; k < 4; ++k)
            correlation[4 * i + k] = sigma(i, k) / std::sqrt(sigma(i, i)) / std::sqrt(sigma(k, k));
    }
    require(LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'L', 4, correlation.data(), 4) == 0,
            "Matched covariance is not positive definite.");
    result.covariance = stored(sigma);
    return result;
}
