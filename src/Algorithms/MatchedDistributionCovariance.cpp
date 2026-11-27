#include "Algorithms/MatchedDistributionCovariance.h"
#include <cmath>
#include <complex>
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"
#include "tps/LinearSolve.hpp"
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
    void require(bool valid, const char* message) {
        if (!valid) throw OpalException("MOMENTS", message);
    }
    HostMatrix host(const MatchedDistributionCovariance::Matrix& stored) {
        HostMatrix result(6, 6);
        for (unsigned i = 0; i < 6; ++i)
            for (unsigned j = 0; j < 6; ++j)
                result(i, j) = stored[i][j];
        return result;
    }
}  // namespace
MatchedDistributionCovariance::Matrix MatchedDistributionCovariance::sectionProjection(
        const ClosedOrbitOptics& optics) {
    const auto& launch = optics.launch;
    const auto section = launch.section.getRotationMatrix();
    const auto frame   = launch.frame().getRotationMatrix();
    HostMatrix r(3, 3), f(3, 3);
    for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j) {
            r(i, j) = section(i, j);
            f(i, j) = frame(i, j);
        }
    const auto q = r * f.transpose();
    const auto u = launch.section.rotateTo(launch.momentum);
    require(std::isfinite(u[2]) && u[2] > 0,
            "MOMENTS reference must cross the section in the positive direction.");
    require(std::isfinite(launch.massEV) && launch.massEV > 0 && std::isfinite(launch.chargeE)
                    && launch.chargeE != 0,
            "Invalid MOMENTS reference species.");
    const double kappa = launch.chargeE * Physics::c / launch.massEV;
    const double bx = optics.transverseMagneticField[0], by = optics.transverseMagneticField[1],
                 bz = optics.normalMagneticField;
    const double gx = kappa * (u[1] * bz - u[2] * by) / u[2];
    const double gy = kappa * (u[2] * bx - u[0] * bz) / u[2];
    Matrix g{};
    for (unsigned i = 0; i < 3; ++i) {
        g[0][2 * i]     = q(0, i) - u[0] / u[2] * q(2, i);
        g[1][2 * i]     = -gx * q(2, i);
        g[1][2 * i + 1] = q(0, i);
        g[2][2 * i]     = q(1, i) - u[1] / u[2] * q(2, i);
        g[3][2 * i]     = -gy * q(2, i);
        g[3][2 * i + 1] = q(1, i);
    }
    g[4][4] = g[5][5] = 1;
    for (const auto& row : g)
        for (double value : row)
            require(std::isfinite(value), "Nonfinite MOMENTS section conversion.");
    return g;
}
MatchedDistributionCovariance::Matrix MatchedDistributionCovariance::compute(
        const ClosedOrbitOptics& optics, const OneTurnMap::Matrix& transverse, double sigmaZ,
        double sigmaPZ) {
    require(std::isfinite(sigmaZ) && sigmaZ > 0 && std::isfinite(sigmaPZ) && sigmaPZ > 0,
            "MOMENTS requires positive finite SIGMAZ [m] and SIGMAPZ [p/(mc)].");
    const auto g = host(sectionProjection(optics));
    HostMatrix model(6, 6, 0.0);
    for (unsigned i = 0; i < 4; ++i)
        for (unsigned j = 0; j < 4; ++j) {
            require(std::isfinite(transverse[i][j]) && transverse[i][j] == transverse[j][i],
                    "Invalid transverse MOMENTS covariance.");
            model(i, j) = transverse[i][j];
        }
    model(4, 4)        = sigmaZ * sigmaZ;
    model(5, 5)        = sigmaPZ * sigmaPZ;
    const auto inverse = tps::LUFactorization<double>(g).inverse();
    auto covariance    = inverse * model * inverse.transpose();
    covariance         = (covariance + covariance.transpose()) * 0.5;
    Matrix result{};
    std::array<double, 36> correlation{};
    for (unsigned i = 0; i < 6; ++i)
        require(std::isfinite(covariance(i, i)) && covariance(i, i) > 0,
                "Nonpositive MOMENTS covariance diagonal.");
    for (unsigned i = 0; i < 6; ++i)
        for (unsigned j = 0; j < 6; ++j) {
            require(std::isfinite(covariance(i, j)), "Nonfinite MOMENTS covariance.");
            result[i][j] = covariance(i, j);
            correlation[6 * i + j] =
                    covariance(i, j) / std::sqrt(covariance(i, i)) / std::sqrt(covariance(j, j));
        }
    require(LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'L', 6, correlation.data(), 6) == 0,
            "MOMENTS covariance is not positive definite.");
    return result;
}
