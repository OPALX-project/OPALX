#include "SpaceCharge/Analytic3D/FieldTensor.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <vector>

#include "Physics/Physics.h"
#include "Utilities/OpalException.h"

// Keep LAPACKE's unused complex declarations portable in CUDA/HIP builds too.
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

namespace opalx::spacecharge::analytic3d {
    namespace {
        constexpr double epsilon = std::numeric_limits<double>::epsilon();

        void require(bool condition, const char* message) {
            if (!condition) throw OpalException("analytic3d::evaluateRestFrameTensor", message);
        }

        struct PrincipalAxes {
            double scale;
            Vector variances;
            double rotation[9];  // Row-major orthogonal matrix, eigenvectors in columns.
        };

        PrincipalAxes principalAxes(const Matrix& covariance) {
            PrincipalAxes axes{};
            for (unsigned i = 0; i < 3; ++i) {
                require(std::isfinite(covariance[i][i]) && covariance[i][i] > 0,
                        "Covariance diagonal must be positive finite square metres.");
                axes.scale = std::max(axes.scale, covariance[i][i]);
            }
            for (unsigned i = 0; i < 3; ++i) {
                for (unsigned j = i; j < 3; ++j) {
                    const double a = covariance[i][j], b = covariance[j][i];
                    require(std::isfinite(a) && std::isfinite(b), "Nonfinite covariance entry.");
                    const double localScale =
                            std::sqrt(covariance[i][i]) * std::sqrt(covariance[j][j]);
                    require(std::abs(a - b) <= 64 * epsilon * localScale,
                            "Spatial covariance must be symmetric.");
                    const double value = 0.5 * (a / axes.scale) + 0.5 * (b / axes.scale);
                    require(std::isfinite(value), "Scaled covariance is not representable.");
                    axes.rotation[3 * i + j] = axes.rotation[3 * j + i] = value;
                }
            }
            require(LAPACKE_dsyev(
                            LAPACK_ROW_MAJOR, 'V', 'U', 3, axes.rotation, 3, axes.variances.data())
                            == 0,
                    "Covariance eigendecomposition failed.");
            for (double value : axes.variances) {
                require(std::isfinite(value) && value > 0,
                        "Spatial covariance must be numerically positive definite.");
            }
            return axes;
        }

        // s/scale = (t/(1-t))^2 maps [0,infinity) to [0,1]. These are the
        // depolarization-factor integrands, not the large dimensional field values.
        Vector integrand(double t, const Vector& variances) {
            const double remaining = 1 - t;
            Vector denominator{};
            double product = 1;
            for (unsigned i = 0; i < 3; ++i) {
                denominator[i] = variances[i] * remaining * remaining + t * t;
                product *= std::sqrt(variances[i] / denominator[i]);
            }
            Vector value{};
            for (unsigned i = 0; i < 3; ++i) {
                value[i] = t * remaining * remaining / denominator[i] * product;
                require(std::isfinite(value[i]), "Ellipsoid integrand is not representable.");
            }
            return value;
        }

        struct Interval {
            double lower, upper;
            Vector integral{}, error{};
        };

        // Embedded Gauss 7 / Kronrod 15 rule. QUADPACK-style rescaling avoids
        // treating a small raw embedded difference as an exact error bound.
        Interval integrate(double lower, double upper, const Vector& variances) {
            constexpr double nodes[8]   = {0.9914553711208126392, 0.9491079123427585245,
                                           0.8648644233597690728, 0.7415311855993944399,
                                           0.5860872354676911303, 0.4058451513773971669,
                                           0.2077849550078984676, 0.0};
            constexpr double kronrod[8] = {0.0229353220105292250, 0.0630920926299785533,
                                           0.1047900103222501838, 0.1406532597155259187,
                                           0.1690047266392679028, 0.1903505780647854099,
                                           0.2044329400752988924, 0.2094821410847278280};
            constexpr double gauss[4]   = {
                    0.1294849661688696933, 0.2797053914892766679, 0.3818300505051189449,
                    0.4179591836734693878};
            const double half = 0.5 * (upper - lower), midpoint = lower + half;
            const Vector center = integrand(midpoint, variances);
            Vector left[7], right[7];
            for (unsigned n = 0; n < 7; ++n) {
                left[n]  = integrand(midpoint - half * nodes[n], variances);
                right[n] = integrand(midpoint + half * nodes[n], variances);
            }
            Interval result{lower, upper};
            for (unsigned i = 0; i < 3; ++i) {
                double k = kronrod[7] * center[i], g = gauss[3] * center[i];
                for (unsigned n = 0; n < 7; ++n) {
                    const double pair = left[n][i] + right[n][i];
                    k += kronrod[n] * pair;
                    if (n % 2 == 1) g += gauss[n / 2] * pair;
                }
                double deviation = kronrod[7] * std::abs(center[i] - 0.5 * k);
                for (unsigned n = 0; n < 7; ++n) {
                    deviation +=
                            kronrod[n]
                            * (std::abs(left[n][i] - 0.5 * k) + std::abs(right[n][i] - 0.5 * k));
                }
                deviation *= half;
                const double rawError = half * std::abs(k - g);
                double error          = rawError;
                if (deviation > 0 && rawError > 0)
                    error = deviation * std::min(1.0, std::pow(200 * rawError / deviation, 1.5));
                result.integral[i] = half * k;
                // All integrands are positive, so the integral is also the L1 norm.
                result.error[i] = std::max(error, 50 * epsilon * result.integral[i]);
            }
            return result;
        }

        struct GeometryIntegral {
            Vector value{}, error{};
            std::size_t evaluations = 0, intervals = 0;
        };

        GeometryIntegral depolarization(
                const Vector& variances, const IntegrationSettings& settings) {
            std::vector<Interval> intervals;
            // Resolve every factor-of-two length scale before error-driven refinement;
            // an unsplit [0,1] rule can miss a very flat or narrow ellipsoid's endpoint layer.
            double lower = 0, length = std::sqrt(variances[0]);
            while (true) {
                const double upper = length < 1 ? length / (1 + length) : 0.5;
                require(intervals.size() < settings.maximumIntervals,
                        "Quadrature interval limit cannot resolve the ellipsoid aspect ratio.");
                intervals.push_back(integrate(lower, upper, variances));
                lower = upper;
                if (length >= 1) break;
                length = std::min(1.0, 2 * length);
            }
            require(intervals.size() < settings.maximumIntervals,
                    "Quadrature interval limit cannot resolve the ellipsoid tail.");
            intervals.push_back(integrate(0.5, 1, variances));
            GeometryIntegral result;
            result.evaluations = 15 * intervals.size();
            while (true) {
                result.value = {};
                result.error = {};
                for (const auto& interval : intervals) {
                    for (unsigned i = 0; i < 3; ++i) {
                        result.value[i] += interval.integral[i];
                        result.error[i] += interval.error[i];
                    }
                }
                bool converged = true;
                for (unsigned i = 0; i < 3; ++i) {
                    require(std::isfinite(result.value[i]) && result.value[i] > 0
                                    && std::isfinite(result.error[i]),
                            "Ellipsoid integral is not representable.");
                    converged &= result.error[i] <= settings.relativeTolerance * result.value[i];
                }
                if (converged) break;
                require(intervals.size() < settings.maximumIntervals,
                        "Ellipsoid quadrature failed to reach the requested relative tolerance.");
                std::size_t worst = 0;
                double worstError = -1;
                for (std::size_t n = 0; n < intervals.size(); ++n) {
                    for (unsigned i = 0; i < 3; ++i) {
                        const double relative = intervals[n].error[i] / result.value[i];
                        if (relative > worstError) {
                            worstError = relative;
                            worst      = n;
                        }
                    }
                }
                const auto interval   = intervals[worst];
                const double midpoint = interval.lower + 0.5 * (interval.upper - interval.lower);
                require(midpoint > interval.lower && midpoint < interval.upper,
                        "Ellipsoid quadrature reached floating-point resolution.");
                intervals[worst] = integrate(interval.lower, midpoint, variances);
                intervals.push_back(integrate(midpoint, interval.upper, variances));
                result.evaluations += 30;
            }
            result.intervals = intervals.size();
            return result;
        }
    }  // namespace

    TensorResult evaluateRestFrameTensor(
            Profile profile, const Matrix& covariance, double totalCharge,
            const IntegrationSettings& settings) {
        require(profile == Profile::UniformEllipsoid || profile == Profile::GaussianRms,
                "Unknown analytic 3D profile.");
        require(std::isfinite(settings.relativeTolerance) && settings.relativeTolerance >= 1e-13
                        && settings.relativeTolerance < 1e-2 && settings.maximumIntervals > 0,
                "Invalid ellipsoid integration tolerance or interval limit.");
        require(std::isfinite(totalCharge), "Total physical bunch charge must be finite coulombs.");
        const auto axes = principalAxes(covariance);
        TensorResult result;
        if (totalCharge == 0) return result;
        const auto geometry      = depolarization(axes.variances, settings);
        const double coefficient = profile == Profile::UniformEllipsoid
                                           ? 3 / (5 * std::sqrt(5.0))
                                           : 1 / (2 * std::sqrt(Physics::pi));
        const double coulomb     = 1 / (4 * Physics::pi * Physics::epsilon_0);
        // Logs avoid overflowing an intermediate product of three RMS lengths,
        // even when the final SI field coefficient is representable.
        double logScale = std::log(coulomb) + std::log(std::abs(totalCharge))
                          + std::log(coefficient) - 1.5 * std::log(axes.scale);
        for (double variance : axes.variances)
            logScale -= 0.5 * std::log(variance);
        Vector principal{}, error{};
        for (unsigned i = 0; i < 3; ++i) {
            principal[i] =
                    std::copysign(std::exp(logScale + std::log(geometry.value[i])), totalCharge);
            error[i] = std::exp(logScale + std::log(geometry.error[i]));
            require(std::isfinite(principal[i]) && principal[i] != 0 && std::isfinite(error[i]),
                    "Electric-field coefficient overflows or underflows double precision.");
        }
        for (unsigned i = 0; i < 3; ++i) {
            for (unsigned j = i; j < 3; ++j) {
                double value = 0, estimate = 0;
                for (unsigned n = 0; n < 3; ++n) {
                    const double product = axes.rotation[3 * i + n] * axes.rotation[3 * j + n];
                    value += product * principal[n];
                    estimate += std::abs(product) * error[n];
                }
                require(std::isfinite(value) && std::isfinite(estimate),
                        "Rotated electric-field tensor overflows double precision.");
                result.gradient[i][j] = result.gradient[j][i] = value;
                result.absoluteError[i][j] = result.absoluteError[j][i] = estimate;
            }
        }
        result.integrandEvaluations = geometry.evaluations;
        result.integrationIntervals = geometry.intervals;
        return result;
    }
}  // namespace opalx::spacecharge::analytic3d
