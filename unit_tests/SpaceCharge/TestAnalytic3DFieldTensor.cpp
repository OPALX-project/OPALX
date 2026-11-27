/**
 * @file TestAnalytic3DFieldTensor.cpp
 * @brief Host analytic 3D tensors: physical normalization, geometry and failure diagnostics.
 *
 * No OPAL input, particles, mesh or random sampling. Analytic sphere/spheroid limits and an
 * independent radial Gaussian force average test the physical model, rather than just comparing
 * two instances of its quadrature. Default-tolerance comparisons allow 2e-10 relative error;
 * cancellation-sensitive off-diagonals use a tensor-scale absolute bound. No Kokkos/MPI runtime
 * is needed by the host evaluator or the value-only LorentzBoost crosscheck.
 */
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <limits>

#include "Physics/Physics.h"
#include "SpaceCharge/Analytic3D/FieldTensor.h"
#include "SpaceCharge/CartesianPIC3D/RelativisticFieldComposer.h"
#include "Utilities/OpalException.h"

namespace opalx::spacecharge::analytic3d {
    namespace {
        constexpr double coulomb     = 1 / (4 * Physics::pi * Physics::epsilon_0);
        constexpr Profile profiles[] = {Profile::UniformEllipsoid, Profile::GaussianRms};

        Matrix diagonal(double x, double y, double z) {
            Matrix result{};
            result[0][0] = x;
            result[1][1] = y;
            result[2][2] = z;
            return result;
        }

        double normalization(Profile profile, double charge, double rmsVolume) {
            return coulomb * charge / rmsVolume
                   * (profile == Profile::UniformEllipsoid ? 3 / (5 * std::sqrt(5.0))
                                                           : 1 / (2 * std::sqrt(Physics::pi)));
        }

        void expectRelative(double actual, double expected, double tolerance = 2e-10) {
            EXPECT_NEAR(actual, expected, tolerance * std::abs(expected));
        }

        Matrix rotation() {
            // Orthogonal, with all three spatial directions mixed.
            return {{{0.36, 0.48, -0.8}, {-0.8, 0.6, 0.0}, {0.48, 0.64, 0.6}}};
        }

        Matrix rotate(const Matrix& matrix, const Matrix& axes) {
            Matrix result{};
            for (unsigned i = 0; i < 3; ++i)
                for (unsigned j = 0; j < 3; ++j)
                    for (unsigned a = 0; a < 3; ++a)
                        for (unsigned b = 0; b < 3; ++b)
                            result[i][j] += axes[i][a] * matrix[a][b] * axes[j][b];
            return result;
        }

        TEST(Analytic3DFieldTensor, UniformSphereUsesSqrtFiveRmsRadius) {
            const double sigma = 0.002, charge = 2e-9;
            const double radius = std::sqrt(5.0) * sigma;
            const auto result   = evaluateRestFrameTensor(
                    Profile::UniformEllipsoid,
                    diagonal(sigma * sigma, sigma * sigma, sigma * sigma), charge);
            const double expected = coulomb * charge / (radius * radius * radius);
            for (unsigned i = 0; i < 3; ++i) {
                expectRelative(result.gradient[i][i], expected);
                EXPECT_LE(result.absoluteError[i][i], 1e-10 * std::abs(result.gradient[i][i]));
                for (unsigned j = 0; j < 3; ++j)
                    if (i != j) EXPECT_DOUBLE_EQ(result.gradient[i][j], 0);
            }
            EXPECT_GT(result.integrandEvaluations, 0u);
            EXPECT_GT(result.integrationIntervals, 0u);
        }

        TEST(Analytic3DFieldTensor, GaussianRmsAgreesWithIndependentRadialForceAverage) {
            // Integrate the exact spherical Gaussian field against its own radial density:
            // <E_x*x>/sigma^2 = kQ/(3*sigma^3) * integral sqrt(2/pi)*t*exp(-t^2/2)
            //                  *[erf(t/sqrt(2))-sqrt(2/pi)*t*exp(-t^2/2)] dt.
            // Composite Simpson on [0,12] has negligible Gaussian tail and <1e-12 relative
            // discretization error here, independently of the production ellipsoid integral.
            constexpr unsigned count = 8192;
            constexpr double end = 12, step = end / count;
            const double root  = std::sqrt(2 / Physics::pi);
            auto weightedField = [root](double t) {
                const double density = root * t * std::exp(-0.5 * t * t);
                return density * (std::erf(t / std::sqrt(2.0)) - density);
            };
            double integral = weightedField(0) + weightedField(end);
            for (unsigned n = 1; n < count; ++n)
                integral += (n % 2 == 0 ? 2 : 4) * weightedField(n * step);
            integral *= step / 3;
            const double sigma = 0.003, charge = -4e-10;
            const auto result = evaluateRestFrameTensor(
                    Profile::GaussianRms, diagonal(sigma * sigma, sigma * sigma, sigma * sigma),
                    charge);
            const double averaged = coulomb * charge * integral / (3 * sigma * sigma * sigma);
            const double core     = coulomb * charge * root / (3 * sigma * sigma * sigma);
            for (unsigned i = 0; i < 3; ++i) {
                expectRelative(result.gradient[i][i], averaged);
                expectRelative(core / result.gradient[i][i], 2 * std::sqrt(2.0));
            }
        }

        TEST(Analytic3DFieldTensor, TriaxialReferenceAndPoissonTrace) {
            // Independent scipy QUADPACK integral in s on [0,infinity), epsabs=2e-14,
            // epsrel=2e-13; estimated errors <=3.4e-14. Reproduction script is in
            // opalx-documents/prs/693/match-694/analytic3d-validation/.
            const Vector factors{0.60286907673166157, 0.28478048169283698, 0.11235044157550395};
            const double charge = 1e-12, sigma = 0.001;
            const auto covariance = diagonal(sigma * sigma, 4 * sigma * sigma, 16 * sigma * sigma);
            for (auto profile : profiles) {
                const auto result  = evaluateRestFrameTensor(profile, covariance, charge);
                const double scale = normalization(profile, charge, 8 * sigma * sigma * sigma);
                double trace       = 0;
                for (unsigned i = 0; i < 3; ++i) {
                    expectRelative(result.gradient[i][i], scale * factors[i]);
                    EXPECT_LE(result.absoluteError[i][i], 1e-10 * std::abs(result.gradient[i][i]));
                    trace += result.gradient[i][i];
                }
                // For a uniform ellipsoid this is rho/epsilon_0, from Gauss's law.
                expectRelative(trace, scale);
            }
        }

        TEST(Analytic3DFieldTensor, ProlateAndOblateAnalyticLimitsAcrossAspectRatios) {
            for (double aspect : {2.0, 10.0, 1000.0, 1e6}) {
                for (bool prolate : {false, true}) {
                    const double longitudinal = prolate ? aspect : 1 / aspect;
                    double longitudinalFactor, transverseFactor;
                    if (prolate) {
                        const double inverse = 1 / aspect;
                        const double e       = std::sqrt(1 - inverse * inverse);
                        longitudinalFactor   = inverse * inverse / (e * e * e)
                                             * (std::log(aspect) + std::log1p(e) - e);
                        transverseFactor = 0.5 * (1 - longitudinalFactor);
                    } else {
                        const double e = std::sqrt(aspect * aspect - 1);
                        // Avoid 1-N_z cancellation for a very flat ellipsoid.
                        transverseFactor   = ((1 + e * e) * std::atan(e) - e) / (2 * e * e * e);
                        longitudinalFactor = 1 - 2 * transverseFactor;
                    }
                    for (auto profile : profiles) {
                        const auto result = evaluateRestFrameTensor(
                                profile, diagonal(1, 1, longitudinal * longitudinal), 1e-12);
                        const double scale = normalization(profile, 1e-12, longitudinal);
                        expectRelative(result.gradient[0][0], scale * transverseFactor);
                        expectRelative(result.gradient[1][1], scale * transverseFactor);
                        expectRelative(result.gradient[2][2], scale * longitudinalFactor);
                    }
                }
            }
        }

        TEST(Analytic3DFieldTensor, RotationsAndAxisPermutationsIncludeLongitudinalCoupling) {
            const auto axes = rotation();
            for (auto profile : profiles) {
                const auto covariance = diagonal(1e-6, 4e-6, 16e-6);
                const auto original   = evaluateRestFrameTensor(profile, covariance, 1e-9);
                const auto rotated =
                        evaluateRestFrameTensor(profile, rotate(covariance, axes), 1e-9);
                const auto expected = rotate(original.gradient, axes);
                const double scale  = original.gradient[0][0];
                for (unsigned i = 0; i < 3; ++i)
                    for (unsigned j = 0; j < 3; ++j) {
                        EXPECT_NEAR(rotated.gradient[i][j], expected[i][j], 2e-10 * scale);
                        EXPECT_DOUBLE_EQ(rotated.gradient[i][j], rotated.gradient[j][i]);
                        EXPECT_GE(rotated.absoluteError[i][j], 0);
                    }
                EXPECT_GT(std::abs(rotated.gradient[0][2]), 0.01 * scale);
                EXPECT_GT(std::abs(rotated.gradient[1][2]), 0.01 * scale);
                const auto permuted =
                        evaluateRestFrameTensor(profile, diagonal(16e-6, 1e-6, 4e-6), 1e-9);
                expectRelative(permuted.gradient[0][0], original.gradient[2][2]);
                expectRelative(permuted.gradient[1][1], original.gradient[0][0]);
                expectRelative(permuted.gradient[2][2], original.gradient[1][1]);
            }
        }

        TEST(Analytic3DFieldTensor, DegenerateAndNearlyEqualAxesAreContinuous) {
            for (auto profile : profiles) {
                for (double difference : {0.0, 1e-12, 1e-8}) {
                    const auto covariance = diagonal(1, 1 + difference, 4);
                    const auto original   = evaluateRestFrameTensor(profile, covariance, 1e-9);
                    const auto result =
                            evaluateRestFrameTensor(profile, rotate(covariance, rotation()), 1e-9);
                    const auto expected = rotate(original.gradient, rotation());
                    for (unsigned i = 0; i < 3; ++i)
                        for (unsigned j = 0; j < 3; ++j)
                            EXPECT_NEAR(
                                    result.gradient[i][j], expected[i][j],
                                    2e-10 * original.gradient[0][0]);
                }
                const auto sphere    = evaluateRestFrameTensor(profile, diagonal(1, 1, 1), 1e-9);
                const auto perturbed = evaluateRestFrameTensor(
                        profile, rotate(diagonal(1, 1 + 1e-10, 1 - 1e-10), rotation()), 1e-9);
                for (unsigned i = 0; i < 3; ++i)
                    for (unsigned j = 0; j < 3; ++j)
                        EXPECT_NEAR(
                                perturbed.gradient[i][j], sphere.gradient[i][j],
                                2e-10 * sphere.gradient[0][0]);
            }
        }

        TEST(Analytic3DFieldTensor, ChargeSizeAndProfileScaling) {
            for (auto profile : profiles) {
                const auto original =
                        evaluateRestFrameTensor(profile, diagonal(1e-6, 4e-6, 9e-6), 1e-9);
                const auto negative =
                        evaluateRestFrameTensor(profile, diagonal(1e-6, 4e-6, 9e-6), -3e-9);
                const auto larger =
                        evaluateRestFrameTensor(profile, diagonal(4e-6, 16e-6, 36e-6), 1e-9);
                for (unsigned i = 0; i < 3; ++i) {
                    expectRelative(negative.gradient[i][i], -3 * original.gradient[i][i]);
                    expectRelative(larger.gradient[i][i], original.gradient[i][i] / 8);
                    // Species charge has the same sign as total bunch charge: both signs defocus.
                    EXPECT_GT(-negative.gradient[i][i], 0);
                }
            }
            const auto uniform =
                    evaluateRestFrameTensor(Profile::UniformEllipsoid, diagonal(1, 4, 9), 1);
            const auto gaussian =
                    evaluateRestFrameTensor(Profile::GaussianRms, diagonal(1, 4, 9), 1);
            for (unsigned i = 0; i < 3; ++i)
                expectRelative(
                        gaussian.gradient[i][i] / uniform.gradient[i][i],
                        5 * std::sqrt(5.0) / (6 * std::sqrt(Physics::pi)));
        }

        TEST(Analytic3DFieldTensor, ExtremeSiScalesDoNotOverflowIntermediateVolumes) {
            for (auto profile : profiles) {
                const auto reference = evaluateRestFrameTensor(profile, diagonal(1, 4, 16), 1);
                for (double size : {1e-100, 1e100}) {
                    // Q and the RMS volume scale together, leaving the tensor unchanged.
                    const auto result = evaluateRestFrameTensor(
                            profile, diagonal(size * size, 4 * size * size, 16 * size * size),
                            size * size * size);
                    for (unsigned i = 0; i < 3; ++i)
                        expectRelative(result.gradient[i][i], reference.gradient[i][i], 1e-11);
                }
            }
        }

        TEST(Analytic3DFieldTensor, QuadratureRefinementAndErrorEstimates) {
            IntegrationSettings loose, tight;
            loose.relativeTolerance = 1e-4;
            tight.relativeTolerance = 1e-12;
            for (auto profile : profiles) {
                const auto a = evaluateRestFrameTensor(profile, diagonal(1, 4, 16), 1, loose);
                const auto b = evaluateRestFrameTensor(profile, diagonal(1, 4, 16), 1, tight);
                EXPECT_GT(b.integrandEvaluations, a.integrandEvaluations);
                for (unsigned i = 0; i < 3; ++i) {
                    EXPECT_LT(b.absoluteError[i][i], a.absoluteError[i][i]);
                    EXPECT_LE(b.absoluteError[i][i], 1e-12 * std::abs(b.gradient[i][i]));
                    EXPECT_NEAR(
                            a.gradient[i][i], b.gradient[i][i],
                            a.absoluteError[i][i] + b.absoluteError[i][i]);
                }
            }
        }

        TEST(Analytic3DFieldTensor, ExistingLorentzConversionGivesSymmetricReferenceForce) {
            const auto labCovariance = rotate(diagonal(1e-6, 4e-6, 9e-6), rotation());
            for (auto profile : profiles) {
                for (double gamma : {1.0, 2.0, 100.0}) {
                    // Local reference velocity is along +z; D=diag(1,1,gamma).
                    const Vector d{1, 1, gamma};
                    Matrix restCovariance{};
                    for (unsigned i = 0; i < 3; ++i)
                        for (unsigned j = 0; j < 3; ++j)
                            restCovariance[i][j] = d[i] * labCovariance[i][j] * d[j];
                    const auto tensor = evaluateRestFrameTensor(profile, restCovariance, 1e-9);
                    FieldCompositionPolicy policy;
                    policy.gamma        = gamma;
                    policy.meanMomentum = {0, 0, std::sqrt((gamma - 1) * (gamma + 1))};
                    const auto boost    = LorentzBoost::from(policy);
                    for (unsigned column = 0; column < 3; ++column) {
                        LorentzBoost::Vector electricRest(0), electricLab, magneticLab;
                        for (unsigned row = 0; row < 3; ++row)
                            electricRest[row] = tensor.gradient[row][column] * d[column];
                        boost.toLab(electricRest, electricLab, magneticLab);
                        const auto magneticForce = cross(boost.velocity, magneticLab);
                        const auto force         = electricLab + magneticForce;
                        for (unsigned row = 0; row < 3; ++row) {
                            const double expected =
                                    d[row] * tensor.gradient[row][column] * d[column] / gamma;
                            // Budget roundoff in BOTH cancellations: gamma*E' minus
                            // (gamma-1)*E'_parallel, and E + v x B. Scaling only by
                            // the final E_z would miss the first cancellation at large gamma.
                            const double arithmeticScale =
                                    gamma * std::abs(electricRest[row]) + std::abs(electricLab[row])
                                    + std::abs(magneticForce[row]) + std::abs(expected);
                            EXPECT_NEAR(
                                    force[row], expected,
                                    16 * std::numeric_limits<double>::epsilon() * arithmeticScale);
                        }
                    }
                }
            }
        }

        TEST(Analytic3DFieldTensor, ZeroChargeIsExactButStillValidatesInputs) {
            const auto result = evaluateRestFrameTensor(Profile::GaussianRms, diagonal(1, 4, 9), 0);
            EXPECT_EQ(result.integrandEvaluations, 0u);
            EXPECT_EQ(result.integrationIntervals, 0u);
            for (unsigned i = 0; i < 3; ++i)
                for (unsigned j = 0; j < 3; ++j) {
                    EXPECT_DOUBLE_EQ(result.gradient[i][j], 0);
                    EXPECT_DOUBLE_EQ(result.absoluteError[i][j], 0);
                }
            EXPECT_THROW(
                    (void)evaluateRestFrameTensor(Profile::GaussianRms, diagonal(1, 0, 1), 0),
                    OpalException);
        }

        TEST(Analytic3DFieldTensor, InvalidInputsAndNumericalFailuresAreReported) {
            const double infinity = std::numeric_limits<double>::infinity();
            const double nan      = std::numeric_limits<double>::quiet_NaN();
            const auto covariance = diagonal(1, 4, 9);
            for (auto profile : profiles) {
                for (double charge : {nan, infinity, -infinity})
                    EXPECT_THROW(
                            (void)evaluateRestFrameTensor(profile, covariance, charge),
                            OpalException);
                for (double diagonalValue : {0.0, -1.0, nan, infinity})
                    EXPECT_THROW(
                            (void)evaluateRestFrameTensor(
                                    profile, diagonal(1, diagonalValue, 1), 1),
                            OpalException);
                for (double offDiagonal : {1.0, 2.0, nan, infinity}) {
                    auto invalid  = diagonal(1, 1, 1);
                    invalid[0][1] = invalid[1][0] = offDiagonal;
                    EXPECT_THROW((void)evaluateRestFrameTensor(profile, invalid, 1), OpalException);
                }
                auto asymmetric  = covariance;
                asymmetric[0][1] = 0.1;
                EXPECT_THROW((void)evaluateRestFrameTensor(profile, asymmetric, 1), OpalException);
                auto roundoff  = covariance;
                roundoff[0][1] = 0.1;
                roundoff[1][0] = 0.1 + 2 * std::numeric_limits<double>::epsilon();
                EXPECT_NO_THROW((void)evaluateRestFrameTensor(profile, roundoff, 1));
                for (double tolerance : {0.0, -1.0, 1e-15, 0.01, nan, infinity}) {
                    IntegrationSettings settings;
                    settings.relativeTolerance = tolerance;
                    EXPECT_THROW(
                            (void)evaluateRestFrameTensor(profile, covariance, 1, settings),
                            OpalException);
                }
                IntegrationSettings exhausted;
                exhausted.maximumIntervals = 1;
                EXPECT_THROW(
                        (void)evaluateRestFrameTensor(profile, covariance, 1, exhausted),
                        OpalException);
                exhausted.maximumIntervals  = 2;
                exhausted.relativeTolerance = 1e-13;
                EXPECT_THROW(
                        (void)evaluateRestFrameTensor(profile, diagonal(1, 1, 1), 1, exhausted),
                        OpalException);
                exhausted.maximumIntervals = 0;
                EXPECT_THROW(
                        (void)evaluateRestFrameTensor(profile, covariance, 1, exhausted),
                        OpalException);
                EXPECT_THROW(
                        (void)evaluateRestFrameTensor(profile, diagonal(1e-300, 1e-300, 1e-300), 1),
                        OpalException);
                EXPECT_THROW(
                        (void)evaluateRestFrameTensor(
                                profile, diagonal(1e300, 1e300, 1e300), 1e-300),
                        OpalException);
            }
            EXPECT_THROW(
                    (void)evaluateRestFrameTensor(static_cast<Profile>(99), covariance, 0),
                    OpalException);
        }
    }  // namespace
}  // namespace opalx::spacecharge::analytic3d
