#include "Algorithms/SpaceChargeMap.h"

#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <utility>

#include "Physics/Physics.h"
#include "Utilities/OpalException.h"
#include "tps/Algebra/Matrix.h"

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

namespace opalx::maps {
    namespace {
        using Matrix        = SpaceChargeMatrix;
        using Vector        = opalx::spacecharge::analytic3d::Vector;
        using SpatialMatrix = opalx::spacecharge::analytic3d::Matrix;
        using HostMatrix    = tps::Matrix<double>;
        using Chart         = SpaceChargeChart;

        void require(bool condition, const char* message) {
            if (!condition) throw OpalException("SpaceChargeMap", message);
        }
        unsigned index(Chart chart) {
            const auto value = static_cast<unsigned>(chart);
            require(value < 3, "Unknown space-charge coordinate chart.");
            return value;
        }
        Matrix identity() {
            Matrix result{};
            for (unsigned i = 0; i < 6; ++i)
                result[i][i] = 1;
            return result;
        }
        void finite(const Matrix& matrix) {
            for (const auto& row : matrix)
                for (double value : row)
                    require(std::isfinite(value), "Nonfinite 6D coordinate matrix.");
        }
        HostMatrix host(const Matrix& matrix) {
            HostMatrix result(6, 6);
            for (unsigned i = 0; i < 6; ++i)
                for (unsigned j = 0; j < 6; ++j)
                    result(i, j) = matrix[i][j];
            return result;
        }
        Matrix stored(const HostMatrix& matrix) {
            Matrix result{};
            for (unsigned i = 0; i < 6; ++i)
                for (unsigned j = 0; j < 6; ++j)
                    result[i][j] = matrix(i, j);
            finite(result);
            return result;
        }
        Matrix validatedCovariance(const Matrix& covariance) {
            Matrix result{};
            Kokkos::Array<double, 36> correlation{};
            for (unsigned i = 0; i < 6; ++i)
                require(std::isfinite(covariance[i][i]) && covariance[i][i] > 0,
                        "Covariance diagonal must be positive and finite.");
            for (unsigned i = 0; i < 6; ++i) {
                for (unsigned j = 0; j < 6; ++j) {
                    const double a = covariance[i][j], b = covariance[j][i];
                    const double scale = std::sqrt(covariance[i][i]) * std::sqrt(covariance[j][j]);
                    require(std::isfinite(a) && std::isfinite(b)
                                    && std::abs(a - b)
                                               <= 64 * std::numeric_limits<double>::epsilon()
                                                          * scale,
                            "Covariance must be finite and symmetric.");
                    result[i][j]           = 0.5 * a + 0.5 * b;
                    correlation[6 * i + j] = result[i][j] / std::sqrt(covariance[i][i])
                                             / std::sqrt(covariance[j][j]);
                    require(std::isfinite(correlation[6 * i + j]), "Covariance scaling overflow.");
                }
            }
            require(LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'L', 6, correlation.data(), 6) == 0,
                    "Covariance is not numerically positive definite.");
            return result;
        }
        Matrix transportCovariance(const HostMatrix& transform, const Matrix& covariance) {
            const auto input = validatedCovariance(covariance);
            return validatedCovariance(stored(transform * host(input) * transform.transpose()));
        }
        tps::Map<double> linearMap(
                const Matrix& matrix, const tps::PhaseSpaceChart& input,
                const tps::PhaseSpaceChart& output, const std::string& id) {
            finite(matrix);
            try {
                tps::MapContext context(input, output, 1, {id});
                tps::VpsMap<double> polynomial(6, 6);
                for (unsigned i = 0; i < 6; ++i) {
                    polynomial[i] = tps::Tps<double>(1, 6);
                    for (unsigned j = 0; j < 6; ++j)
                        polynomial[i][j + 1] = matrix[i][j];
                }
                return tps::Map<double>(std::move(context), std::move(polynomial));
            } catch (const std::invalid_argument& error) {
                throw OpalException("SpaceChargeMap", error.what());
            }
        }
        Vector rotateField(const Vector& field, const tps::ReferenceFrame& frame) {
            Vector result{};
            for (double value : field)
                require(std::isfinite(value), "Nonfinite external reference field.");
            for (unsigned i = 0; i < 3; ++i) {
                for (unsigned j = 0; j < 3; ++j)
                    result[i] += frame.axes[3 * i + j] * field[j];
                require(std::isfinite(result[i]), "External field rotation overflow.");
            }
            return result;
        }
        Vector crossProduct(const Vector& a, const Vector& b) {
            return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                    a[0] * b[1] - a[1] * b[0]};
        }
    }  // namespace

    SpaceChargeMatrix transportSpaceChargeCovariance(
            const tps::Map<double>& map, const Matrix& covariance) {
        require(map.context().order() == 1 && map.context().input().coordinates.size() == 6
                        && map.context().output().coordinates.size() == 6,
                "Covariance transport requires a first-order 6D map.");
        for (unsigned i = 0; i < 6; ++i)
            require(map.polynomial()[i][0] == 0,
                    "Covariance transport requires a zero affine term.");
        return transportCovariance(map.jacobian(std::array<double, 6>{}), covariance);
    }

    SpaceChargeCoordinates::SpaceChargeCoordinates(const SpaceChargeReference& reference)
        : massEV_m(reference.massEV), chargeE_m(reference.chargeE) {
        require(std::isfinite(massEV_m) && massEV_m > 0 && std::isfinite(chargeE_m),
                "Invalid reference species charge or rest energy.");
        charts_m[0]       = linearTransferMapChart(reference.orbit);
        const auto& frame = charts_m[0].frame;
        Vector u{};
        for (unsigned i = 0; i < 3; ++i)
            for (unsigned j = 0; j < 3; ++j)
                u[i] += frame.axes[3 * i + j] * reference.orbit.momentum(j);
        const double momentum = *charts_m[0].referenceMomentum;
        require(std::isfinite(u[2]) && u[2] / momentum > 1e-12,
                "Space-charge reference must positively cross a nongrazing section (n_s > 1e-12).");
        gamma_m            = std::hypot(1.0, momentum);
        beta_m             = momentum / gamma_m;
        const double kappa = chargeE_m / massEV_m * Physics::c;
        require(std::isfinite(kappa) && (chargeE_m == 0 || kappa != 0) && beta_m > 0,
                "Reference force normalization is not representable.");
        for (unsigned i = 0; i < 3; ++i)
            direction_m[i] = u[i] / momentum;
        const auto electric      = rotateField(reference.electricFieldLab, frame);
        const auto magnetic      = rotateField(reference.magneticFieldLab, frame);
        const auto magneticForce = crossProduct(direction_m, magnetic);
        Vector force{}, normalizedB{};
        for (unsigned i = 0; i < 3; ++i) {
            force[i]       = kappa * (electric[i] / (beta_m * Physics::c) + magneticForce[i]);
            normalizedB[i] = kappa * magnetic[i];
            require(std::isfinite(normalizedB[i]), "Canonical gauge normalization overflow.");
        }
        auto& toTime    = toMechanical_m[0];
        auto& toSection = fromMechanical_m[0];
        toTime[0][0] = toTime[2][2] = 1;
        toSection[0][0] = toSection[2][2] = 1;
        // Write the diagonal projections as sums of squares. 1-n_x^2 loses
        // the small positive derivative when the momentum is almost along x.
        const auto& n = direction_m;
        const Vector slopeX{n[1] * n[1] + n[2] * n[2], -n[0] * n[1], -n[0] * n[2]};
        const Vector slopeY{-n[0] * n[1], n[0] * n[0] + n[2] * n[2], -n[1] * n[2]};
        for (unsigned i = 0; i < 3; ++i) {
            toTime[2 * i][4]        = direction_m[i];
            toTime[2 * i + 1][1]    = u[2] * slopeX[i];
            toTime[2 * i + 1][3]    = u[2] * slopeY[i];
            toTime[2 * i + 1][4]    = force[i];
            toTime[2 * i + 1][5]    = u[i];
            toSection[1][2 * i + 1] = ((i == 0 ? 1.0 : 0.0) - (i == 2 ? u[0] / u[2] : 0.0)) / u[2];
            toSection[3][2 * i + 1] = ((i == 1 ? 1.0 : 0.0) - (i == 2 ? u[1] / u[2] : 0.0)) / u[2];
            toSection[5][2 * i + 1] = direction_m[i] / momentum;
        }
        toSection[0][4] = -u[0] / u[2];
        toSection[2][4] = -u[1] / u[2];
        toSection[4][4] = 1 / direction_m[2];
        for (unsigned row : {1u, 3u, 5u}) {
            double projectedForce = 0;
            for (unsigned i = 0; i < 3; ++i)
                projectedForce += toSection[row][2 * i + 1] * force[i];
            toSection[row][4] = -projectedForce / direction_m[2];
        }
        toMechanical_m[1] = fromMechanical_m[1] = identity();
        toMechanical_m[2] = fromMechanical_m[2] = identity();
        const SpatialMatrix crossB{
                {{0, -normalizedB[2], normalizedB[1]},
                 {normalizedB[2], 0, -normalizedB[0]},
                 {-normalizedB[1], normalizedB[0], 0}}};
        for (unsigned i = 0; i < 3; ++i)
            for (unsigned j = 0; j < 3; ++j) {
                fromMechanical_m[2][2 * i + 1][2 * j] = 0.5 * crossB[i][j];
                toMechanical_m[2][2 * i + 1][2 * j]   = -0.5 * crossB[i][j];
            }
        for (unsigned kind = 0; kind < 3; ++kind) {
            finite(toMechanical_m[kind]);
            finite(fromMechanical_m[kind]);
        }
        charts_m[1]                = charts_m[0];
        charts_m[1].name           = "common-time mechanical perturbations";
        charts_m[1].coordinates    = {{"x", "m", 1},       {"ux", "p/(mc)", 1}, {"y", "m", 1},
                                      {"uy", "p/(mc)", 1}, {"z", "m", 1},       {"uz", "p/(mc)", 1}};
        charts_m[1].reference      = {0, u[0], 0, u[1], 0, u[2]};
        charts_m[1].parameters     = {{"species_rest_energy", "eV", massEV_m}};
        charts_m[2]                = charts_m[1];
        charts_m[2].name           = "common-time canonical perturbations in local symmetric gauge";
        charts_m[2].convention     = tps::CoordinateConvention::canonical;
        charts_m[2].canonicalPairs = {{0, 1}, {2, 3}, {4, 5}};
        for (unsigned i = 0; i < 3; ++i) {
            charts_m[2].coordinates[2 * i + 1].name = std::string("pi_") + "xyz"[i];
            charts_m[2].parameters.push_back(
                    {std::string("kappa_B_") + "xyz"[i], "1/m", normalizedB[i]});
        }
    }

    const tps::PhaseSpaceChart& SpaceChargeCoordinates::chart(Chart kind) const {
        return charts_m[index(kind)];
    }

    SpaceChargeMatrix SpaceChargeCoordinates::jacobian(Chart from, Chart to) const {
        const auto input = index(from), output = index(to);
        if (input == output) return identity();
        if (from == Chart::CommonTimeMechanical) return fromMechanical_m[output];
        if (to == Chart::CommonTimeMechanical) return toMechanical_m[input];
        return stored(host(fromMechanical_m[output]) * host(toMechanical_m[input]));
    }

    tps::Map<double> SpaceChargeCoordinates::map(
            Chart from, Chart to, const std::string& id) const {
        return linearMap(jacobian(from, to), chart(from), chart(to), id);
    }

    SpaceChargeMatrix SpaceChargeCoordinates::convertCovariance(
            const Matrix& covariance, Chart from, Chart to) const {
        const auto transform = host(jacobian(from, to));
        return transportCovariance(transform, covariance);
    }

    SpaceChargeKick makeSpaceChargeKick(
            const SpaceChargeCoordinates& coordinates, const SpaceChargeMatrix& covariance,
            Chart chart, opalx::spacecharge::analytic3d::Profile profile, double totalCharge,
            double pathLength, const std::string& id,
            const opalx::spacecharge::analytic3d::IntegrationSettings& settings) {
        require(std::isfinite(pathLength), "SC path length must be finite metres.");
        const auto mechanical =
                coordinates.convertCovariance(covariance, chart, Chart::CommonTimeMechanical);
        SpatialMatrix stretch{}, inverseForceBoost{}, restCovariance{};
        const auto& n      = coordinates.direction();
        const double gamma = coordinates.gamma();
        for (unsigned i = 0; i < 3; ++i)
            for (unsigned j = 0; j < 3; ++j) {
                stretch[i][j] = (i == j ? 1.0 : 0.0) + (gamma - 1) * n[i] * n[j];
                inverseForceBoost[i][j] =
                        (i == j ? 1 / gamma : 0.0) + (1 - 1 / gamma) * n[i] * n[j];
            }
        for (unsigned i = 0; i < 3; ++i)
            for (unsigned j = i; j < 3; ++j) {
                double value = 0;
                for (unsigned a = 0; a < 3; ++a)
                    for (unsigned b = 0; b < 3; ++b)
                        value += stretch[i][a] * mechanical[2 * a][2 * b] * stretch[j][b];
                restCovariance[i][j] = restCovariance[j][i] = value;
            }
        const auto tensor = opalx::spacecharge::analytic3d::evaluateRestFrameTensor(
                profile, restCovariance, totalCharge, settings);
        SpatialMatrix gradient{}, error{};
        for (unsigned i = 0; i < 3; ++i)
            for (unsigned j = i; j < 3; ++j) {
                double value = 0, estimate = 0;
                for (unsigned a = 0; a < 3; ++a)
                    for (unsigned b = 0; b < 3; ++b) {
                        const double coefficient = 0.5
                                                   * (inverseForceBoost[i][a] * stretch[j][b]
                                                      + inverseForceBoost[j][a] * stretch[i][b]);
                        value += coefficient * tensor.gradient[a][b];
                        estimate += std::abs(coefficient) * tensor.absoluteError[a][b];
                    }
                require(std::isfinite(value) && std::isfinite(estimate),
                        "SC force tensor overflow.");
                gradient[i][j] = gradient[j][i] = value;
                error[i][j] = error[j][i] = estimate;
            }
        Matrix kick = identity();
        if (pathLength != 0 && coordinates.chargeE() != 0 && totalCharge != 0) {
            // Scale each final coefficient: h/beta can overflow even when the
            // final impulse is finite after multiplication by a small field.
            const double logStrength =
                    std::log(std::abs(pathLength)) + std::log(std::abs(coordinates.chargeE()))
                    - std::log(coordinates.beta()) - std::log(coordinates.massEV());
            const double sign =
                    std::copysign(1.0, pathLength) * std::copysign(1.0, coordinates.chargeE());
            for (unsigned i = 0; i < 3; ++i)
                for (unsigned j = 0; j < 3; ++j) {
                    if (gradient[i][j] == 0) continue;
                    const double magnitude =
                            std::exp(logStrength + std::log(std::abs(gradient[i][j])));
                    require(std::isfinite(magnitude) && magnitude > 0,
                            "SC kick coefficient overflows or underflows double precision.");
                    kick[2 * i + 1][2 * j] = sign * std::copysign(magnitude, gradient[i][j]);
                }
            finite(kick);
            if (chart == Chart::SectionSlopes) {
                const auto toSection =
                        host(coordinates.jacobian(Chart::CommonTimeMechanical, chart));
                const auto toTime = host(coordinates.jacobian(chart, Chart::CommonTimeMechanical));
                kick              = stored(toSection * host(kick) * toTime);
            }
            // At common time, r is unchanged: the symmetric-gauge momentum shift
            // commutes with this kick, so mechanical and canonical shears are identical.
        }
        return {linearMap(kick, coordinates.chart(chart), coordinates.chart(chart), id), tensor,
                restCovariance, gradient, error};
    }
}  // namespace opalx::maps
