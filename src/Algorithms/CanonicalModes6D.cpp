#include "Algorithms/CanonicalModes6D.h"
// Match the existing real-LAPACKE C++ interface on host, CUDA and HIP compilers.
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
#include <algorithm>
#include <cmath>
#include <limits>
#include "Utilities/OpalException.h"

namespace opalx::maps::modes6d {
    namespace {
        using Dense          = std::array<double, 36>;
        using LongMatrix     = std::array<std::array<long double, 6>, 6>;
        constexpr unsigned n = 6;
        void require(bool condition, const std::string& message) {
            if (!condition) throw OpalException("CanonicalModes6D", message);
        }
        void checkInfo(lapack_int info, const char* routine) {
            require(info == 0, std::string(routine) + " failed with INFO=" + std::to_string(info));
        }
        double checked(long double x) {
            const double value = static_cast<double>(x);
            require(std::isfinite(value) && (x == 0 || value != 0),
                    "Numerical result is not representable as double.");
            return value;
        }
        void validate(const Settings& s) {
            for (double d : s.scales)
                require(std::isfinite(d) && d > 0, "Invalid coordinate scale.");
            for (double t :
                 {s.modulusTolerance, s.eigenResidualTolerance, s.symplecticTolerance,
                  s.covarianceTolerance, s.minimumBasisRcond, s.minimumRelativeSymplecticNorm})
                require(std::isfinite(t) && t > 0 && t < 1, "Invalid numerical tolerance.");
            require(std::isfinite(s.minimumPhaseSeparation) && s.minimumPhaseSeparation > 0
                            && s.minimumPhaseSeparation < 0.25,
                    "Invalid phase separation.");
            require(std::isfinite(s.minimumSubspaceOverlap) && s.minimumSubspaceOverlap > 0
                            && s.minimumSubspaceOverlap <= 1
                            && std::isfinite(s.minimumAssignmentGap) && s.minimumAssignmentGap > 0
                            && s.minimumAssignmentGap < 3,
                    "Invalid mode assignment policy.");
        }
        void finite(const Matrix& a) {
            for (const auto& row : a)
                for (double x : row)
                    require(std::isfinite(x), "Nonfinite matrix.");
        }
        struct Scaling {
            long double action = 0;
            std::array<double, 3> weight{};
            explicit Scaling(const Settings& s) {
                for (unsigned k = 0; k < 3; ++k)
                    action = std::max(
                            action,
                            static_cast<long double>(s.scales[2 * k]) * s.scales[2 * k + 1]);
                for (unsigned k = 0; k < 3; ++k)
                    weight[k] =
                            checked(static_cast<long double>(s.scales[2 * k]) * s.scales[2 * k + 1]
                                    / action);
            }
        };
        // Columns a,b in numerical coordinates. G=D^T J D / max_k(D_q D_p).
        long double area(const Dense& a, unsigned i, unsigned j, const Scaling& scale) {
            if (i == j) return 0;  // Exact skew diagonal; avoid a contracted x*y-y*x residual.
            long double value = 0;
            for (unsigned k = 0; k < 3; ++k)
                value +=
                        scale.weight[k]
                        * (static_cast<long double>(a[n * (2 * k) + i]) * a[n * (2 * k + 1) + j]
                           - static_cast<long double>(a[n * (2 * k + 1) + i]) * a[n * (2 * k) + j]);
            return value;
        }
        double canonicalResidual(const Matrix& basis) {
            finite(basis);
            long double sum = 0;
            for (unsigned i = 0; i < n; ++i)
                for (unsigned j = i + 1; j < n; ++j) {
                    long double value = 0;
                    for (unsigned k = 0; k < 3; ++k)
                        value += static_cast<long double>(basis[2 * k][i]) * basis[2 * k + 1][j]
                                 - static_cast<long double>(basis[2 * k + 1][i]) * basis[2 * k][j];
                    if (i % 2 == 0 && j == i + 1) value -= 1;
                    sum += 2 * value * value;
                }
            return checked(std::sqrt(sum / 6));
        }
        Matrix symmetricSPD(const Matrix& input) {
            finite(input);
            Matrix result = input;
            Dense correlation{};
            for (unsigned i = 0; i < n; ++i)
                require(input[i][i] > 0, "Covariance diagonal must be positive.");
            for (unsigned i = 0; i < n; ++i)
                for (unsigned j = 0; j < n; ++j) {
                    const long double d = std::sqrt(static_cast<long double>(input[i][i]))
                                          * std::sqrt(static_cast<long double>(input[j][j]));
                    require(std::abs(static_cast<long double>(input[i][j]) - input[j][i]) / d
                                    <= 64 * std::numeric_limits<double>::epsilon(),
                            "Covariance must be symmetric.");
                    result[i][j] =
                            checked((static_cast<long double>(input[i][j]) + input[j][i]) / 2);
                    correlation[n * i + j] = checked(result[i][j] / d);
                }
            checkInfo(
                    LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'L', n, correlation.data(), n),
                    "Covariance Cholesky (SPD required)");
            return result;
        }
        double covarianceResidual(const Matrix& a, const Matrix& reference) {
            long double sum = 0, norm = 0;
            for (unsigned i = 0; i < n; ++i)
                for (unsigned j = 0; j < n; ++j) {
                    const long double d = std::sqrt(static_cast<long double>(reference[i][i]))
                                          * std::sqrt(static_cast<long double>(reference[j][j]));
                    const long double r = (static_cast<long double>(a[i][j]) - reference[i][j]) / d;
                    const long double x = reference[i][j] / d;
                    sum += r * r;
                    norm += x * x;
                }
            return checked(std::sqrt(sum / norm));
        }
        Matrix covariance(const Matrix& basis, const Emittances& eps, const Settings& s) {
            require(canonicalResidual(basis) <= s.symplecticTolerance, "Basis is not canonical.");
            for (double e : eps)
                require(std::isfinite(e) && e > 0, "Emittances must be positive and finite.");
            LongMatrix weighted{};
            for (unsigned i = 0; i < n; ++i)
                for (unsigned k = 0; k < n; ++k)
                    weighted[i][k] = basis[i][k] * std::sqrt(static_cast<long double>(eps[k / 2]));
            Matrix result{};
            for (unsigned i = 0; i < n; ++i)
                for (unsigned j = i; j < n; ++j) {
                    long double value = 0;
                    for (unsigned k = 0; k < n; ++k)
                        value += weighted[i][k] * weighted[j][k];
                    result[j][i] = result[i][j] = checked(value);
                }
            return symmetricSPD(result);
        }
        Matrix transported(const Matrix& m, const Matrix& sigma) {
            LongMatrix left{};
            Matrix result{};
            for (unsigned i = 0; i < n; ++i)
                for (unsigned j = 0; j < n; ++j)
                    for (unsigned k = 0; k < n; ++k)
                        left[i][j] += static_cast<long double>(m[i][k]) * sigma[k][j];
            for (unsigned i = 0; i < n; ++i)
                for (unsigned j = i; j < n; ++j) {
                    long double value = 0;
                    for (unsigned k = 0; k < n; ++k)
                        value += left[i][k] * m[j][k];
                    result[j][i] = result[i][j] = checked(value);
                }
            return result;
        }
        void stable(const Analysis& a, const Settings& s) {
            require(a.status == Status::Stable && a.modes.size() == 3,
                    "Three resolved stable modes required.");
            finite(a.map);
            require(canonicalResidual(a.basis) <= s.symplecticTolerance,
                    "Invalid canonical mode basis.");
            // Check the public value snapshot before consuming it, including its branch signs.
            for (unsigned k = 0; k < 3; ++k) {
                const auto& mode = a.modes[k];
                require((mode.kreinSign == 1 || mode.kreinSign == -1) && std::isfinite(mode.phase)
                                && mode.phase > 0 && mode.phase < 0.5
                                && std::isfinite(mode.eigenvalue.real())
                                && std::isfinite(mode.eigenvalue.imag())
                                && mode.eigenvalue.imag() > 0
                                && std::abs(
                                           mode.phase
                                           - std::arg(mode.eigenvalue) / (2 * std::acos(-1.0)))
                                           <= s.eigenResidualTolerance
                                && std::abs(std::abs(mode.eigenvalue) - 1) <= s.modulusTolerance,
                        "Invalid mode metadata.");
                long double error = 0, norm = 0;
                for (unsigned i = 0; i < n; ++i) {
                    std::complex<long double> mv{};
                    for (unsigned j = 0; j < n; ++j)
                        mv += static_cast<long double>(a.map[i][j])
                              * std::complex<long double>(
                                      a.basis[j][2 * k], mode.kreinSign * a.basis[j][2 * k + 1]);
                    const std::complex<long double> v(
                            a.basis[i][2 * k], mode.kreinSign * a.basis[i][2 * k + 1]);
                    const auto lv       = std::complex<long double>(mode.eigenvalue) * v;
                    const auto mvScaled = mv / static_cast<long double>(s.scales[i]);
                    const auto lvScaled = lv / static_cast<long double>(s.scales[i]);
                    error += std::norm(mvScaled - lvScaled);
                    norm += std::norm(mvScaled) + std::norm(lvScaled);
                }
                require(norm > 0 && std::isfinite(norm) && std::isfinite(error)
                                && std::sqrt(error / norm) <= s.symplecticTolerance,
                        "Mode basis does not represent the supplied map.");
            }
        }
    }  // namespace

    Analysis analyze(const Matrix& input, const Settings& s) {
        validate(s);
        finite(input);
        Scaling scale(s);
        Dense original{}, a{}, vr{};
        double matrixNorm = 0;
        for (unsigned i = 0; i < n; ++i)
            for (unsigned j = 0; j < n; ++j) {
                original[n * i + j] =
                        checked(static_cast<long double>(input[i][j]) * s.scales[j] / s.scales[i]);
                matrixNorm = std::hypot(matrixNorm, original[n * i + j]);
            }
        require(std::isfinite(matrixNorm), "Scaled matrix norm overflow.");
        Analysis result;
        result.map       = input;
        long double symp = 0, metricNorm = 0;
        for (double w : scale.weight)
            metricNorm += 2 * static_cast<long double>(w) * w;
        for (unsigned i = 0; i < n; ++i)
            for (unsigned j = 0; j < n; ++j) {
                auto r = area(original, i, j, scale);
                if (i % 2 == 0 && j == i + 1) r -= scale.weight[i / 2];
                if (j % 2 == 0 && i == j + 1) r += scale.weight[j / 2];
                symp += r * r;
            }
        result.symplecticResidual = checked(std::sqrt(symp / metricNorm));
        a                         = original;
        std::array<double, n> wr{}, wi{};
        double unused = 0;
        checkInfo(
                LAPACKE_dgeev(
                        LAPACK_ROW_MAJOR, 'N', 'V', n, a.data(), n, wr.data(), wi.data(), &unused,
                        1, vr.data(), n),
                "LAPACKE_dgeev");
        bool growing = false, offCircle = false, marginal = false;
        std::vector<unsigned> pairs;
        for (unsigned j = 0; j < n; ++j) {
            const std::complex<double> lambda(wr[j], wi[j]);
            require(std::isfinite(wr[j]) && std::isfinite(wi[j]), "Nonfinite eigenvalue.");
            result.eigenvalues[j] = lambda;
            growing |= std::abs(lambda) > 1 + s.modulusTolerance;
            offCircle |= std::abs(std::abs(lambda) - 1) > s.modulusTolerance;
            const double phase = std::abs(std::arg(lambda)) / (2 * std::acos(-1.0));
            marginal |= wi[j] == 0 || phase <= s.minimumPhaseSeparation
                        || 0.5 - phase <= s.minimumPhaseSeparation;
            std::array<std::complex<double>, n> v{};
            double vn = 0;
            for (unsigned i = 0; i < n; ++i) {
                if (wi[j] > 0) {
                    require(j + 1 < n, "Malformed conjugate pair.");
                    v[i] = {vr[n * i + j], vr[n * i + j + 1]};
                } else if (wi[j] < 0) {
                    require(j > 0, "Malformed conjugate pair.");
                    v[i] = {vr[n * i + j - 1], -vr[n * i + j]};
                } else
                    v[i] = vr[n * i + j];
                vn = std::hypot(vn, std::abs(v[i]));
            }
            require(vn > 0 && std::isfinite(vn), "Invalid eigenvector norm.");
            double error = 0;
            for (unsigned i = 0; i < n; ++i) {
                auto r = -lambda * v[i] / vn;
                for (unsigned k = 0; k < n; ++k)
                    r += original[n * i + k] * (v[k] / vn);
                error = std::hypot(error, std::abs(r));
            }
            const double denominator = std::max(matrixNorm, std::abs(lambda));
            if (denominator > 0) error /= denominator;
            require(std::isfinite(error) && error <= s.eigenResidualTolerance,
                    "Eigenpair residual exceeds tolerance.");
            result.maximumEigenResidual = std::max(result.maximumEigenResidual, error);
            if (wi[j] > 0) pairs.push_back(j);
        }
        a = vr;
        std::array<double, n> singular{};
        std::array<double, n - 1> superb{};
        checkInfo(
                LAPACKE_dgesvd(
                        LAPACK_ROW_MAJOR, 'N', 'N', n, n, a.data(), n, singular.data(), &unused, 1,
                        &unused, 1, superb.data()),
                "LAPACKE_dgesvd");
        result.basisReciprocalCondition = singular[0] > 0 ? singular[n - 1] / singular[0] : 0;
        require(std::isfinite(result.basisReciprocalCondition), "Nonfinite eigenbasis condition.");
        if (result.symplecticResidual > s.symplecticTolerance)
            result.status = Status::NonSymplectic;
        else if (growing)
            result.status = Status::Unstable;
        else if (offCircle)
            result.status = Status::NonUnitCircle;
        else if (marginal || pairs.size() != 3)
            result.status = Status::Marginal;
        else if (result.basisReciprocalCondition < s.minimumBasisRcond)
            result.status = Status::IllConditioned;
        else {
            std::sort(pairs.begin(), pairs.end(), [&](unsigned i, unsigned j) {
                return std::arg(result.eigenvalues[i]) < std::arg(result.eigenvalues[j]);
            });
            for (unsigned k = 1; k < 3; ++k)
                if (std::abs(
                            std::arg(result.eigenvalues[pairs[k]])
                            - std::arg(result.eigenvalues[pairs[k - 1]]))
                            / (2 * std::acos(-1.0))
                    <= s.minimumPhaseSeparation) {
                    result.status = Status::Degenerate;
                    return result;
                }
            for (unsigned k = 0; k < 3; ++k) {
                const unsigned j             = pairs[k];
                const long double signedArea = area(vr, j, j + 1, scale);
                long double absolute         = 0;
                for (unsigned p = 0; p < 3; ++p)
                    absolute += scale.weight[p]
                                * (std::abs(
                                           static_cast<long double>(vr[n * (2 * p) + j])
                                           * vr[n * (2 * p + 1) + j + 1])
                                   + std::abs(
                                           static_cast<long double>(vr[n * (2 * p + 1) + j])
                                           * vr[n * (2 * p) + j + 1]));
                if (absolute == 0
                    || std::abs(signedArea) / absolute < s.minimumRelativeSymplecticNorm) {
                    result.status = Status::IllConditioned;
                    result.basis  = {};
                    result.modes.clear();
                    return result;
                }
                const int sign                = signedArea > 0 ? 1 : -1;
                const long double denominator = std::sqrt(scale.action * std::abs(signedArea));
                for (unsigned i = 0; i < n; ++i) {
                    result.basis[i][2 * k] = checked(
                            s.scales[i] * static_cast<long double>(vr[n * i + j]) / denominator);
                    result.basis[i][2 * k + 1] =
                            checked(sign * s.scales[i] * static_cast<long double>(vr[n * i + j + 1])
                                    / denominator);
                }
                result.modes.push_back(
                        {result.eigenvalues[j],
                         std::arg(result.eigenvalues[j]) / (2 * std::acos(-1.0)), sign});
            }
            result.normalizationResidual = canonicalResidual(result.basis);
            result.status                = result.normalizationResidual <= s.symplecticTolerance
                                                   ? Status::Stable
                                                   : Status::IllConditioned;
        }
        if (result.status != Status::Stable) {
            result.basis = {};
            result.modes.clear();
        }
        return result;
    }

    Analysis analyze(const tps::Map<double>& map, const Settings& s) {
        const auto& context = map.context();
        const auto& chart   = context.input();
        require(context.order() == 1 && chart.coordinates.size() == n && chart == context.output()
                        && chart.convention == tps::CoordinateConvention::canonical
                        && chart.canonicalPairs
                                   == std::vector<std::pair<
                                           std::size_t, std::size_t>>{{0, 1}, {2, 3}, {4, 5}},
                "Requires a first-order same-chart interleaved canonical 6D TPS map.");
        for (unsigned i = 0; i < n; ++i) {
            require(chart.coordinates[i].unit == (i % 2 ? "p/(mc)" : "m"),
                    "Invalid canonical coordinate units.");
            const auto& component = map.polynomial()[i];
            for (int j = 0; j < component.getSize(); ++j)
                require(std::isfinite(component[j]), "Nonfinite TPS coefficient.");
            require(component[0] == 0, "Affine map requires an explicit reference recentering.");
        }
        const auto jacobian = map.jacobian(std::array<double, n>{});
        Matrix physical{};
        for (unsigned i = 0; i < n; ++i)
            for (unsigned j = 0; j < n; ++j)
                physical[i][j] =
                        checked(static_cast<long double>(jacobian(i, j))
                                * chart.coordinates[i].scale / chart.coordinates[j].scale);
        return analyze(physical, s);
    }

    Matrix invariantCovariance(const Analysis& a, const Emittances& eps, const Settings& s) {
        validate(s);
        stable(a, s);
        auto result = covariance(a.basis, eps, s);
        require(covarianceResidual(transported(a.map, result), result) <= s.covarianceTolerance,
                "Constructed covariance is not periodic within tolerance.");
        return result;
    }

    Williamson williamson(const Matrix& input, const Settings& s) {
        validate(s);
        const auto sigma = symmetricSPD(input);
        const Scaling scale(s);
        long double covarianceScale = 0;
        for (unsigned i = 0; i < n; ++i)
            covarianceScale = std::max(
                    covarianceScale,
                    sigma[i][i] / (static_cast<long double>(s.scales[i]) * s.scales[i]));
        Dense l{};
        for (unsigned i = 0; i < n; ++i)
            for (unsigned j = 0; j < n; ++j)
                l[n * i + j] = checked(
                        sigma[i][j]
                        / (static_cast<long double>(s.scales[i]) * s.scales[j] * covarianceScale));
        checkInfo(LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'L', n, l.data(), n), "LAPACKE_dpotrf");
        for (unsigned i = 0; i < n; ++i)
            for (unsigned j = i + 1; j < n; ++j)
                l[n * i + j] = 0;
        Dense skew{}, q{};
        for (unsigned i = 0; i < n; ++i)
            for (unsigned j = i + 1; j < n; ++j) {
                skew[n * i + j] = checked(area(l, i, j, scale));
                skew[n * j + i] = -skew[n * i + j];
            }
        std::array<double, n> wr{}, wi{};
        lapack_int selected = 0;
        checkInfo(
                LAPACKE_dgees(
                        LAPACK_ROW_MAJOR, 'V', 'N', nullptr, n, skew.data(), n, &selected,
                        wr.data(), wi.data(), q.data(), n),
                "LAPACKE_dgees");
        std::array<unsigned, 3> pairs{};
        unsigned count = 0;
        for (unsigned j = 0; j < n;) {
            require(j + 1 < n && wi[j] > 0 && wi[j + 1] < 0 && std::isfinite(wi[j]),
                    "Unresolved positive symplectic eigenemittance.");
            pairs[count++] = j;
            j += 2;
        }
        require(count == 3, "Three Williamson pairs required.");
        std::sort(pairs.begin(), pairs.end(), [&](unsigned i, unsigned j) {
            return wi[i] < wi[j];
        });
        Williamson result;
        for (unsigned k = 0; k < 3; ++k) {
            unsigned j               = pairs[k];
            const long double action = scale.action * wi[j];
            result.emittances[k]     = checked(covarianceScale * action);
            const int orientation    = skew[n * j + j + 1] > 0 ? 1 : -1;
            for (unsigned i = 0; i < n; ++i)
                for (unsigned side = 0; side < 2; ++side) {
                    long double value = 0;
                    for (unsigned r = 0; r < n; ++r)
                        value += static_cast<long double>(l[n * i + r]) * q[n * r + j + side];
                    result.basis[i][2 * k + side] = checked(
                            s.scales[i] * value / std::sqrt(action) * (side ? orientation : 1));
                }
        }
        result.normalizationResidual = canonicalResidual(result.basis);
        require(result.normalizationResidual <= s.symplecticTolerance,
                "Williamson basis normalization failed.");
        const auto reconstructed      = covariance(result.basis, result.emittances, s);
        result.reconstructionResidual = covarianceResidual(reconstructed, sigma);
        require(result.reconstructionResidual <= s.covarianceTolerance,
                "Williamson covariance reconstruction failed.");
        return result;
    }

    Matrix withEmittances(const Williamson& w, const Emittances& target, const Settings& s) {
        validate(s);
        return covariance(w.basis, target, s);
    }

    Analysis trackModes(const Analysis& current, const Analysis& reference, const Settings& s) {
        validate(s);
        stable(current, s);
        stable(reference, s);
        // Two-pass modified Gram-Schmidt makes the overlap independent of eigenvector phase
        // and normalization. Both bases use the SAME numerical coordinate metric.
        auto orthonormal = [&](const Matrix& basis) {
            LongMatrix q{};
            for (unsigned k = 0; k < 3; ++k)
                for (unsigned side = 0; side < 2; ++side) {
                    const unsigned j = 2 * k + side;
                    for (unsigned i = 0; i < n; ++i)
                        q[i][j] = basis[i][j] / static_cast<long double>(s.scales[i]);
                    for (unsigned pass = 0; side && pass < 2; ++pass) {
                        long double dot = 0;
                        for (unsigned i = 0; i < n; ++i)
                            dot += q[i][j] * q[i][j - 1];
                        for (unsigned i = 0; i < n; ++i)
                            q[i][j] -= dot * q[i][j - 1];
                    }
                    long double norm = 0;
                    for (unsigned i = 0; i < n; ++i)
                        norm = std::hypot(norm, q[i][j]);
                    require(norm > 0 && std::isfinite(norm), "Unresolved mode subspace.");
                    for (unsigned i = 0; i < n; ++i)
                        q[i][j] /= norm;
                }
            return q;
        };
        const auto a = orthonormal(reference.basis), b = orthonormal(current.basis);
        std::array<std::array<double, 3>, 3> overlap{};
        for (unsigned k = 0; k < 3; ++k)
            for (unsigned l = 0; l < 3; ++l)
                for (unsigned i = 0; i < 2; ++i)
                    for (unsigned j = 0; j < 2; ++j) {
                        long double dot = 0;
                        for (unsigned r = 0; r < n; ++r)
                            dot += a[r][2 * k + i] * b[r][2 * l + j];
                        overlap[k][l] += static_cast<double>(dot * dot / 2);
                    }
        std::array<unsigned, 3> permutation{0, 1, 2}, best{};
        double first = -1, second = -1;
        do {
            double score = 0;
            for (unsigned k = 0; k < 3; ++k)
                score += overlap[k][permutation[k]];
            if (score > first) {
                second = first;
                first  = score;
                best   = permutation;
            } else
                second = std::max(second, score);
        } while (std::next_permutation(permutation.begin(), permutation.end()));
        require(first - second >= s.minimumAssignmentGap, "Ambiguous mode subspace assignment.");
        Analysis result = current;
        for (unsigned k = 0; k < 3; ++k) {
            const unsigned l = best[k];
            require(overlap[k][l] >= s.minimumSubspaceOverlap,
                    "Insufficient mode subspace overlap.");
            require(current.modes[l].kreinSign == reference.modes[k].kreinSign,
                    "Mode Krein sign changed.");
            result.modes[k] = current.modes[l];
            for (unsigned i = 0; i < n; ++i)
                for (unsigned side = 0; side < 2; ++side)
                    result.basis[i][2 * k + side] = current.basis[i][2 * l + side];
        }
        return result;
    }
}  // namespace opalx::maps::modes6d
