// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#include "Algorithms/LinearMapEigenAnalysis.h"
// Only real LAPACKE routines are called. Use the documented C++ type overrides
// for unused complex declarations so host/device C++ compilers need not parse
// the C99 _Complex extension. Keep these macros local to this translation unit.
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include <lapacke.h>
#undef lapack_complex_float
#undef lapack_complex_double
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include "Utilities/OpalException.h"

namespace {
    void require(bool valid, const std::string& message) {
        if (!valid) throw OpalException("LinearMapEigenAnalysis", message);
    }
    void checkInfo(lapack_int info, const char* routine) {
        require(info == 0, std::string(routine) + " failed with INFO=" + std::to_string(info));
    }
}  // namespace

LinearMapEigenAnalysis::Result LinearMapEigenAnalysis::analyze(
        const Matrix& input, const Settings& s) {
    require(std::isfinite(s.modulusTolerance) && s.modulusTolerance > 0 && s.modulusTolerance < 1
                    && std::isfinite(s.residualTolerance) && s.residualTolerance > 0
                    && std::isfinite(s.nearIntegerTolerance) && s.nearIntegerTolerance > 0
                    && s.nearIntegerTolerance < 0.5 && std::isfinite(s.minimumBasisRcond)
                    && s.minimumBasisRcond > 0 && s.minimumBasisRcond < 1,
            "Invalid eigenanalysis tolerances.");
    for (double scale : s.scales)
        require(std::isfinite(scale) && scale > 0, "Invalid coordinate scale.");
    std::array<double, 16> original{}, a{}, vr{};
    double matrixNorm = 0;
    for (unsigned i = 0; i < 4; ++i)
        for (unsigned j = 0; j < 4; ++j) {
            require(std::isfinite(input[i][j]), "Nonfinite input matrix.");
            original[4 * i + j] = input[i][j] * s.scales[j] / s.scales[i];
            require(std::isfinite(original[4 * i + j]), "Nonfinite scaled matrix.");
            matrixNorm = std::hypot(matrixNorm, original[4 * i + j]);
        }
    require(std::isfinite(matrixNorm), "Scaled matrix norm overflow.");
    a = original;
    std::array<double, 4> wr{}, wi{};
    double unused = 0;
    checkInfo(
            LAPACKE_dgeev(
                    LAPACK_ROW_MAJOR, 'N', 'V', 4, a.data(), 4, wr.data(), wi.data(), &unused, 1,
                    vr.data(), 4),
            "LAPACKE_dgeev");
    Result result;
    bool growing = false, offCircle = false, realUnit = false;
    const double twoPi = 2 * std::acos(-1.0);
    for (unsigned j = 0; j < 4; ++j) {
        const Complex lambda(wr[j], wi[j]);
        require(std::isfinite(wr[j]) && std::isfinite(wi[j]), "Nonfinite eigenvalue.");
        result.eigenvalues[j] = lambda;
        const double modulus  = std::abs(lambda);
        const bool unit       = std::abs(modulus - 1) <= s.modulusTolerance;
        growing               = growing || modulus > 1 + s.modulusTolerance;
        offCircle             = offCircle || !unit;
        realUnit              = realUnit || (unit && wi[j] == 0);
        result.nearInteger =
                result.nearInteger
                || (unit && std::abs(std::arg(lambda)) / twoPi <= s.nearIntegerTolerance);
        auto& v   = result.eigenvectors[j];
        double vn = 0;
        for (unsigned i = 0; i < 4; ++i) {
            if (wi[j] > 0) {
                require(j + 1 < 4, "Malformed conjugate pair.");
                v[i] = Complex(vr[4 * i + j], vr[4 * i + j + 1]);
            } else if (wi[j] < 0) {
                require(j > 0, "Malformed conjugate pair.");
                v[i] = Complex(vr[4 * i + j - 1], -vr[4 * i + j]);
            } else
                v[i] = vr[4 * i + j];
            vn = std::hypot(vn, std::abs(v[i]));
        }
        require(std::isfinite(vn) && vn > 0, "Invalid eigenvector.");
        for (auto& value : v)
            value /= vn;
        double residual = 0;
        for (unsigned i = 0; i < 4; ++i) {
            Complex r = -lambda * v[i];
            for (unsigned k = 0; k < 4; ++k)
                r += original[4 * i + k] * v[k];
            residual = std::hypot(residual, std::abs(r));
        }
        const double denominator = std::max(matrixNorm, modulus);
        const double relative    = denominator > 0 ? residual / denominator : residual;
        require(std::isfinite(relative) && relative <= s.residualTolerance,
                "Eigenpair residual exceeds tolerance.");
        result.maximumResidual = std::max(result.maximumResidual, relative);
        if (wi[j] > 0) {
            Mode mode;
            mode.eigenvalue = lambda;
            if (unit) {
                mode.tune              = std::arg(lambda) / twoPi;
                mode.complementaryTune = 1 - *mode.tune;
                mode.nearInteger       = *mode.tune <= s.nearIntegerTolerance;
            }
            result.modes.push_back(mode);
        }
    }
    // Detect defective/nearly dependent eigenvector bases, including repeated
    // complex pairs away from +/-1. This is a diagnostic, not a proof of stability.
    std::array<double, 4> singular{};
    std::array<double, 3> superb{};
    checkInfo(
            LAPACKE_dgesvd(
                    LAPACK_ROW_MAJOR, 'N', 'N', 4, 4, vr.data(), 4, singular.data(), &unused, 1,
                    &unused, 1, superb.data()),
            "LAPACKE_dgesvd");
    result.basisReciprocalCondition = singular[0] > 0 ? singular[3] / singular[0] : 0;
    require(std::isfinite(result.basisReciprocalCondition),
            "Nonfinite eigenvector basis condition.");
    result.stability = growing     ? Stability::Unstable
                       : offCircle ? Stability::NonUnitCircle
                       : (realUnit || result.nearInteger
                          || result.basisReciprocalCondition < s.minimumBasisRcond)
                               ? Stability::Marginal
                               : Stability::Stable;
    std::sort(result.modes.begin(), result.modes.end(), [](const Mode& x, const Mode& y) {
        return std::arg(x.eigenvalue) < std::arg(y.eigenvalue);
    });
    return result;
}

void LinearMapEigenAnalysis::writeReport(std::ostream& output, const Result& result) {
    std::ostringstream report;
    report << std::setprecision(16) << "One-turn eigenanalysis: ";
    switch (result.stability) {
        case Stability::Stable:
            report << "STABLE (spectral test)";
            break;
        case Stability::Unstable:
            report << "UNSTABLE (growing eigenvalue)";
            break;
        case Stability::Marginal:
            report << "MARGINAL (resonance or ill-conditioned eigenbasis)";
            break;
        case Stability::NonUnitCircle:
            report << "NON-UNIT-CIRCLE (not a conservative stable spectrum)";
            break;
    }
    report << "\nMaximum relative eigenpair residual: " << result.maximumResidual
           << "\nScaled packed-basis reciprocal condition: " << result.basisReciprocalCondition
           << "\nNear-integer warning: " << (result.nearInteger ? "yes" : "no") << '\n';
    for (unsigned i = 0; i < 4; ++i)
        report << "lambda[" << i << "] = " << result.eigenvalues[i]
               << ", modulus = " << std::abs(result.eigenvalues[i]) << '\n';
    for (unsigned i = 0; i < result.modes.size(); ++i) {
        const auto& mode = result.modes[i];
        report << "Mode " << i + 1 << ": ";
        if (mode.tune)
            report << "fractional phase " << *mode.tune << ", conjugate branch "
                   << *mode.complementaryTune << (mode.nearInteger ? " (near integer)" : "");
        else
            report << "off unit circle; no stable tune reported";
        report << '\n';
    }
    report << "Modes may be coupled; no x/y assignment. Integer tune and oriented branch are "
              "undetermined.\n";
    output << report.str();
}
