// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#ifndef OPAL_LINEAR_MAP_EIGEN_ANALYSIS_H
#define OPAL_LINEAR_MAP_EIGEN_ANALYSIS_H
#include <complex>
#include <iosfwd>
#include <optional>
#include <vector>
#include "Algorithms/OneTurnMap.h"

/** LAPACKE analysis of a real nonsymmetric four-dimensional one-turn map.
 * Analyze the complete matrix, including transverse coupling. Coordinates must
 * use the same entrance/exit frame at a closed orbit. Kinetic coordinates need
 * not obey the canonical-J symplectic test. This class does not establish orbit
 * closure or integration accuracy, and does not modify the input or pusher.
 *
 * A positive-imaginary eigenvalue defines one representative fractional phase
 * q=arg(lambda)/(2*pi) in (0,1/2). Its conjugate gives 1-q. Neither integer tune
 * nor oriented branch is determined here; modes are not labelled x/y or r/y.
 * The packed eigenvector basis condition and residuals are measured after the
 * similarity scaling S^-1 M S. Host-only reference LAPACK, no MPI/GPU work.
 */
class LinearMapEigenAnalysis {
public:
    using Matrix  = OneTurnMap::Matrix;
    using Complex = std::complex<double>;
    enum class Stability { Stable, Unstable, Marginal, NonUnitCircle };
    struct Settings {
        OneTurnMap::Coordinates scales{1, 1, 1, 1};  ///< Metres and p/(mc) scales.
        double modulusTolerance     = 1e-8;   ///< Numerical unit-circle band, not physics accuracy.
        double residualTolerance    = 1e-10;  ///< Relative scaled eigenpair residual ceiling.
        double nearIntegerTolerance = 1e-6;   ///< Fractional turns.
        double minimumBasisRcond    = 1e-10;  ///< Below this, no stable classification.
    };
    struct Mode {
        Complex eigenvalue;                       ///< Positive-imaginary representative.
        std::optional<double> tune;               ///< Only populated for a unit-circle pair.
        std::optional<double> complementaryTune;  ///< 1-tune, not a second mode.
        bool nearInteger = false;
    };
    struct Result {
        Stability stability = Stability::Marginal;
        std::array<Complex, 4> eigenvalues{};
        /// Right eigenvectors of S^-1 M S; vectors stored by column index first.
        std::array<std::array<Complex, 4>, 4> eigenvectors{};
        std::vector<Mode> modes;  ///< Sorted by positive-imaginary phase; no plane identity.
        double maximumResidual          = 0;
        double basisReciprocalCondition = 0;  ///< SVD of real-packed eigenvector basis.
        bool nearInteger                = false;
    };
    /// Throws OpalException for invalid input, LAPACK failure or inaccurate eigenpairs.
    static Result analyze(const Matrix&, const Settings&);
    /// Human-readable diagnostic report; phases are not integer or plane-labelled tunes.
    static void writeReport(std::ostream&, const Result&);
};
#endif
