// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#ifndef OPAL_SPECTRAL_TUNES_H
#define OPAL_SPECTRAL_TUNES_H
#include <string>
#include <vector>
#include "Algorithms/ExternalFieldRayTracker.h"
class OpalBeamline;
class PartData;
/** Serial two-ray coasting tunes, independent of map construction.
 * Old SEO's uniform radial separation/vertical motion are analyzed with
 * mean-subtracted Lomb--Scargle power: oversampling 4, high-frequency factor
 * 0.8, maximum grid peak without fitting. N samples and legacy turn count L
 * give bins k/(4L), k=1..floor(1.6N). Sample index replaces physical time;
 * fixed sampling is required. This is neither an FFT nor a closed-orbit solve.
 * Host rays reuse ExternalFieldRayTracker and its support subdivision, no RF
 * or collective fields, no bunch/GPU kernel or transfer-map eigenvalues.
 *
 * For centered samples d and unbiased sample variance v, each power is
 * \f$P=d^T A(A^T A)^{-1}A^T d/(2v)\f$, with columns
 * \f$A_i=(\cos\omega i,\sin\omega i)\f$ and
 * \f$\omega=2\pi k/[4(N-1)]\f$. This equals the legacy phase-shifted
 * Lomb formula. Trigonometric recurrence is reanchored every 256 samples;
 * summation order differs from old OPAL, without changing its frequency grid.
 */
class SpectralTunes {
public:
    struct Peak {
        double tune, power, spacing;
    };
    struct Spectrum {
        Peak peak;
        std::vector<double> tune, power;
    };
    /// Requires at least 12 finite samples, positive variance and positive turns.
    static Spectrum analyze(const std::vector<double>&, unsigned legacyTurns);
    struct Settings {
        unsigned turns = 100, stepsPerTurn = 720, sampleEvery = 50;
        double dt = 0, radialOffset = 0.005, verticalOffset = 0.005;
        std::string sector = "SM0";
        ExternalFieldRayTracker::IntegrationMethod integrator =
                ExternalFieldRayTracker::IntegrationMethod::RK4;
    };
    /** Launch triples: kinetic GeV, radius m, pr/(mc), in the named sector's
     * entrance radial plane. Positive tangential momentum. Writes tune/signal/
     * spectrum CSVs. Geometry supplies machine centre and rotated local axes.
     */
    static void run(
            OpalBeamline&, const PartData&, const std::vector<double>&, const Settings&,
            const std::string& base);
};
#endif
