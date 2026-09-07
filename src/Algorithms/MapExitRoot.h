#ifndef OPALX_MAP_EXIT_ROOT_H
#define OPALX_MAP_EXIT_ROOT_H

#include <algorithm>
#include <cmath>
#include "Utilities/OpalException.h"

namespace map_exit_detail {
    /** Internal sign-bracket search for a private ray's arrival time (seconds).
     * Interpolation proposes only a time; evaluate must integrate the complete trial.
     * A poor secant step forces bisection next, halving the bracket every two trials.
     * Stop at an exact zero, the absolute time tolerance, or adjacent floating-point
     * times. The returned time is always an evaluated point.
     */
    template <class Evaluate>
    double locate(double a, double b, double fa, double fb, double tolerance, Evaluate evaluate) {
        const auto fail = []() {
            throw OpalException("MapExitRoot::locate", "Invalid or unconverged exit-plane bracket.");
        };
        if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(fa)
            || !std::isfinite(fb) || !std::isfinite(tolerance) || tolerance <= 0.0) fail();
        if (a > b) { std::swap(a, b); std::swap(fa, fb); }
        if (fa == 0.0) return a;
        if (fb == 0.0) return b;
        if (a == b || !std::isfinite(b - a) || std::signbit(fa) == std::signbit(fb)) fail();
        bool bisect = false;
        for (unsigned iteration = 0; iteration < 96; ++iteration) {
            const double width = b - a;
            if (width <= tolerance || std::nextafter(a, b) == b)
                return std::abs(fa) < std::abs(fb) ? a : b;
            double trial = a + 0.5 * width;
            if (!bisect) {
                const double scale = std::max(std::abs(fa), std::abs(fb));
                const double fraction = -(fa / scale) / (fb / scale - fa / scale);
                const double candidate = a + fraction * width;
                // Avoid repeatedly proposing an endpoint when its residual is tiny.
                const double guarded = std::clamp(candidate, a + 0.5 * tolerance,
                                                   b - 0.5 * tolerance);
                if (std::isfinite(guarded) && guarded > a && guarded < b) trial = guarded;
            }
            const double value = evaluate(trial);
            if (!std::isfinite(value)) fail();
            if (value == 0.0) return trial;
            if (std::signbit(value) == std::signbit(fa)) { a = trial; fa = value; }
            else { b = trial; fb = value; }
            if (b - a <= tolerance) return trial;
            bisect = !bisect && b - a > 0.5 * width;
        }
        fail();
        return a;
    }
}
#endif
