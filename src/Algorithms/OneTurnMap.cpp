// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#include "Algorithms/OneTurnMap.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <utility>
#include "Utilities/OpalException.h"

namespace {
    void require(bool valid, const char* message) {
        if (!valid) throw OpalException("OneTurnMap", message);
    }
    void checkState(const OneTurnMap::State& ray) {
        for (unsigned d = 0; d < 3; ++d)
            require(std::isfinite(ray.position(d)) && std::isfinite(ray.momentum(d))
                            && std::isfinite(ray.positionCorrection(d)),
                    "Nonfinite ray state.");
        require(std::isfinite(ray.time) && std::isfinite(ray.timeCorrection)
                        && std::isfinite(ray.pathLength) && std::isfinite(ray.pathLengthCorrection),
                "Nonfinite ray clock/path.");
    }
}  // namespace

OneTurnMap::OneTurnMap(
        const ExternalFieldRayTracker& tracker, const CoordinateSystemTrafo& section,
        const Settings& settings)
    : OneTurnMap(
              [&tracker](const State& ray, double dt, std::vector<Step>* accepted) {
                  return tracker.advance(ray, dt, accepted);
              },
              section, settings, tracker.maximumTimeStep()) {}

OneTurnMap::OneTurnMap(
        Advance advance, const CoordinateSystemTrafo& section, const Settings& settings,
        double maximumStep)
    : advance_m(std::move(advance)),
      section_m(section),
      settings_m(settings),
      maximumStep_m(maximumStep) {
    require(bool(advance_m), "Missing ray stepper.");
    require(std::isfinite(maximumStep_m) && maximumStep_m > 0, "Invalid internal step cap.");
    require(std::isfinite(settings.momentum) && settings.momentum > 0
                    && std::isfinite(settings.momentum * settings.momentum),
            "Invalid momentum.");
    require(std::isfinite(settings.dt) && settings.dt > 0 && settings.maxSteps > 0
                    && std::isfinite(settings.time),
            "Invalid integration controls.");
    require(std::isfinite(settings.maxPath) && settings.maxPath > 0
                    && std::isfinite(settings.armDistance) && settings.armDistance > 0,
            "Invalid return search bounds.");
}

double OneTurnMap::distance(const State& ray) const {
    // Apply the small compensated correction after the frame transformation.
    const auto position   = section_m.transformTo(ray.position);
    const auto correction = section_m.rotateTo(ray.positionCorrection);
    return position(2) - correction(2);
}

OneTurnMap::State OneTurnMap::locate(const State& before, double duration) const {
    double lo = 0, hi = duration;
    State end = advance_m(before, hi, nullptr);
    require(distance(before) < 0 && distance(end) >= 0, "Return is not bracketed.");
    // Time resolution follows the shared ray tracker's relative step scale.
    for (unsigned iteration = 0; iteration < 80; ++iteration) {
        const double mid = lo + (hi - lo) / 2;
        if (mid == lo || mid == hi || hi - lo <= 1e-12 * settings_m.dt) break;
        State trial = advance_m(before, mid, nullptr);
        checkState(trial);
        if (distance(trial) >= 0) {
            hi  = mid;
            end = trial;
        } else
            lo = mid;
    }
    return end;
}

OneTurnMap::Result OneTurnMap::operator()(const Coordinates& u) const {
    for (double value : u)
        require(std::isfinite(value), "Nonfinite launch coordinate.");
    const double pz2 = settings_m.momentum * settings_m.momentum - u[1] * u[1] - u[3] * u[3];
    require(std::isfinite(pz2) && pz2 > 0, "Launch has no real forward longitudinal momentum.");
    State ray;
    ray.position = section_m.transformFrom(Vector_t<double, 3>(u[0], u[2], 0));
    ray.momentum = section_m.rotateFrom(Vector_t<double, 3>(u[1], u[3], std::sqrt(pz2)));
    ray.time     = settings_m.time;
    checkState(ray);
    bool armed = false;
    std::vector<Step> accepted;
    for (unsigned step = 0; step < settings_m.maxSteps; ++step) {
        accepted.clear();
        const State end = advance_m(ray, settings_m.dt, &accepted);
        require(!accepted.empty(), "Stepper returned no accepted intervals.");
        for (const auto& interval : accepted) {
            require(std::isfinite(interval.duration) && interval.duration > 0,
                    "Stepper returned an invalid interval duration.");
            require(!interval.hitMaterial, "Ray intercepted material.");
            checkState(interval.end);
            const double before = distance(ray), after = distance(interval.end);
            if (before < -settings_m.armDistance || after < -settings_m.armDistance) armed = true;
            if (armed && before < 0 && after >= 0) {
                State crossing = locate(ray, interval.duration);
                const auto p   = section_m.rotateTo(crossing.momentum);
                require(p(2) > 0, "Return is not forward directed.");
                require(crossing.pathLength - crossing.pathLengthCorrection <= settings_m.maxPath,
                        "Return exceeds path budget.");
                const Vector_t<double, 3> r = section_m.transformTo(crossing.position)
                                              - section_m.rotateTo(crossing.positionCorrection);
                return {{r(0), p(0), r(1), p(1)}, crossing, step + 1, r(2)};
            }
            require(interval.end.pathLength - interval.end.pathLengthCorrection
                            <= settings_m.maxPath,
                    "No directed return within path budget.");
            ray = interval.end;
        }
        ray = end;
        checkState(ray);
    }
    throw OpalException("OneTurnMap", "No directed return within MAXSTEPS.");
}

OneTurnMap::Matrix OneTurnMap::jacobian(const Coordinates& u, const Coordinates& steps) const {
    for (double h : steps)
        require(std::isfinite(2 * h) && h > 0, "Invalid finite-difference step.");
    Matrix matrix{};
    for (unsigned column = 0; column < 4; ++column) {
        auto plus = u, minus = u;
        plus[column] += steps[column];
        minus[column] -= steps[column];
        require(plus[column] != u[column] && minus[column] != u[column],
                "Finite-difference step is below coordinate resolution.");
        const auto a = (*this)(plus).coordinates, b = (*this)(minus).coordinates;
        for (unsigned row = 0; row < 4; ++row) {
            matrix[row][column] = (a[row] - b[row]) / (2 * steps[column]);
            require(std::isfinite(matrix[row][column]), "Nonfinite return-map derivative.");
        }
    }
    return matrix;
}

OneTurnMap::FullTunes OneTurnMap::fullTunes(
        const Coordinates& u, const Coordinates& steps, const Matrix& matrix) const {
    FullTunes result;
    try {
        const double pi = std::acos(-1.0);
        std::array<double, 2> alpha{}, beta{}, endpoint{};
        // Dimensionless similarity scaling: one metre, p0, one metre, p0.
        const Coordinates scale{1, settings_m.momentum, 1, settings_m.momentum};
        double diagonal = 1, coupling = 0;
        for (unsigned i = 0; i < 4; ++i)
            for (unsigned j = 0; j < 4; ++j) {
                const double value = std::abs(matrix[i][j] * scale[j] / scale[i]);
                if (i / 2 == j / 2)
                    diagonal = std::max(diagonal, value);
                else
                    coupling = std::max(coupling, value);
            }
        require(coupling <= 1e-7 * diagonal, "coupled transverse map.");
        for (unsigned plane = 0; plane < 2; ++plane) {
            const unsigned i = 2 * plane;
            const double determinant =
                    matrix[i][i] * matrix[i + 1][i + 1] - matrix[i][i + 1] * matrix[i + 1][i];
            require(determinant > 0 && std::abs(std::sqrt(determinant) - 1) <= 1e-8,
                    "non-unit-circle or inaccurate plane map.");
            const double modulus = std::sqrt(determinant);
            const double cosine  = (matrix[i][i] + matrix[i + 1][i + 1]) / (2 * modulus);
            require(std::abs(cosine) < 1 - 1e-10, "unstable or resonant plane map.");
            const double sine = std::copysign(std::sqrt(1 - cosine * cosine), matrix[i][i + 1]);
            beta[plane]       = matrix[i][i + 1] / (modulus * sine);
            require(beta[plane] > 0, "degenerate plane map.");
            alpha[plane]    = (matrix[i][i] - matrix[i + 1][i + 1]) / (2 * modulus * sine);
            endpoint[plane] = std::atan2(sine, cosine);
        }
        const auto returned   = (*this)(u);
        const double duration = returned.ray.time - returned.ray.timeCorrection - settings_m.time;
        const auto windingAtStep = [&](double phaseDt) {
            FullTunes result;
            std::array<double, 2> phase{}, previous{};
            auto launch = [&](const Coordinates& v) {
                State ray;
                const double pz2 =
                        settings_m.momentum * settings_m.momentum - v[1] * v[1] - v[3] * v[3];
                require(pz2 > 0, "Invalid full-tune shadow launch.");
                ray.position = section_m.transformFrom(Vector_t<double, 3>(v[0], v[2], 0));
                ray.momentum =
                        section_m.rotateFrom(Vector_t<double, 3>(v[1], v[3], std::sqrt(pz2)));
                ray.time = settings_m.time;
                return ray;
            };
            State reference = launch(u);
            std::array<State, 8> shadows;
            for (unsigned j = 0; j < 4; ++j) {
                require(std::isfinite(steps[j]) && steps[j] > 0, "Invalid full-tune FDSTEP.");
                auto plus = u, minus = u;
                plus[j] += steps[j];
                minus[j] -= steps[j];
                require(plus[j] != u[j] && minus[j] != u[j], "Unresolved full-tune FDSTEP.");
                shadows[2 * j]     = launch(plus);
                shadows[2 * j + 1] = launch(minus);
            }
            auto sample = [&]() {
                checkState(reference);
                const auto p                = section_m.rotateTo(reference.momentum);
                const Vector_t<double, 3> r = section_m.transformTo(reference.position)
                                              - section_m.rotateTo(reference.positionCorrection);
                require(std::abs(r(1)) < 1e-7 && std::abs(p(1)) < 1e-7 * settings_m.momentum,
                        "reference is not in the section's horizontal plane.");
                const double horizontal = std::hypot(p(0), p(2));
                require(horizontal > 0, "Degenerate full-tune moving frame.");
                // Normal to the reference velocity removes first-order arrival-time
                // displacement. Vertical is the section's fixed local y axis.
                const Vector_t<double, 3> normal(p(2) / horizontal, 0, -p(0) / horizontal);
                std::array<std::array<double, 4>, 2> projected{};
                double direct = 1, crossPlane = 0;
                for (unsigned j = 0; j < 4; ++j) {
                    const auto& a = shadows[2 * j];
                    const auto& b = shadows[2 * j + 1];
                    checkState(a);
                    checkState(b);
                    const Vector_t<double, 3> delta = section_m.rotateTo(
                            (a.position - b.position)
                            - (a.positionCorrection - b.positionCorrection));
                    projected[0][j] = dot(delta, normal) / (2 * steps[j]);
                    projected[1][j] = delta(1) / (2 * steps[j]);
                    for (unsigned plane = 0; plane < 2; ++plane) {
                        const double value = std::abs(projected[plane][j] * scale[j]);
                        if (plane == j / 2)
                            direct = std::max(direct, value);
                        else
                            crossPlane = std::max(crossPlane, value);
                    }
                }
                require(crossPlane <= 1e-7 * direct, "coupled transport along reference orbit.");
                for (unsigned plane = 0; plane < 2; ++plane) {
                    std::array<double, 2> derivative{};
                    for (unsigned column = 0; column < 2; ++column) {
                        const unsigned j   = 2 * plane + column;
                        derivative[column] = projected[plane][j];
                    }
                    const std::complex<double> z =
                            derivative[0]
                            + std::complex<double>(-alpha[plane], 1) / beta[plane] * derivative[1];
                    require(std::isfinite(std::abs(z)) && std::abs(z) > 1e-8,
                            "unresolved shadow phase.");
                    const double angle = std::arg(z);
                    if (result.samples) {
                        const double increment = std::remainder(angle - previous[plane], 2 * pi);
                        result.maximumPhaseStep =
                                std::max(result.maximumPhaseStep, std::abs(increment));
                        require(std::abs(increment) < pi / 2,
                                "phase sampling too coarse; reduce DT.");
                        phase[plane] += increment;
                    }
                    previous[plane] = angle;
                }
                ++result.samples;
            };
            sample();
            double elapsed = 0;
            for (unsigned step = 0; elapsed < duration; ++step) {
                require(step < settings_m.maxSteps + 1, "Full-tune traversal exceeds MAXSTEPS.");
                // Uniform nested grids can alias the same integer oscillations.
                // Alternate incommensurate intervals on the coarse grids; keep
                // the original DT traversal unchanged when refinement falls back.
                const double fraction =
                        phaseDt > settings_m.dt && step % 2 ? (std::sqrt(5.) - 1) / 2 : 1;
                const double dt = std::min(fraction * phaseDt, duration - elapsed);
                reference       = advance_m(reference, dt, nullptr);
                for (auto& ray : shadows)
                    ray = advance_m(ray, dt, nullptr);
                elapsed = std::min(duration, elapsed + dt);
                sample();
            }
            for (unsigned plane = 0; plane < 2; ++plane) {
                require(std::abs(std::remainder(phase[plane] - endpoint[plane], 2 * pi)) < 1e-4,
                        "tracked phase disagrees with return-map phase.");
                // Use the return matrix for the fractional phase and the continuous
                // trajectory only for the missing integer/oriented branch.
                const double winding = std::round((phase[plane] - endpoint[plane]) / (2 * pi));
                result.tunes[plane]  = winding + endpoint[plane] / (2 * pi);
            }
            result.available = true;
            result.diagnostic =
                    "Uncoupled planar phase winding; validate DT and FDSTEP convergence.";
            return result;
        };
        // Only winding is obtained from these rays; fractional phases always
        // come from the accurate endpoint matrix. Bound coarsening and compare
        // successive traversals, retaining the original DT as the fallback.
        double phaseDt = std::max(settings_m.dt, std::min(16 * settings_m.dt, duration / 256));
        // Require ample headroom to pay for two confirmation traversals.
        // A large nominal DT already split by the tracker cannot provide that gain.
        if (std::min(phaseDt, maximumStep_m) < 8 * settings_m.dt) phaseDt = settings_m.dt;
        FullTunes previous;
        for (;;) {
            try {
                auto candidate = windingAtStep(phaseDt);
                if (phaseDt == settings_m.dt
                    || (previous.available && candidate.tunes == previous.tunes))
                    return candidate;
                previous = std::move(candidate);
            } catch (const OpalException&) {
                if (phaseDt == settings_m.dt) throw;
                previous.available = false;
            }
            phaseDt = std::max(settings_m.dt, phaseDt / 2);
        }
    } catch (const OpalException& error) {
        result.diagnostic = error.what();
        for (char& c : result.diagnostic)
            if (c == '\n' || c == '\r' || c == '\t') c = ' ';
    }
    return result;
}
