// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#include "Algorithms/OneTurnMap.h"
#include <cmath>
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
              section, settings) {}

OneTurnMap::OneTurnMap(
        Advance advance, const CoordinateSystemTrafo& section, const Settings& settings)
    : advance_m(std::move(advance)), section_m(section), settings_m(settings) {
    require(bool(advance_m), "Missing ray stepper.");
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
