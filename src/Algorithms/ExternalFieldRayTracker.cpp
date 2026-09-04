// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#include "Algorithms/ExternalFieldRayTracker.h"
#include <algorithm>
#include <cmath>
#include "Algorithms/PartData.h"
#include "Elements/OpalBeamline.h"
#include "Utilities/OpalException.h"

ExternalFieldRayTracker::ExternalFieldRayTracker(OpalBeamline& beamline, const PartData& reference)
    : beamline_m(beamline), reference_m(reference) {
    for (const auto& element : beamline_m.getElements()) {
        if (element->getType() == ElementType::MARKER || element->getType() == ElementType::MONITOR)
            continue;
        double begin = 0.0, end = 0.0;
        element->getFieldExtent(begin, end);
        const double length = std::abs(end - begin);
        if (length > 0.0) maximumStep_m = std::min(maximumStep_m, length / (4.0 * Physics::c));
    }
}

ExternalFieldRayTracker::Step ExternalFieldRayTracker::step(
        const State& initial, const double dt, const FieldEvaluator& fields) const {
    Step result;
    result.midpoint = result.end = initial;
    result.duration              = dt;
    if (dt == 0.0) return result;
    auto& ray                  = result.midpoint;
    Vector_t<double, 3> scaled = ray.position / (Physics::c * dt);
    integrator_m.push(scaled, ray.momentum, dt);
    ray.position = scaled * Physics::c * dt;
    ray.time += 0.5 * dt;
    result.hitMaterial = fields(ray, result.electric, result.magnetic);
    result.end         = ray;
    if (result.hitMaterial) return result;
    scaled = ray.position / (Physics::c * dt);
    integrator_m.kick(
            scaled, result.end.momentum, result.electric, result.magnetic, dt, reference_m.getM(),
            reference_m.getQ());
    integrator_m.push(scaled, result.end.momentum, dt);
    result.end.position = scaled * Physics::c * dt;
    result.end.time     = initial.time + dt;
    return result;
}

ExternalFieldRayTracker::State ExternalFieldRayTracker::advance(
        const State& initial, const double dt, std::vector<Step>* accepted) const {
    if (!std::isfinite(dt)) {
        throw OpalException(
                "ExternalFieldRayTracker::advance", "The ray time step must be finite.");
    }
    if (dt == 0.0) return initial;
    return advanceRecursive(initial, dt, 1.0e-12 * std::abs(dt), accepted, 0);
}

ExternalFieldRayTracker::State ExternalFieldRayTracker::advanceRecursive(
        const State& initial, const double dt, const double tolerance, std::vector<Step>* accepted,
        const unsigned depth) const {
    const double timeFloor = 64.0 * std::numeric_limits<double>::epsilon()
                             * std::max(1.0, euclidean_norm(initial.position)) / Physics::c;
    const bool resolved = std::abs(dt) <= std::max(tolerance, timeFloor);
    const auto split    = [&]() {
        if (depth >= 64) {
            throw OpalException(
                    "ExternalFieldRayTracker::advance",
                    "Field-boundary refinement did not converge.");
        }
        const auto middle = advanceRecursive(initial, 0.5 * dt, tolerance, accepted, depth + 1);
        return advanceRecursive(middle, 0.5 * dt, tolerance, accepted, depth + 1);
    };
    if (!resolved && std::abs(dt) > maximumStep_m) return split();

    const auto initialSet = beamline_m.getElements(initial.position);
    auto middleSet        = initialSet;
    const auto trial = step(initial, dt, [&](const State& ray, auto& electric, auto& magnetic) {
        middleSet = beamline_m.getElements(ray.position);
        for (const auto& element : middleSet) {
            if (element->getType() == ElementType::MARKER
                || element->getType() == ElementType::MONITOR)
                continue;
            const auto localR = beamline_m.transformToLocalCS(element, ray.position);
            const auto localP = beamline_m.rotateToLocalCS(element, ray.momentum);
            Vector_t<double, 3> localE(0.0), localB(0.0);
            if (element->applyToReferenceParticle(localR, localP, ray.time, localE, localB))
                return true;
            electric += beamline_m.rotateFromLocalCS(element, localE);
            magnetic += beamline_m.rotateFromLocalCS(element, localB);
        }
        return false;
    });
    if (!resolved
        && (middleSet != initialSet || middleSet != beamline_m.getElements(trial.end.position)))
        return split();
    if (trial.hitMaterial) {
        throw OpalException(
                "ExternalFieldRayTracker::advance",
                "A ray reached material while evaluating external fields.");
    }
    if (accepted) accepted->push_back(trial);
    return trial.end;
}
