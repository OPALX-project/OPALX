// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#include "Algorithms/ExternalFieldRayTracker.h"
#include "Algorithms/PartData.h"
#include "Elements/OpalBeamline.h"

ExternalFieldRayTracker::ExternalFieldRayTracker(OpalBeamline& beamline, const PartData& reference)
    : beamline_m(beamline), reference_m(reference) {}

ExternalFieldRayTracker::Step ExternalFieldRayTracker::step(
        const State& initial, const double dt, const FieldEvaluator& fields) const {
    Step result;
    result.midpoint = result.end = initial;
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
        const State& initial, const double dt) const {
    return step(initial, dt,
                [this](const State& ray, auto& electric, auto& magnetic) {
                    beamline_m.getFieldAt(ray.position, ray.momentum, ray.time, electric, magnetic);
                    return false;
                })
            .end;
}
