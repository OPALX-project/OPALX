// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#include "Algorithms/ExternalFieldRayTracker.h"
#include "Algorithms/CompensatedSum.h"
#include "Algorithms/RungeKuttaTableau.h"
#include <algorithm>
#include <cmath>
#include "Algorithms/PartData.h"
#include "Elements/OpalBeamline.h"
#include "Utilities/OpalException.h"

namespace {
    using State = ExternalFieldRayTracker::State;
    using Fields = ExternalFieldRayTracker::FieldEvaluator;

    // Each k holds h*f for (r[3], u[3], s), avoiding sums of large SI derivatives.
    template <unsigned Stages>
    State rkAdvance(const State& initial, double h, const Fields& fields,
                    const PartData& reference, const external_field_rk::Tableau<Stages>& table,
                    bool& hitMaterial) {
        std::array<std::array<double, 7>, Stages> k{};
        const auto combine = [&](const auto& weights, unsigned count, double fraction) {
            State ray = initial;
            std::array<double, 7> increment{}, correction{};
            for (unsigned j = 0; j < count; ++j)
                for (unsigned d = 0; d < 7; ++d)
                    compensated::add(weights[j] * k[j][d], increment[d], correction[d]);
            for (unsigned d = 0; d < 3; ++d) {
                compensated::add(increment[d] - correction[d], ray.position(d),
                                 ray.positionCorrection(d));
                ray.momentum(d) += increment[d + 3] - correction[d + 3];
            }
            compensated::add(increment[6] - correction[6], ray.pathLength, ray.pathLengthCorrection);
            compensated::add(fraction * h, ray.time, ray.timeCorrection);
            return ray;
        };
        const double force = h * reference.getQ() * Physics::c / reference.getM();
        for (unsigned i = 0; i < Stages; ++i) {
            const State ray = combine(table.a[i], i, table.c[i]);
            Vector_t<double, 3> electric(0.0), magnetic(0.0);
            if (fields(ray, electric, magnetic)) {
                hitMaterial = true;
                return ray;
            }
            const double gamma = std::sqrt(1.0 + dot(ray.momentum, ray.momentum));
            const Vector_t<double, 3> velocity = (Physics::c / gamma) * ray.momentum;
            const Vector_t<double, 3> du = force * (electric + cross(velocity, magnetic));
            for (unsigned d = 0; d < 3; ++d) {
                k[i][d] = h * velocity(d);
                k[i][d + 3] = du(d);
            }
            k[i][6] = h * euclidean_norm(velocity);
        }
        return combine(table.b, Stages, 1.0);
    }
}

ExternalFieldRayTracker::IntegrationMethod ExternalFieldRayTracker::parseIntegrationMethod(
        const std::string& name) {
    if (name == "BORIS" || name == "LF2") return IntegrationMethod::BORIS;
    if (name == "RK4") return IntegrationMethod::RK4;
    if (name == "DOP853") return IntegrationMethod::DOP853;
    throw OpalException(
            "ExternalFieldRayTracker::parseIntegrationMethod",
            "Unsupported external-field ray integrator '" + name + "'. Use BORIS (alias LF2), RK4 or DOP853.");
}

std::string ExternalFieldRayTracker::integrationMethodName(const IntegrationMethod method) {
    switch (method) {
        case IntegrationMethod::BORIS:
            return "BORIS";
        case IntegrationMethod::RK4:
            return "RK4";
        case IntegrationMethod::DOP853:
            return "DOP853";
    }
    throw OpalException(
            "ExternalFieldRayTracker::integrationMethodName", "Unsupported ray integrator enum.");
}

ExternalFieldRayTracker::ExternalFieldRayTracker(
        OpalBeamline& beamline, const PartData& reference, const IntegrationMethod method)
    : beamline_m(beamline), reference_m(reference), integrationMethod_m(method) {
    integrationMethodName(method);  // Validate even when no ray will be advanced.
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
    if (!std::isfinite(dt))
        throw OpalException("ExternalFieldRayTracker::step", "The ray time step must be finite.");
    switch (integrationMethod_m) {
        case IntegrationMethod::BORIS:
            return borisStep(initial, dt, fields);
        case IntegrationMethod::RK4:
        case IntegrationMethod::DOP853:
            return rungeKuttaStep(initial, dt, fields);
    }
    throw OpalException("ExternalFieldRayTracker::step", "Unsupported ray integrator enum.");
}

ExternalFieldRayTracker::Step ExternalFieldRayTracker::rungeKuttaStep(
        const State& initial, const double dt, const FieldEvaluator& fields) const {
    Step result;
    result.midpoint = result.end = initial;
    result.duration = dt;
    if (dt == 0.0) return result;
    const auto integrate = [&](double h) {
        if (integrationMethod_m == IntegrationMethod::RK4)
            return rkAdvance(initial, h, fields, reference_m, external_field_rk::rk4, result.hitMaterial);
        return rkAdvance(initial, h, fields, reference_m, external_field_rk::dop853, result.hitMaterial);
    };
    result.end = integrate(dt);
    if (result.hitMaterial) return result;
    // An independent half-duration solve supplies a same-order diagnostic midpoint.
    // It is not composed into the endpoint: nominal DT retains its usual meaning.
    result.midpoint = integrate(0.5 * dt);
    if (result.hitMaterial) return result;
    result.hitMaterial = fields(result.midpoint, result.electric, result.magnetic);
    return result;
}

ExternalFieldRayTracker::Step ExternalFieldRayTracker::borisStep(
        const State& initial, const double dt, const FieldEvaluator& fields) const {
    Step result;
    result.midpoint = result.end = initial;
    result.duration              = dt;
    if (dt == 0.0) return result;
    const auto halfDrift = [dt](State& ray) {
        const double momentum2 = dot(ray.momentum, ray.momentum);
        const double factor = (0.5 * Physics::c * dt) / std::sqrt(1.0 + momentum2);
        for (unsigned component = 0; component < 3; ++component)
            compensated::add(factor * ray.momentum(component), ray.position(component),
                             ray.positionCorrection(component));
        compensated::add(factor * std::sqrt(momentum2), ray.pathLength,
                         ray.pathLengthCorrection);
        compensated::add(0.5 * dt, ray.time, ray.timeCorrection);
    };
    auto& ray = result.midpoint;
    halfDrift(ray);
    result.hitMaterial = fields(ray, result.electric, result.magnetic);
    result.end         = ray;
    if (result.hitMaterial) return result;
    // BorisPusher::kick does not use its position argument. Keep the shared kick;
    // only the host drift/bookkeeping changes, not TRACK's particle kernel.
    integrator_m.kick(
            ray.position, result.end.momentum, result.electric, result.magnetic, dt, reference_m.getM(),
            reference_m.getQ());
    halfDrift(result.end);
    return result;
}

ExternalFieldRayTracker::State ExternalFieldRayTracker::advanceToPathLength(
        const State& initial, const double dt, const double target) const {
    const auto distance = [target](const State& ray) {
        return compensated::difference(ray.pathLength, ray.pathLengthCorrection, target, 0.0);
    };
    if (!std::isfinite(target))
        throw OpalException("ExternalFieldRayTracker::advanceToPathLength", "Non-finite path target.");
    State trial = advance(initial, dt);
    const double direction = std::copysign(1.0, dt);
    if (direction * distance(initial) > 0.0 || direction * distance(trial) < 0.0)
        throw OpalException("ExternalFieldRayTracker::advanceToPathLength", "Path target is not bracketed.");
    if (distance(initial) == 0.0) return initial;
    if (distance(trial) == 0.0) return trial;
    const double timeTolerance = std::max(
            1.e-12 * std::abs(dt), 64.0 * std::numeric_limits<double>::epsilon()
                    * std::max(1.0, euclidean_norm(initial.position)) / Physics::c);
    double lower = 0.0, upper = dt;
    for (unsigned iteration = 0; iteration < 64; ++iteration) {
        const double middle = 0.5 * (lower + upper);
        trial = advance(initial, middle);
        if (distance(trial) == 0.0) break;
        if (direction * distance(trial) >= 0.0) upper = middle;
        else lower = middle;
        if (std::abs(upper - lower) <= timeTolerance) break;
    }
    return trial;
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
    bool crossedSupport = false;
    const auto trial = step(initial, dt, [&](const State& ray, auto& electric, auto& magnetic) {
        const auto stageSet = beamline_m.getElements(ray.position);
        crossedSupport = crossedSupport || stageSet != initialSet;
        for (const auto& element : stageSet) {
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
        && (crossedSupport || initialSet != beamline_m.getElements(trial.end.position)))
        return split();
    if (trial.hitMaterial) {
        throw OpalException(
                "ExternalFieldRayTracker::advance",
                "A ray reached material while evaluating external fields.");
    }
    if (accepted) accepted->push_back(trial);
    return trial.end;
}
