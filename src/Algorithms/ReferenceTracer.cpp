//
// Class ReferenceTracer
//   Traces the reference particle of the primary beam once through the lattice before tracking.
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPALX.
//
// OPALX is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPALX. If not, see <https://www.gnu.org/licenses/>.
//
#include "Algorithms/ReferenceTracer.h"

#include "AbstractObjects/OpalData.h"
#include "Algorithms/CavityAutophaser.h"
#include "Algorithms/CompensatedSum.h"
#include "Algorithms/PartData.h"
#include "Elements/OpalBeamline.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utilities/Util.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <sstream>

extern Inform* gmsg;

ReferenceTracer::ReferenceTracer(
        const PartData& reference, OpalBeamline& beamline, const StepSizeConfig& stepSizes,
        const double dt, const double ringPeriod)
    : reference_m(reference),
      beamline_m(beamline),
      stepSizes_m(stepSizes),
      dt_m(dt),
      ringPeriod_m(ringPeriod),
      rayTracker_m(beamline, reference),
      designPathFrequency_m(std::numeric_limits<std::size_t>::max()) {
    if (!(std::isfinite(dt_m) && dt_m > 0.0))
        throw OpalException("ReferenceTracer", "The reference pass needs a positive DT.");
}

ReferenceTracer::RayState ReferenceTracer::execute(RayState start, const double bunchTail) {
    ray_m       = start;
    labBounds_m = beamline_m.getElementMapping().getLabBounds();
    openDesignPath();

    // Step back over the bunch tail so that the design path also covers it.
    const double tail = std::min(start.pathLength, std::max(0.0, bunchTail));
    ElementSet fields;
    while (std::abs(start.pathLength - ray_m.pathLength) < tail && step(-dt_m, fields)) {
        ++step_m;
    }
    step_m *= -1;

    startPathLength_m               = ray_m.pathLength;
    endPathLength_m                 = stepSizes_m.getFinalSStop() + 2.0 * std::max(0.0, bunchTail);
    maxForwardSteps_m               = stepSizes_m.getNumStepsFinestResolution();
    unsigned long long forwardSteps = 0;
    bool hitMaterial                = false;
    do {
        hitMaterial = !step(dt_m, fields);
        ++step_m;
        ++forwardSteps;
    } while (!hitMaterial && !reachedEnd(fields) && forwardSteps < maxForwardSteps_m);

    designPath_m.close();
    *gmsg << level2 << "* Reference pass: path length " << startPathLength_m << " m to "
          << ray_m.pathLength << " m in " << forwardSteps << " steps of " << dt_m << " s"
          << (hitMaterial ? "; stopped where the reference hit an element." : ".") << endl;
    return ray_m;
}

bool ReferenceTracer::step(const double dt, ElementSet& fields) {
    // The Boris step samples the fields at the midpoint of its first drift.
    Vector_t<double, 3> midpoint = ray_m.position;
    auto correction              = ray_m.positionCorrection;
    const double factor          = 0.5 * Physics::c * dt / Util::getGamma(ray_m.momentum);
    for (unsigned d = 0; d < 3; ++d)
        compensated::add(factor * ray_m.momentum(d), midpoint(d), correction(d));
    fields = beamline_m.getElements(midpoint);

    if (dt > 0.0 && fields != previousFields_m) {
        checkElementLengths(fields);
        autophaseCavities(fields);
        previousFields_m = fields;
    }

    const auto result =
            rayTracker_m.step(ray_m, dt, [&](const RayState& ray, auto& electric, auto& magnetic) {
                for (const auto& element : fields) {
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
    beamline_m.getElementMapping().checkReferencePoint(
            result.midpoint.position, result.midpoint.pathLength);
    if (result.hitMaterial) return false;

    writeDesignPath(result, fields, dt);
    ray_m = result.end;
    if (dt > 0.0) {
        for (const auto& element : fields)
            visitedElements_m.insert(element->getName());
    }
    return true;
}

void ReferenceTracer::checkElementLengths(const ElementSet& fields) {
    while (!stepSizes_m.reachedEnd() && ray_m.pathLength > stepSizes_m.getSStop()) {
        ++stepSizes_m;
    }
    if (stepSizes_m.reachedEnd()) return;

    const double driftLength = Physics::c * std::abs(stepSizes_m.getdT())
                               * euclidean_norm(ray_m.momentum) / Util::getGamma(ray_m.momentum);
    for (const auto& field : fields) {
        // Markers and monitors need no field integration.
        if (field->getType() == ElementType::MARKER || field->getType() == ElementType::MONITOR) {
            continue;
        }
        double fieldBegin = 0.0;
        double fieldEnd   = 0.0;
        field->getFieldExtent(fieldBegin, fieldEnd);

        const double length = std::abs(fieldEnd - fieldBegin);
        const int numSteps  = field->getRequiredNumberOfTimeSteps();
        if (length < numSteps * driftLength) {
            throw OpalException(
                    "ReferenceTracer::checkElementLengths",
                    "The time step is too long compared to the length of the\n"
                    "element '" + field->getName() + "'\n"
                    "The field-support extent of the element is: " + std::to_string(length)
                            + "\n"
                              "The distance the particles drift in "
                            + std::to_string(numSteps)
                            + " time step(s) is: " + std::to_string(numSteps * driftLength));
        }
    }
}

void ReferenceTracer::autophaseCavities(const ElementSet& fields) {
    if (Options::autoPhase == 0) return;
    for (const auto& element : fields) {
        if ((element->getType() == ElementType::TRAVELINGWAVE
             || element->getType() == ElementType::RFCAVITY)
            && visitedElements_m.find(element->getName()) == visitedElements_m.end()) {
            const auto initialR = beamline_m.transformToLocalCS(element, ray_m.position);
            const auto initialP = beamline_m.rotateToLocalCS(element, ray_m.momentum);
            CavityAutophaser autophaser(reference_m, element);
            autophaser.getPhaseAtMaxEnergy(initialR, initialP, ray_m.time, dt_m);
        }
    }
}

bool ReferenceTracer::reachedEnd(const ElementSet& fields) const {
    if (ringPeriod_m > 0.0) return std::abs(ray_m.pathLength - startPathLength_m) >= ringPeriod_m;
    if (!fields.empty()) return false;
    if (ray_m.pathLength >= endPathLength_m) return true;
    return !labBounds_m.isInside(ray_m.position)
           && !rayHitsBox(labBounds_m, ray_m.position, ray_m.momentum);
}

bool ReferenceTracer::rayHitsBox(
        const BoundingBox& box, const Vector_t<double, 3>& position,
        const Vector_t<double, 3>& direction) {
    const auto [lower, upper] = box.getCorners();
    double enter              = 0.0;
    double leave              = std::numeric_limits<double>::infinity();
    for (unsigned d = 0; d < 3; ++d) {
        if (!(lower(d) <= upper(d))) return false;
        if (direction(d) == 0.0) {
            if (position(d) < lower(d) || position(d) > upper(d)) return false;
            continue;
        }
        const double t1 = (lower(d) - position(d)) / direction(d);
        const double t2 = (upper(d) - position(d)) / direction(d);
        enter           = std::max(enter, std::min(t1, t2));
        leave           = std::min(leave, std::max(t1, t2));
        if (enter > leave) return false;
    }
    return true;
}

void ReferenceTracer::openDesignPath() {
    auto opal = OpalData::getInstance();
    if (ippl::Comm->rank() != 0) return;

    const std::string fileName = Util::combineFilePath(
            {opal->getAuxiliaryOutputDirectory(), opal->getInputBasename() + "_DesignPath.dat"});
    if (opal->getOpenMode() == OpalData::OpenMode::WRITE || !std::filesystem::exists(fileName)) {
        designPath_m.open(fileName);
        designPath_m << "#" << std::setw(17) << "1 - s" << std::setw(18) << "2 - Rx"
                     << std::setw(18) << "3 - Ry" << std::setw(18) << "4 - Rz" << std::setw(18)
                     << "5 - Px" << std::setw(18) << "6 - Py" << std::setw(18) << "7 - Pz"
                     << std::setw(18) << "8 - Efx" << std::setw(18) << "9 - Efy" << std::setw(18)
                     << "10 - Efz" << std::setw(18) << "11 - Bfx" << std::setw(18) << "12 - Bfy"
                     << std::setw(18) << "13 - Bfz" << std::setw(18) << "14 - Ekin" << std::setw(18)
                     << "15 - t" << std::endl;
    } else {
        designPath_m.open(fileName, std::ios_base::app);
    }
    designPathFrequency_m = std::max(1.0, std::round(1e-11 / std::abs(dt_m)));
}

void ReferenceTracer::writeDesignPath(
        const ExternalFieldRayTracker::Step& step, const ElementSet& fields, const double dt) {
    if (!designPath_m.is_open()) return;
    const double pathLength = step.midpoint.pathLength;
    const bool forwardInside =
            ray_m.pathLength > 0.0 && (ringPeriod_m > 0.0 || ray_m.pathLength < endPathLength_m);
    if (!((forwardInside || dt < 0.0) && step_m % designPathFrequency_m == 0)) return;

    std::string names("\t");
    for (const auto& element : fields)
        names += element->getName() + ", ";
    const auto& r   = step.midpoint.position;
    const auto& p   = step.midpoint.momentum;
    const auto& Ef  = step.electric;
    const auto& Bf  = step.magnetic;
    const auto cell = [this](const double value) -> std::ostream& {
        designPath_m << std::setw(18) << std::setprecision(8) << value;
        return designPath_m;
    };
    cell(pathLength);
    for (unsigned d = 0; d < 3; ++d)
        cell(r(d));
    for (unsigned d = 0; d < 3; ++d)
        cell(p(d));
    for (unsigned d = 0; d < 3; ++d)
        cell(Ef(d));
    for (unsigned d = 0; d < 3; ++d)
        cell(Bf(d));
    cell(reference_m.getM() * (std::sqrt(dot(ray_m.momentum, ray_m.momentum) + 1) - 1)
         * Units::eV2MeV);
    cell(step.midpoint.time * Units::s2ns) << names << std::endl;
}
