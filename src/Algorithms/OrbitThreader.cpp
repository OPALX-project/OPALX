//
// Class OrbitThreader
//
// This class determines the design path by tracking the reference particle through
// the 3D lattice.
//
// All rights reserved
//
// This file is part of OPALX.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#include "Algorithms/OrbitThreader.h"

#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/TravelingWave.h"
#include "AbstractObjects/OpalData.h"
#include "Algorithms/CavityAutophaser.h"
#include "BasicActions/Option.h"
#include "BeamlineCore/MarkerRep.h"
#include "Physics/Units.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utilities/Util.h"

#include <filesystem>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

#define HITMATERIAL 0x80000000
#define EOL 0x40000000
#define EVERYTHINGFINE 0x00000000
extern Inform* gmsg;

OrbitThreader::OrbitThreader(
        const PartData& ref, const Vector_t<double, 3>& r, const Vector_t<double, 3>& p, double s,
        double maxDiffZBunch, double t, double dt, StepSizeConfig stepSizes, OpalBeamline& bl)
    : r_m(r),
      p_m(p),
      pathLength_m(s),
      time_m(t),
      dt_m(dt),
      stepSizes_m(stepSizes),
      zstop_m(stepSizes.getFinalZStop() + std::copysign(1.0, dt) * 2 * maxDiffZBunch),
      itsOpalBeamline_m(bl),
      errorFlag_m(0),
      integrator_m{},
      reference_m(ref) {
    auto opal = OpalData::getInstance();
    if (ippl::Comm->rank() == 0 && !opal->isOptimizerRun()) {
        std::string fileName = Util::combineFilePath(
                {opal->getAuxiliaryOutputDirectory(),
                 opal->getInputBasename() + "_DesignPath.dat"});
        if (opal->getOpenMode() == OpalData::OpenMode::WRITE
            || !std::filesystem::exists(fileName)) {
            logger_m.open(fileName);
            logger_m << "#" << std::setw(17) << "1 - s" << std::setw(18) << "2 - Rx"
                     << std::setw(18) << "3 - Ry" << std::setw(18) << "4 - Rz" << std::setw(18)
                     << "5 - Px" << std::setw(18) << "6 - Py" << std::setw(18) << "7 - Pz"
                     << std::setw(18) << "8 - Efx" << std::setw(18) << "9 - Efy" << std::setw(18)
                     << "10 - Efz" << std::setw(18) << "11 - Bfx" << std::setw(18) << "12 - Bfy"
                     << std::setw(18) << "13 - Bfz" << std::setw(18) << "14 - Ekin" << std::setw(18)
                     << "15 - t" << std::endl;
        } else {
            logger_m.open(fileName, std::ios_base::app);
        }

        loggingFrequency_m = std::max(1.0, std::round(1e-11 / std::abs(dt_m)));
    } else {
        loggingFrequency_m = std::numeric_limits<size_t>::max();
    }
    pathLengthRange_m = stepSizes_m.getPathLengthRange();
    pathLengthRange_m.enlargeIfOutside(pathLength_m);
    pathLengthRange_m.enlargeIfOutside(zstop_m);

    stepRange_m.enlargeIfOutside(0);
    stepRange_m.enlargeIfOutside(stepSizes_m.getNumStepsFinestResolution());
    distTrackBack_m = std::min(pathLength_m, std::max(0.0, maxDiffZBunch));
    computeBoundingBox();
}

void OrbitThreader::logDesignPathRow(
        double pathLength, const Vector_t<double, 3>& position, const Vector_t<double, 3>& momentum,
        const Vector_t<double, 3>& electricField, const Vector_t<double, 3>& magneticField,
        double time, const std::string& activeElementNames) {
    logger_m << std::setw(18) << std::setprecision(8) << pathLength << std::setw(18)
             << std::setprecision(8) << position(0) << std::setw(18) << std::setprecision(8)
             << position(1) << std::setw(18) << std::setprecision(8) << position(2) << std::setw(18)
             << std::setprecision(8) << momentum(0) << std::setw(18) << std::setprecision(8)
             << momentum(1) << std::setw(18) << std::setprecision(8) << momentum(2) << std::setw(18)
             << std::setprecision(8) << electricField(0) << std::setw(18) << std::setprecision(8)
             << electricField(1) << std::setw(18) << std::setprecision(8) << electricField(2)
             << std::setw(18) << std::setprecision(8) << magneticField(0) << std::setw(18)
             << std::setprecision(8) << magneticField(1) << std::setw(18) << std::setprecision(8)
             << magneticField(2) << std::setw(18) << std::setprecision(8)
             << reference_m.getM() * (sqrt(dot(momentum, momentum) + 1) - 1) * Units::eV2MeV
             << std::setw(18) << std::setprecision(8) << time * Units::s2ns << activeElementNames
             << std::endl;
}

void OrbitThreader::checkElementLengths(const std::set<std::shared_ptr<Component>>& fields) {
    while (!stepSizes_m.reachedEnd() && pathLength_m > stepSizes_m.getZStop()) {
        ++stepSizes_m;
    }
    if (stepSizes_m.reachedEnd()) {
        return;
    }
    double driftLength =
            Physics::c * std::abs(stepSizes_m.getdT()) * euclidean_norm(p_m) / Util::getGamma(p_m);
    for (const std::shared_ptr<Component>& field : fields) {
        double fieldBegin = 0.0;
        double fieldEnd   = 0.0;
        field->getFieldExtend(fieldBegin, fieldEnd);
        double length = fieldEnd - fieldBegin;
        int numSteps  = field->getRequiredNumberOfTimeSteps();

        if (length < numSteps * driftLength) {
            throw OpalException("OrbitThreader::checkElementLengths",
                                "The time step is too long compared to the field-support extent of the\n"
                                "element '" + field->getName() + "'\n" +
                                "The field-support extent is: " + std::to_string(length) + "\n"
                                "The distance the particles drift in " + std::to_string(numSteps) +
                                " time step(s) is: " + std::to_string(numSteps * driftLength));
        }
    }
}

void OrbitThreader::execute() {
    double initialPathLength = pathLength_m;

    auto allElements = itsOpalBeamline_m.getElementByType(ElementType::ANY);
    std::set<std::string> visitedElements;

    trackBack();
    updateBoundingBoxWithCurrentPosition();
    pathLengthRange_m.enlargeIfOutside(pathLength_m);

    Vector_t<double, 3> nextR = r_m / (Physics::c * dt_m);
    integrator_m.push(nextR, p_m, dt_m);
    nextR = nextR * Physics::c * dt_m;
    setDesignEnergy(allElements, visitedElements);

    auto elementSet = itsOpalBeamline_m.getElements(nextR);
    std::set<std::shared_ptr<Component>> intersection, currentSet;
    errorFlag_m = EVERYTHINGFINE;

    if (ippl::Comm->rank() == 0 && !OpalData::getInstance()->isOptimizerRun()) {
        Vector_t<double, 3> Ef(0.0), Bf(0.0);
        std::string names("\t");
        for (const auto& element : elementSet) {
            Vector_t<double, 3> localR = itsOpalBeamline_m.transformToFieldLocalCS(element, r_m);
            Vector_t<double, 3> localP = itsOpalBeamline_m.rotateToFieldLocalCS(element, p_m);
            Vector_t<double, 3> localE(0.0), localB(0.0);

            if (!element->applyToReferenceParticle(localR, localP, time_m, localE, localB)) {
                names += element->getName() + ", ";
                Ef += itsOpalBeamline_m.rotateFromLocalCS(element, localE);
                Bf += itsOpalBeamline_m.rotateFromLocalCS(element, localB);
            }
        }
        logDesignPathRow(pathLength_m, r_m, p_m, Ef, Bf, time_m, names);
    }

    *gmsg << "* OrbitThreader dt_m= " << dt_m << endl;

    do {
        checkElementLengths(elementSet);
        if (containsCavity(elementSet)) {
            autophaseCavities(elementSet, visitedElements);
        }

        double initialS              = pathLength_m;
        Vector_t<double, 3> initialR = r_m;
        Vector_t<double, 3> initialP = p_m;
        double maxDistance           = computeDriftLengthToBoundingBox(elementSet, r_m, p_m);

        integrate(elementSet, maxDistance);

        *gmsg << "* OrbitThreader maxDistance= " << maxDistance << endl;
        *gmsg << "* OrbitThreader #elements  = " << elementSet.size() << endl;

        registerElement(elementSet, initialS, initialR, initialP);

        if (errorFlag_m == HITMATERIAL) {
            // Shouldn't be reached because reference particle
            // isn't stopped by collimators
            pathLength_m += std::copysign(1.0, dt_m);
        }

        imap_m.add(initialS, pathLength_m, elementSet);

        IndexMap::value_t::const_iterator it        = elementSet.begin();
        const IndexMap::value_t::const_iterator end = elementSet.end();
        for (; it != end; ++it) {
            visitedElements.insert((*it)->getName());
        }

        setDesignEnergy(allElements, visitedElements);

        currentSet = elementSet;
        if (errorFlag_m == EVERYTHINGFINE) {
            nextR = r_m / (Physics::c * dt_m);
            integrator_m.push(nextR, p_m, dt_m);
            nextR = nextR * Physics::c * dt_m;

            elementSet = itsOpalBeamline_m.getElements(nextR);
        }
        intersection.clear();
        std::set_intersection(
                currentSet.begin(), currentSet.end(), elementSet.begin(), elementSet.end(),
                std::inserter(intersection, intersection.begin()));
    } while (errorFlag_m != EOL && stepRange_m.isInside(currentStep_m)
             && !(pathLengthRange_m.isOutside(pathLength_m) && intersection.empty()
                  && !(elementSet.empty() || currentSet.empty())));

    imap_m.tidyUp(zstop_m);
    *gmsg << level1 << "\n" << imap_m << endl;
    validateVisitedElements(allElements, visitedElements, initialPathLength);
    imap_m.saveSDDS(initialPathLength);
    processElementRegister();
}

void OrbitThreader::integrate(const IndexMap::value_t& activeSet, double /*maxDrift*/) {
    CoordinateSystemTrafo labToBeamline = itsOpalBeamline_m.getCSTrafoLab2Local();
    Vector_t<double, 3> nextR;
    do {
        errorFlag_m = EVERYTHINGFINE;

        IndexMap::value_t::const_iterator it        = activeSet.begin();
        const IndexMap::value_t::const_iterator end = activeSet.end();
        Vector_t<double, 3> oldR                    = r_m;
        Vector_t<double, 3> oldP                    = p_m;
        const double oldTime                        = time_m;
        const double oldPathLength                  = pathLength_m;

        r_m /= Physics::c * dt_m;
        integrator_m.push(r_m, p_m, dt_m);
        r_m = r_m * Physics::c * dt_m;

        Vector_t<double, 3> Ef(0.0), Bf(0.0);
        std::string names("\t");
        for (; it != end; ++it) {
            Vector_t<double, 3> localR = itsOpalBeamline_m.transformToFieldLocalCS(*it, r_m);
            Vector_t<double, 3> localP = itsOpalBeamline_m.rotateToFieldLocalCS(*it, p_m);
            Vector_t<double, 3> localE(0.0), localB(0.0);

            if ((*it)->applyToReferenceParticle(
                        localR, localP, time_m + 0.5 * dt_m, localE, localB)) {
                errorFlag_m = HITMATERIAL;
                return;
            }
            names += (*it)->getName() + ", ";

            Ef += itsOpalBeamline_m.rotateFromLocalCS(*it, localE);
            Bf += itsOpalBeamline_m.rotateFromLocalCS(*it, localB);
        }

        const bool shouldLog = currentStep_m % loggingFrequency_m == 0 && ippl::Comm->rank() == 0
                               && !OpalData::getInstance()->isOptimizerRun();

        r_m /= Physics::c * dt_m;
        integrator_m.kick(r_m, p_m, Ef, Bf, dt_m, reference_m.getM(), reference_m.getQ());
        integrator_m.push(r_m, p_m, dt_m);
        r_m = r_m * Physics::c * dt_m;

        const Vector<double, 3> d = r_m - oldR;

        pathLength_m += std::copysign(euclidean_norm(d), dt_m);

        if (shouldLog) {
            if (dt_m > 0.0 && oldPathLength < zstop_m && pathLength_m > zstop_m) {
                const double alpha = (zstop_m - oldPathLength) / (pathLength_m - oldPathLength);
                logDesignPathRow(
                        zstop_m, oldR + alpha * (r_m - oldR), oldP + alpha * (p_m - oldP), Ef, Bf,
                        oldTime + alpha * dt_m, names);
            } else if ((pathLength_m > 0.0 && pathLength_m <= zstop_m) || dt_m < 0.0) {
                logDesignPathRow(pathLength_m, r_m, p_m, Ef, Bf, time_m + dt_m, names);
            }
        }

        ++currentStep_m;
        time_m += dt_m;

        nextR = r_m / (Physics::c * dt_m);
        integrator_m.push(nextR, p_m, dt_m);
        nextR = nextR * Physics::c * dt_m;

        if (activeSet.empty()
            && (pathLengthRange_m.isOutside(pathLength_m)
                || stepRange_m.isOutside(currentStep_m))) {
            errorFlag_m = EOL;
            globalBoundingBox_m.enlargeToContainPosition(r_m);
            return;
        }

    } while (activeSet == itsOpalBeamline_m.getElements(nextR));
}

bool OrbitThreader::containsCavity(const IndexMap::value_t& activeSet) {
    IndexMap::value_t::const_iterator it        = activeSet.begin();
    const IndexMap::value_t::const_iterator end = activeSet.end();

    for (; it != end; ++it) {
        if ((*it)->getType() == ElementType::TRAVELINGWAVE
            || (*it)->getType() == ElementType::RFCAVITY) {
            return true;
        }
    }

    return false;
}

void OrbitThreader::autophaseCavities(
        const IndexMap::value_t& activeSet, const std::set<std::string>& visitedElements) {
    if (Options::autoPhase == 0) return;

    IndexMap::value_t::const_iterator it        = activeSet.begin();
    const IndexMap::value_t::const_iterator end = activeSet.end();

    for (; it != end; ++it) {
        if (((*it)->getType() == ElementType::TRAVELINGWAVE
             || (*it)->getType() == ElementType::RFCAVITY)
            && visitedElements.find((*it)->getName()) == visitedElements.end()) {
            Vector_t<double, 3> initialR = itsOpalBeamline_m.transformToFieldLocalCS(*it, r_m);
            Vector_t<double, 3> initialP = itsOpalBeamline_m.rotateToFieldLocalCS(*it, p_m);

            CavityAutophaser ap(reference_m, *it);
            ap.getPhaseAtMaxEnergy(initialR, initialP, time_m, dt_m);
        }
    }
}

double OrbitThreader::getMaxDesignEnergy(const IndexMap::value_t& elementSet) const {
    IndexMap::value_t::const_iterator it        = elementSet.begin();
    const IndexMap::value_t::const_iterator end = elementSet.end();

    double designEnergy = 0.0;
    for (; it != end; ++it) {
        if ((*it)->getType() == ElementType::TRAVELINGWAVE
            || (*it)->getType() == ElementType::RFCAVITY) {
            const RFCavity* element = static_cast<const RFCavity*>((*it).get());
            designEnergy            = std::max(designEnergy, element->getDesignEnergy());
        }
    }

    return designEnergy;
}

void OrbitThreader::trackBack() {
    dt_m *= -1;
    ValueRange<double> tmpRange;
    std::swap(tmpRange, pathLengthRange_m);
    double initialPathLength = pathLength_m;

    Vector_t<double, 3> nextR = r_m / (Physics::c * dt_m);
    integrator_m.push(nextR, p_m, dt_m);
    nextR = nextR * Physics::c * dt_m;

    while (std::abs(initialPathLength - pathLength_m) < distTrackBack_m) {
        auto elementSet = itsOpalBeamline_m.getElements(nextR);

        double maxDrift = computeDriftLengthToBoundingBox(elementSet, r_m, -p_m);
        maxDrift        = std::min(maxDrift, distTrackBack_m);
        integrate(elementSet, maxDrift);

        nextR = r_m / (Physics::c * dt_m);
        integrator_m.push(nextR, p_m, dt_m);
        nextR = nextR * Physics::c * dt_m;
    }
    std::swap(tmpRange, pathLengthRange_m);
    currentStep_m *= -1;
    dt_m *= -1;

    stepRange_m.enlargeIfOutside(currentStep_m);
}

void OrbitThreader::registerElement(
        const IndexMap::value_t& elementSet, double start, const Vector_t<double, 3>& R,
        const Vector_t<double, 3>& P) {
    IndexMap::value_t::const_iterator it        = elementSet.begin();
    const IndexMap::value_t::const_iterator end = elementSet.end();

    for (; it != end; ++it) {
        bool found = false;
        auto prior = elementRegistry_m.equal_range(*it);

        for (auto pit = prior.first; pit != prior.second; ++pit) {
            if (std::abs((*pit).second.endField_m - start) < 1e-10) {
                found                    = true;
                (*pit).second.endField_m = pathLength_m;
                break;
            }
        }

        if (found) continue;

        Vector_t<double, 3> initialR = itsOpalBeamline_m.transformToFieldLocalCS(*it, R);
        Vector_t<double, 3> initialP = itsOpalBeamline_m.rotateToFieldLocalCS(*it, P);
        double elementEdge           = start - initialR(2) * euclidean_norm(initialP) / initialP(2);

        elementPosition ep = {start, pathLength_m, elementEdge};
        elementRegistry_m.insert(std::make_pair(*it, ep));
    }
}

void OrbitThreader::processElementRegister() {
    using registry_key_t = std::shared_ptr<Component>;
    using registry_map_t = std::map<
            registry_key_t, std::set<elementPosition, elementPositionComp>,
            std::owner_less<registry_key_t>>;
    registry_map_t tmpRegistry;

    for (auto it = elementRegistry_m.begin(); it != elementRegistry_m.end(); ++it) {
        const registry_key_t& element = (*it).first;
        elementPosition& ep           = (*it).second;

        auto prior = tmpRegistry.find(element);
        if (prior == tmpRegistry.end()) {
            std::set<elementPosition, elementPositionComp> tmpSet;
            tmpSet.insert(ep);
            tmpRegistry.insert(std::make_pair(element, tmpSet));
            continue;
        }

        std::set<elementPosition, elementPositionComp>& set = (*prior).second;
        set.insert(ep);
    }

    actionRangeRegistrationModel_m.clear();
    std::vector<ReferencePathSegment> registeredSegments;
    for (auto& [element, set] : tmpRegistry) {
        std::queue<std::pair<double, double>> range;

        for (auto sit = set.begin(); sit != set.end(); ++sit) {
            range.push(std::make_pair((*sit).elementEdge_m, (*sit).endField_m));
            registeredSegments.push_back(ReferencePathSegment(
                    (*sit).startField_m, (*sit).endField_m,
                    ReferencePathSegment::element_set_t{element}, (*sit).elementEdge_m));
        }
        element->setActionRange(range);
    }

    std::sort(
            registeredSegments.begin(), registeredSegments.end(),
            [](const ReferencePathSegment& lhs, const ReferencePathSegment& rhs) {
                if (lhs.getBegin() < rhs.getBegin()) {
                    return true;
                }
                if (lhs.getBegin() > rhs.getBegin()) {
                    return false;
                }
                return lhs.getEnd() < rhs.getEnd();
            });

    for (const auto& segment : registeredSegments) {
        actionRangeRegistrationModel_m.addSegment(segment);
    }
}

void OrbitThreader::setDesignEnergy(
        FieldList& allElements, const std::set<std::string>& visitedElements) {
    double kineticEnergyeV = reference_m.getM() * (sqrt(dot(p_m, p_m) + 1.0) - 1.0);

    FieldList::iterator it        = allElements.begin();
    const FieldList::iterator end = allElements.end();
    for (; it != end; ++it) {
        std::shared_ptr<Component> element = (*it).getElement();
        if (visitedElements.find(element->getName()) == visitedElements.end()
            && !(element->getType() == ElementType::RFCAVITY
                 || element->getType() == ElementType::TRAVELINGWAVE)) {
            element->setDesignEnergy(kineticEnergyeV);
        }
    }
}

void OrbitThreader::validateVisitedElements(
        const FieldList& allElements, const std::set<std::string>& visitedElements,
        double initialPathLength) const {
    const double lowerTrackedPathLength = std::min(initialPathLength, zstop_m) - 1.0e-9;
    const double upperTrackedPathLength = std::max(initialPathLength, zstop_m) + 1.0e-9;

    auto requiresReferencePathVisit = [](const std::shared_ptr<const Component>& element) {
        switch (element->getType()) {
            case ElementType::BEAMLINE:
            case ElementType::MARKER:
            case ElementType::PROBE:
            case ElementType::SOURCE:
            case ElementType::VACUUM:
                return false;
            default:
                break;
        }

        double fieldBegin = 0.0;
        double fieldEnd   = 0.0;
        element->getFieldExtend(fieldBegin, fieldEnd);
        return std::abs(fieldEnd - fieldBegin) > 1.0e-12
               || std::abs(element->getElementLength()) > 1.0e-12;
    };

    std::ostringstream diagnostic;
    bool foundUnvisitedElement = false;

    for (auto it = allElements.begin(); it != allElements.end(); ++it) {
        const std::shared_ptr<const Component> element = (*it).getElement();
        if (!requiresReferencePathVisit(element)
            || visitedElements.find(element->getName()) != visitedElements.end()) {
            continue;
        }

        const PlacedElement placed =
                itsOpalBeamline_m.getPlacedElement(std::const_pointer_cast<Component>(element));
        const Vector_t<double, 3> entry = placed.getNominalEntryTransform().getOrigin();
        const Vector_t<double, 3> body  = placed.getNominalBodyTransform().getOrigin();
        const Vector_t<double, 3> exit  = placed.getNominalExitTransform().getOrigin();

        const std::array<double, 3> zCoordinates = {entry(2), body(2), exit(2)};
        const bool intersectsTrackedWindow =
                std::any_of(zCoordinates.begin(), zCoordinates.end(), [&](double z) {
                    return z >= lowerTrackedPathLength && z <= upperTrackedPathLength;
                });
        if (!intersectsTrackedWindow) {
            continue;
        }

        foundUnvisitedElement = true;
        diagnostic << "Element '" << element->getName() << "' (" << element->getTypeString()
                   << ") is part of the beam line but is never reached by the traced reference "
                      "path.\n"
                   << "  nominal entry = (" << entry(0) << ", " << entry(1) << ", " << entry(2)
                   << ")\n"
                   << "  nominal body  = (" << body(0) << ", " << body(1) << ", " << body(2)
                   << ")\n"
                   << "  nominal exit  = (" << exit(0) << ", " << exit(1) << ", " << exit(2)
                   << ")\n";
    }

    if (foundUnvisitedElement) {
        throw OpalException(
                "OrbitThreader::validateVisitedElements",
                "One or more placed elements are not intersected by the traced reference path.\n"
                "Adjust the input placement so the design trajectory reaches the element body or "
                "field support.\n"
                        + diagnostic.str());
    }
}

void OrbitThreader::computeBoundingBox() {
    FieldList allElements         = itsOpalBeamline_m.getElementByType(ElementType::ANY);
    FieldList::iterator it        = allElements.begin();
    const FieldList::iterator end = allElements.end();
    for (; it != end; ++it) {
        if (it->getElement()->getType() == ElementType::MARKER) {
            continue;
        }
        BoundingBox other = it->getBoundingBoxInLabCoords();
        globalBoundingBox_m.enlargeToContainBoundingBox(other);
    }
    updateBoundingBoxWithCurrentPosition();
}

void OrbitThreader::updateBoundingBoxWithCurrentPosition() {
    Vector_t<double, 3> dR                       = Physics::c * dt_m * p_m / Util::getGamma(p_m);
    std::array<Vector_t<double, 3>, 2> positions = {r_m - 10 * dR, r_m + 10 * dR};

    for (const Vector_t<double, 3>& pos : positions) {
        globalBoundingBox_m.enlargeToContainPosition(pos);
    }
}

double OrbitThreader::computeDriftLengthToBoundingBox(
        const std::set<std::shared_ptr<Component>>& elements, const Vector_t<double, 3>& position,
        const Vector_t<double, 3>& direction) const {
    if (elements.empty()
        || (elements.size() == 1 && (*elements.begin())->getType() == ElementType::DRIFT)) {
        std::optional<Vector_t<double, 3>> intersectionPoint =
                globalBoundingBox_m.getIntersectionPoint(position, direction);
        if (intersectionPoint) {
            const Vector_t<double, 3> r = intersectionPoint.value() - position;
            return euclidean_norm(r);
        }
        return 10;
    }

    return std::numeric_limits<double>::max();
}
