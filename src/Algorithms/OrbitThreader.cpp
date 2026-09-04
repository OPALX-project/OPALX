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
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <vector>

#define HITMATERIAL 0x80000000
#define EOL 0x40000000
#define EVERYTHINGFINE 0x00000000
extern Inform* gmsg;

namespace {
    constexpr double mapDiagnosticTolerance = 1.0e-6;
}

OrbitThreader::OrbitThreader(
        const PartData& ref, const Vector_t<double, 3>& r, const Vector_t<double, 3>& p, double s,
        double maxDiffZBunch, double t, double dt, StepSizeConfig stepSizes, OpalBeamline& bl,
        bool isDesignBeam, double period)
    : r_m(r),
      p_m(p),
      pathLength_m(s),
      time_m(t),
      dt_m(dt),
      stepSizes_m(stepSizes),
      sStop_m(period > 0.0
                      ? s + std::copysign(period, dt)
                      : stepSizes.getFinalSStop() + std::copysign(1.0, dt) * 2 * maxDiffZBunch),
      period_m(period),
      itsOpalBeamline_m(bl),
      isDesignBeam_m(isDesignBeam),
      errorFlag_m(0),
      integrator_m{},
      reference_m(ref),
      rayTracker_m(bl, ref) {
    if (period_m > 0.0) {
        imap_m.setPeriod(pathLength_m, period_m);
    }
    auto opal = OpalData::getInstance();
    // Only the design beam writes the _DesignPath.dat trajectory log; secondary species
    // must not open (and truncate) it.
    if (isDesignBeam_m && ippl::Comm->rank() == 0 && !opal->isOptimizerRun()) {
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
    pathLengthRange_m.enlargeIfOutside(sStop_m);

    stepRange_m.enlargeIfOutside(0);
    stepRange_m.enlargeIfOutside(stepSizes_m.getNumStepsFinestResolution());
    distTrackBack_m = std::min(pathLength_m, std::max(0.0, maxDiffZBunch));
    computeBoundingBox();
}

void OrbitThreader::checkElementLengths(const std::set<std::shared_ptr<ElementBase>>& fields) {
    while (!stepSizes_m.reachedEnd() && pathLength_m > stepSizes_m.getSStop()) {
        ++stepSizes_m;
    }
    if (stepSizes_m.reachedEnd()) {
        return;
    }
    double driftLength =
            Physics::c * std::abs(stepSizes_m.getdT()) * euclidean_norm(p_m) / Util::getGamma(p_m);
    for (const std::shared_ptr<ElementBase>& field : fields) {
        // Markers and monitors carry no field integration requirement. Their finite display
        // geometry must not constrain the reference-particle time step at the RING seam.
        if (field->getType() == ElementType::MARKER || field->getType() == ElementType::MONITOR) {
            continue;
        }
        double fieldBegin = 0.0;
        double fieldEnd   = 0.0;
        field->getFieldExtent(fieldBegin, fieldEnd);

        const double length = std::abs(fieldEnd - fieldBegin);
        const int numSteps  = field->getRequiredNumberOfTimeSteps();

        if (length < numSteps * driftLength) {
            throw OpalException("OrbitThreader::checkElementLengths",
                                "The time step is too long compared to the length of the\n"
                                "element '" + field->getName() + "'\n" +
                                "The field-support extent of the element is: "
                                        + std::to_string(length) + "\n"
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

    const bool calculateMaps = isDesignBeam_m && Options::enableLinearTransferMaps;
    if (calculateMaps) {
        for (const auto& element : itsOpalBeamline_m.getElements()) {
            element->clearLinearTransferMaps();
            element->setOverlapping(false);
        }
        referenceSamples_m.clear();
        combinedLinearTransferMap_m.reset();
        combinedDeterminantResidual_m.reset();
        combinedSymplecticResidual_m.reset();
        transferMapStartPathLength_m = initialPathLength;
        collectReferenceSamples_m    = true;
        recordReferenceSample();
    }

    Vector_t<double, 3> nextR = r_m / (Physics::c * dt_m);
    integrator_m.push(nextR, p_m, dt_m);
    nextR = nextR * Physics::c * dt_m;
    if (isDesignBeam_m) {
        setDesignEnergy(allElements, visitedElements);
    }

    auto elementSet = itsOpalBeamline_m.getElements(nextR);
    std::set<std::shared_ptr<ElementBase>> intersection, currentSet;
    errorFlag_m = EVERYTHINGFINE;

    *gmsg << "* OrbitThreader dt_m= " << dt_m << endl;

    do {
        checkElementLengths(elementSet);
        if (isDesignBeam_m && containsCavity(elementSet)) {
            autophaseCavities(elementSet, visitedElements);
        }

        double initialS    = pathLength_m;
        double maxDistance = computeDriftLengthToBoundingBox(elementSet, r_m, p_m);

        integrate(elementSet, maxDistance);

        if (errorFlag_m == HITMATERIAL) {
            // Shouldn't be reached because reference particle
            // isn't stopped by collimators
            pathLength_m += std::copysign(1.0, dt_m);
        }

        const double finalS = reachedThreadingEnd() ? sStop_m : pathLength_m;
        imap_m.add(initialS, finalS, elementSet);

        // Store overlap participation on the runtime occurrences during the reference pass.
        // Ignore the backwards pre-roll before the requested map origin.
        if (calculateMaps && finalS > transferMapStartPathLength_m && elementSet.size() > 1) {
            for (const auto& element : elementSet) {
                element->setOverlapping(true);
            }
        }

        IndexMap::value_t::const_iterator it        = elementSet.begin();
        const IndexMap::value_t::const_iterator end = elementSet.end();
        for (; it != end; ++it) {
            visitedElements.insert((*it)->getName());
        }

        if (isDesignBeam_m) {
            setDesignEnergy(allElements, visitedElements);
        }

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
    } while (errorFlag_m != EOL
             && (collectReferenceSamples_m || period_m > 0.0 || stepRange_m.isInside(currentStep_m))
             && !(pathLengthRange_m.isOutside(pathLength_m) && intersection.empty()
                  && !(elementSet.empty() || currentSet.empty())));

    imap_m.tidyUp(sStop_m);
    collectReferenceSamples_m = false;
    if (calculateMaps) {
        LinearTransferMapBuilder builder(itsOpalBeamline_m, reference_m, dt_m);
        auto result = builder.build(referenceSamples_m, transferMapStartPathLength_m);
        std::map<std::shared_ptr<ElementBase>, std::size_t,
                 std::owner_less<std::shared_ptr<ElementBase>>> passes;
        for (auto& segment : result.segments) {
            for (const auto& element : segment.owners) {
                auto attached = segment.map;
                attached.pass = passes[element]++;
                element->addLinearTransferMap(std::move(attached));
            }
        }
        combinedLinearTransferMap_m   = result.combined;
        combinedDeterminantResidual_m = result.determinantResidual;
        combinedSymplecticResidual_m  = result.symplecticResidual;
        printCombinedLinearTransferMap();
    }
    *gmsg << level1 << "\n" << imap_m << endl;
    if (isDesignBeam_m) {
        // Geometry SDDS dump is a design-beam output; secondary species only build their map.
        imap_m.saveSDDS(initialPathLength);
    }
    logger_m.close();
}

void OrbitThreader::integrate(const IndexMap::value_t& activeSet, double /*maxDrift*/) {
    CoordinateSystemTrafo labToBeamline = itsOpalBeamline_m.getCSTrafoLab2Local();
    Vector_t<double, 3> nextR;
    do {
        errorFlag_m = EVERYTHINGFINE;

        const Vector_t<double, 3> oldR = r_m;
        const Vector_t<double, 3> oldP = p_m;
        const double oldTime = time_m;
        const double oldPathLength = pathLength_m;
        std::string names("\t");
        const auto step = rayTracker_m.step(
                {r_m, p_m, time_m}, dt_m,
                [&](const RayState& ray, auto& electric, auto& magnetic) {
                    for (const auto& element : activeSet) {
                        const auto localR = itsOpalBeamline_m.transformToLocalCS(element, ray.position);
                        const auto localP = itsOpalBeamline_m.rotateToLocalCS(element, ray.momentum);
                        Vector_t<double, 3> localE(0.0), localB(0.0);
                        if (element->applyToReferenceParticle(
                                    localR, localP, ray.time, localE, localB)) return true;
                        names += element->getName() + ", ";
                        electric += itsOpalBeamline_m.rotateFromLocalCS(element, localE);
                        magnetic += itsOpalBeamline_m.rotateFromLocalCS(element, localB);
                    }
                    return false;
                });
        r_m = step.midpoint.position;
        if (step.hitMaterial) {
            errorFlag_m = HITMATERIAL;
            return;
        }
        const auto& Ef = step.electric;
        const auto& Bf = step.magnetic;

        if (((pathLength_m > 0.0 && pathLength_m < sStop_m) || dt_m < 0.0)
            && currentStep_m % loggingFrequency_m == 0 && ippl::Comm->rank() == 0
            && !OpalData::getInstance()->isOptimizerRun()) {
            const Vector<double, 3> d = r_m - oldR;

            logger_m << std::setw(18) << std::setprecision(8)
                     << pathLength_m + std::copysign(euclidean_norm(d), dt_m) << std::setw(18)
                     << std::setprecision(8) << r_m(0) << std::setw(18) << std::setprecision(8)
                     << r_m(1) << std::setw(18) << std::setprecision(8) << r_m(2) << std::setw(18)
                     << std::setprecision(8) << p_m(0) << std::setw(18) << std::setprecision(8)
                     << p_m(1) << std::setw(18) << std::setprecision(8) << p_m(2) << std::setw(18)
                     << std::setprecision(8) << Ef(0) << std::setw(18) << std::setprecision(8)
                     << Ef(1) << std::setw(18) << std::setprecision(8) << Ef(2) << std::setw(18)
                     << std::setprecision(8) << Bf(0) << std::setw(18) << std::setprecision(8)
                     << Bf(1) << std::setw(18) << std::setprecision(8) << Bf(2) << std::setw(18)
                     << std::setprecision(8)
                     << reference_m.getM() * (sqrt(dot(p_m, p_m) + 1) - 1) * Units::eV2MeV
                     << std::setw(18) << std::setprecision(8) << (time_m + 0.5 * dt_m) * Units::s2ns
                     << names << std::endl;
        }

        r_m = step.end.position;
        p_m = step.end.momentum;

        const Vector<double, 3> d = r_m - oldR;

        pathLength_m += std::copysign(euclidean_norm(d), dt_m);

        ++currentStep_m;
        time_m += dt_m;
        if (collectReferenceSamples_m) {
            if (reachedThreadingEnd() && pathLength_m != oldPathLength) {
                const double fraction  = (sStop_m - oldPathLength) / (pathLength_m - oldPathLength);
                const RayState clipped = rayTracker_m.advance({oldR, oldP, oldTime}, fraction * dt_m);
                LinearTransferMapReference state =
                        LinearTransferMapBuilder::transportFrame(referenceSamples_m.back().state, clipped.momentum);
                state.position   = clipped.position;
                state.momentum   = clipped.momentum;
                state.time       = clipped.time;
                state.pathLength = sStop_m;
                referenceSamples_m.push_back({state});
            } else {
                recordReferenceSample();
            }
        }

        if (reachedThreadingEnd()) {
            errorFlag_m = EOL;
            globalBoundingBox_m.enlargeToContainPosition(r_m);
            return;
        }

        nextR = r_m / (Physics::c * dt_m);
        integrator_m.push(nextR, p_m, dt_m);
        nextR = nextR * Physics::c * dt_m;

        if (activeSet.empty() && !collectReferenceSamples_m
            && (pathLengthRange_m.isOutside(pathLength_m)
                || (period_m <= 0.0 && stepRange_m.isOutside(currentStep_m)))) {
            errorFlag_m = EOL;
            globalBoundingBox_m.enlargeToContainPosition(r_m);
            return;
        }

    } while (activeSet == itsOpalBeamline_m.getElements(nextR));
}

bool OrbitThreader::reachedThreadingEnd() const {
    const bool stopAtSStop = period_m > 0.0 || collectReferenceSamples_m;
    return stopAtSStop && (dt_m > 0.0 ? pathLength_m >= sStop_m : pathLength_m <= sStop_m);
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
            Vector_t<double, 3> initialR = itsOpalBeamline_m.transformToLocalCS(*it, r_m);
            Vector_t<double, 3> initialP = itsOpalBeamline_m.rotateToLocalCS(*it, p_m);

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

void OrbitThreader::setDesignEnergy(
        ElementList& allElements, const std::set<std::string>& visitedElements) {
    double kineticEnergyeV = reference_m.getM() * (sqrt(dot(p_m, p_m) + 1.0) - 1.0);

    ElementList::iterator it        = allElements.begin();
    const ElementList::iterator end = allElements.end();
    for (; it != end; ++it) {
        std::shared_ptr<ElementBase> element = (*it);
        if (visitedElements.find(element->getName()) == visitedElements.end()
            && !(element->getType() == ElementType::RFCAVITY
                 || element->getType() == ElementType::TRAVELINGWAVE)) {
            element->setDesignEnergy(kineticEnergyeV);
        }
    }
}

void OrbitThreader::computeBoundingBox() {
    ElementList allElements         = itsOpalBeamline_m.getElementByType(ElementType::ANY);
    ElementList::iterator it        = allElements.begin();
    const ElementList::iterator end = allElements.end();
    for (; it != end; ++it) {
        if ((*it)->getType() == ElementType::MARKER) {
            continue;
        }
        BoundingBox other = (*it)->getBoundingBoxInLabCoords();
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
        const std::set<std::shared_ptr<ElementBase>>& elements, const Vector_t<double, 3>& position,
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

void OrbitThreader::recordReferenceSample() {
    LinearTransferMapReference state;
    state = referenceSamples_m.empty()
                    ? LinearTransferMapBuilder::initialFrame(itsOpalBeamline_m, p_m)
                    : LinearTransferMapBuilder::transportFrame(referenceSamples_m.back().state, p_m);
    state.position   = r_m;
    state.momentum   = p_m;
    state.time       = time_m;
    state.pathLength = pathLength_m;
    referenceSamples_m.push_back({state});
}

void OrbitThreader::printCombinedLinearTransferMap() const {
    if (!combinedLinearTransferMap_m || ippl::Comm->rank() != 0) return;
    *gmsg << level1
          << (period_m > 0.0 ? "\n* Combined one-turn linear transfer map"
                             : "\n* Combined linear transfer map")
          << "  (x, x', y, y', zeta, delta):\n";
    const auto& map = *combinedLinearTransferMap_m;
    for (int row = 0; row < 6; ++row) {
        *gmsg << "  ";
        for (int column = 0; column < 6; ++column) {
            *gmsg << std::setw(15) << std::setprecision(7) << std::scientific << map(row, column);
        }
        *gmsg << "\n";
    }
    const double determinantError = *combinedDeterminantResidual_m;
    const double symplecticError  = *combinedSymplecticResidual_m;
    *gmsg << "* Volume-preservation diagnostic: "
          << (determinantError <= mapDiagnosticTolerance ? "PASS" : "FAIL") << "\n"
          << "*   |det(M) - 1| = " << std::setprecision(7) << std::scientific << determinantError
          << "  (tolerance " << mapDiagnosticTolerance << ")\n";
    *gmsg << "* Canonical-form symplecticity diagnostic: "
          << (symplecticError <= mapDiagnosticTolerance ? "PASS" : "FAIL") << "\n"
          << "*   max_ij |(M^T J M - J)_ij| = " << symplecticError << "  (tolerance "
          << mapDiagnosticTolerance << ")\n"
          << "*   The reported slopes and mechanical momenta are not globally canonical;\n"
          << "*   this test is diagnostic unless the entrance/exit coordinates are canonical.\n";
    *gmsg << std::defaultfloat << endl;
}
