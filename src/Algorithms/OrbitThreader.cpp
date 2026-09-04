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
    constexpr std::array<double, 6> transferMapSteps{1.0e-3, 1.0e-3, 1.0e-3,
                                                     1.0e-3, 1.0e-3, 1.0e-3};
    constexpr double boundaryTolerance      = 1.0e-12;
    constexpr double mapDiagnosticTolerance = 1.0e-6;

    Vector_t<double, 3> crossProduct(const Vector_t<double, 3>& a, const Vector_t<double, 3>& b) {
        return Vector_t<double, 3>(
                a(1) * b(2) - a(2) * b(1), a(2) * b(0) - a(0) * b(2), a(0) * b(1) - a(1) * b(0));
    }

    Vector_t<double, 3> normalized(const Vector_t<double, 3>& value) {
        const double norm = euclidean_norm(value);
        if (!(norm > 0.0)) {
            throw OpalException("OrbitThreader::normalized", "Cannot normalize a zero vector.");
        }
        return value / norm;
    }

    Vector_t<double, 3> rotateRodrigues(
            const Vector_t<double, 3>& value, const Vector_t<double, 3>& axis, const double angle) {
        return value * std::cos(angle) + crossProduct(axis, value) * std::sin(angle)
               + axis * dot(axis, value) * (1.0 - std::cos(angle));
    }

    double invert6x6(const matrix6x6_t& input, matrix6x6_t& inverse) {
        double augmented[6][12]{};
        double normInput = 0.0;
        for (int row = 0; row < 6; ++row) {
            double rowSum = 0.0;
            for (int col = 0; col < 6; ++col) {
                augmented[row][col]     = input(row, col);
                augmented[row][col + 6] = row == col ? 1.0 : 0.0;
                rowSum += std::abs(input(row, col));
            }
            normInput = std::max(normInput, rowSum);
        }

        for (int col = 0; col < 6; ++col) {
            int pivot = col;
            for (int row = col + 1; row < 6; ++row) {
                if (std::abs(augmented[row][col]) > std::abs(augmented[pivot][col])) pivot = row;
            }
            if (std::abs(augmented[pivot][col]) < 1.0e-18) {
                throw OpalException(
                        "OrbitThreader::makeLinearTransferMap",
                        "The finite-difference input matrix is singular.");
            }
            if (pivot != col) {
                for (int j = 0; j < 12; ++j)
                    std::swap(augmented[pivot][j], augmented[col][j]);
            }
            const double diagonal = augmented[col][col];
            for (int j = 0; j < 12; ++j)
                augmented[col][j] /= diagonal;
            for (int row = 0; row < 6; ++row) {
                if (row == col) continue;
                const double factor = augmented[row][col];
                for (int j = 0; j < 12; ++j)
                    augmented[row][j] -= factor * augmented[col][j];
            }
        }

        double normInverse = 0.0;
        for (int row = 0; row < 6; ++row) {
            double rowSum = 0.0;
            for (int col = 0; col < 6; ++col) {
                inverse(row, col) = augmented[row][col + 6];
                rowSum += std::abs(inverse(row, col));
            }
            normInverse = std::max(normInverse, rowSum);
        }
        return normInput * normInverse;
    }

    double symplecticResidual(const matrix6x6_t& map) {
        matrix6x6_t j(0.0);
        for (int pair = 0; pair < 3; ++pair) {
            j(2 * pair, 2 * pair + 1) = 1.0;
            j(2 * pair + 1, 2 * pair) = -1.0;
        }
        const matrix6x6_t residual = prod(get_transpose(map), prod(j, map));
        double maximum             = 0.0;
        for (int row = 0; row < 6; ++row) {
            for (int col = 0; col < 6; ++col) {
                maximum = std::max(maximum, std::abs(residual(row, col) - j(row, col)));
            }
        }
        return maximum;
    }

    double determinant(const matrix6x6_t& map) {
        double work[6][6]{};
        for (int row = 0; row < 6; ++row) {
            for (int col = 0; col < 6; ++col) {
                work[row][col] = map(row, col);
            }
        }

        double result = 1.0;
        for (int col = 0; col < 6; ++col) {
            int pivot = col;
            for (int row = col + 1; row < 6; ++row) {
                if (std::abs(work[row][col]) > std::abs(work[pivot][col])) pivot = row;
            }
            if (work[pivot][col] == 0.0) return 0.0;
            if (pivot != col) {
                for (int j = col; j < 6; ++j) {
                    std::swap(work[pivot][j], work[col][j]);
                }
                result = -result;
            }
            const double diagonal = work[col][col];
            result *= diagonal;
            for (int row = col + 1; row < 6; ++row) {
                const double factor = work[row][col] / diagonal;
                for (int j = col + 1; j < 6; ++j) {
                    work[row][j] -= factor * work[col][j];
                }
            }
        }
        return result;
    }

    double determinantResidual(const matrix6x6_t& map) { return std::abs(determinant(map) - 1.0); }
}  // namespace

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
      reference_m(ref) {
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
        transferMapSegments_m.clear();
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

        const double finalS = reachedPeriodicEnd() ? sStop_m : pathLength_m;
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
    } while (errorFlag_m != EOL && (period_m > 0.0 || stepRange_m.isInside(currentStep_m))
             && !(pathLengthRange_m.isOutside(pathLength_m) && intersection.empty()
                  && !(elementSet.empty() || currentSet.empty())));

    imap_m.tidyUp(sStop_m);
    collectReferenceSamples_m = false;
    if (calculateMaps) {
        calculateLinearTransferMaps();
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

        IndexMap::value_t::const_iterator it        = activeSet.begin();
        const IndexMap::value_t::const_iterator end = activeSet.end();
        Vector_t<double, 3> oldR                    = r_m;
        const Vector_t<double, 3> oldP              = p_m;
        const double oldTime                        = time_m;
        const double oldPathLength                  = pathLength_m;

        r_m /= Physics::c * dt_m;
        integrator_m.push(r_m, p_m, dt_m);
        r_m = r_m * Physics::c * dt_m;

        Vector_t<double, 3> Ef(0.0), Bf(0.0);
        std::string names("\t");
        for (; it != end; ++it) {
            Vector_t<double, 3> localR = itsOpalBeamline_m.transformToLocalCS(*it, r_m);
            Vector_t<double, 3> localP = itsOpalBeamline_m.rotateToLocalCS(*it, p_m);
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

        r_m /= Physics::c * dt_m;
        integrator_m.kick(r_m, p_m, Ef, Bf, dt_m, reference_m.getM(), reference_m.getQ());
        integrator_m.push(r_m, p_m, dt_m);
        r_m = r_m * Physics::c * dt_m;

        const Vector<double, 3> d = r_m - oldR;

        pathLength_m += std::copysign(euclidean_norm(d), dt_m);

        ++currentStep_m;
        time_m += dt_m;
        if (collectReferenceSamples_m) {
            if (reachedPeriodicEnd() && pathLength_m != oldPathLength) {
                const double fraction  = (sStop_m - oldPathLength) / (pathLength_m - oldPathLength);
                const RayState clipped = advanceRay({oldR, oldP, oldTime}, fraction * dt_m);
                LinearTransferMapReference state =
                        transportFrame(referenceSamples_m.back().state, clipped.momentum);
                state.position   = clipped.position;
                state.momentum   = clipped.momentum;
                state.time       = clipped.time;
                state.pathLength = sStop_m;
                referenceSamples_m.push_back({state});
            } else {
                recordReferenceSample();
            }
        }

        if (reachedPeriodicEnd()) {
            errorFlag_m = EOL;
            globalBoundingBox_m.enlargeToContainPosition(r_m);
            return;
        }

        nextR = r_m / (Physics::c * dt_m);
        integrator_m.push(nextR, p_m, dt_m);
        nextR = nextR * Physics::c * dt_m;

        if (activeSet.empty()
            && (pathLengthRange_m.isOutside(pathLength_m)
                || (period_m <= 0.0 && stepRange_m.isOutside(currentStep_m)))) {
            errorFlag_m = EOL;
            globalBoundingBox_m.enlargeToContainPosition(r_m);
            return;
        }

    } while (activeSet == itsOpalBeamline_m.getElements(nextR));
}

bool OrbitThreader::reachedPeriodicEnd() const {
    return period_m > 0.0 && (dt_m > 0.0 ? pathLength_m >= sStop_m : pathLength_m <= sStop_m);
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

LinearTransferMapReference OrbitThreader::transportFrame(
        const LinearTransferMapReference& frame, const Vector_t<double, 3>& momentum) {
    LinearTransferMapReference result  = frame;
    const Vector_t<double, 3> newS     = normalized(momentum);
    const Vector_t<double, 3> rotation = crossProduct(frame.sAxis, newS);
    const double sine                  = euclidean_norm(rotation);
    const double cosine                = std::max(-1.0, std::min(1.0, dot(frame.sAxis, newS)));

    Vector_t<double, 3> newX = frame.xAxis;
    if (sine > 1.0e-14) {
        const Vector_t<double, 3> axis = rotation / sine;
        newX = rotateRodrigues(frame.xAxis, axis, std::atan2(sine, cosine));
    } else if (cosine < 0.0) {
        // A 180-degree reversal has no unique minimum-rotation axis.  Preserve the old
        // transverse x direction, which is already perpendicular to both tangents.
        newX = frame.xAxis;
    }

    newX         = normalized(newX - newS * dot(newX, newS));
    result.sAxis = newS;
    result.xAxis = newX;
    result.yAxis = normalized(crossProduct(newS, newX));
    return result;
}

void OrbitThreader::recordReferenceSample() {
    LinearTransferMapReference state;
    if (referenceSamples_m.empty()) {
        state.sAxis               = normalized(p_m);
        Vector_t<double, 3> xAxis = itsOpalBeamline_m.getCSTrafoLab2Local().rotateFrom(
                Vector_t<double, 3>(1.0, 0.0, 0.0));
        xAxis -= state.sAxis * dot(xAxis, state.sAxis);
        if (euclidean_norm(xAxis) < 1.0e-12) {
            xAxis = itsOpalBeamline_m.getCSTrafoLab2Local().rotateFrom(
                    Vector_t<double, 3>(0.0, 1.0, 0.0));
            xAxis -= state.sAxis * dot(xAxis, state.sAxis);
        }
        state.xAxis = normalized(xAxis);
        state.yAxis = normalized(crossProduct(state.sAxis, state.xAxis));
    } else {
        state = transportFrame(referenceSamples_m.back().state, p_m);
    }
    state.position   = r_m;
    state.momentum   = p_m;
    state.time       = time_m;
    state.pathLength = pathLength_m;
    referenceSamples_m.push_back({state});
}

OrbitThreader::RayState OrbitThreader::advanceRay(const RayState& state, const double dt) {
    RayState result                    = state;
    Vector_t<double, 3> scaledPosition = result.position / (Physics::c * dt);
    integrator_m.push(scaledPosition, result.momentum, dt);
    result.position = scaledPosition * Physics::c * dt;

    Vector_t<double, 3> electric(0.0), magnetic(0.0);
    itsOpalBeamline_m.getFieldAt(
            result.position, result.momentum, result.time + 0.5 * dt, electric, magnetic);

    scaledPosition = result.position / (Physics::c * dt);
    integrator_m.kick(
            scaledPosition, result.momentum, electric, magnetic, dt, reference_m.getM(),
            reference_m.getQ());
    integrator_m.push(scaledPosition, result.momentum, dt);
    result.position = scaledPosition * Physics::c * dt;
    result.time += dt;
    return result;
}

LinearTransferMapReference OrbitThreader::refineBoundary(
        const std::shared_ptr<ElementBase>& element, const LinearTransferMapReference& before,
        const LinearTransferMapReference& after, const bool entering) {
    const auto isInside = [&](const RayState& ray) {
        return element->isInside(itsOpalBeamline_m.transformToLocalCS(element, ray.position));
    };

    const RayState start{before.position, before.momentum, before.time};
    double lower   = 0.0;
    double upper   = after.time - before.time;
    RayState trial = start;
    for (int iteration = 0; iteration < 60; ++iteration) {
        const double middle = 0.5 * (lower + upper);
        trial               = advanceRay(start, middle);
        if (isInside(trial) == entering) {
            upper = middle;
        } else {
            lower = middle;
        }
        if (std::abs(upper - lower) <= boundaryTolerance * std::abs(dt_m)) break;
    }
    const double fraction               = (after.time != before.time)
                                                  ? (trial.time - before.time) / (after.time - before.time)
                                                  : 0.0;
    LinearTransferMapReference boundary = transportFrame(before, trial.momentum);
    boundary.position                   = trial.position;
    boundary.momentum                   = trial.momentum;
    boundary.time                       = trial.time;
    boundary.pathLength = before.pathLength + fraction * (after.pathLength - before.pathLength);
    return boundary;
}

std::array<double, 6> OrbitThreader::coordinates(
        const RayState& ray, const LinearTransferMapReference& reference) {
    const Vector_t<double, 3> displacement = ray.position - reference.position;
    const double longitudinalMomentum      = dot(ray.momentum, reference.sAxis);
    if (std::abs(longitudinalMomentum) < 1.0e-14) {
        throw OpalException(
                "OrbitThreader::coordinates", "A shadow ray has zero longitudinal momentum.");
    }
    const double referenceMomentum = euclidean_norm(reference.momentum);
    const double beta = referenceMomentum / std::sqrt(1.0 + referenceMomentum * referenceMomentum);
    return {dot(displacement, reference.xAxis),
            dot(ray.momentum, reference.xAxis) / longitudinalMomentum,
            dot(displacement, reference.yAxis),
            dot(ray.momentum, reference.yAxis) / longitudinalMomentum,
            -beta * Physics::c * (ray.time - reference.time),
            euclidean_norm(ray.momentum) / referenceMomentum - 1.0};
}

OrbitThreader::RayState OrbitThreader::rayFromCoordinates(
        const std::array<double, 6>& coordinate, const LinearTransferMapReference& reference) {
    RayState ray;
    ray.position =
            reference.position + coordinate[0] * reference.xAxis + coordinate[2] * reference.yAxis;
    const double referenceMomentum = euclidean_norm(reference.momentum);
    const double momentum          = referenceMomentum * (1.0 + coordinate[5]);
    Vector_t<double, 3> direction =
            coordinate[1] * reference.xAxis + coordinate[3] * reference.yAxis + reference.sAxis;
    ray.momentum      = momentum * normalized(direction);
    const double beta = referenceMomentum / std::sqrt(1.0 + referenceMomentum * referenceMomentum);
    ray.time          = reference.time - coordinate[4] / (beta * Physics::c);
    return ray;
}

OrbitThreader::RayState OrbitThreader::trackRayToExit(
        const RayState& initial, const LinearTransferMapReference& exit,
        const double referenceFlightTime) {
    RayState previous                        = initial;
    Vector_t<double, 3> previousDisplacement = previous.position - exit.position;
    double previousDistance                  = dot(previousDisplacement, exit.sAxis);
    const double expectedDirection           = dot(exit.momentum, exit.sAxis) >= 0.0 ? 1.0 : -1.0;
    const auto crossed = [&](const double oldDistance, const double newDistance) {
        return expectedDirection * oldDistance <= 0.0 && expectedDirection * newDistance >= 0.0;
    };

    const std::size_t referenceSteps =
            static_cast<std::size_t>(std::ceil(std::abs(referenceFlightTime / dt_m)));
    const std::size_t maximumSteps = std::max<std::size_t>(100, 4 * referenceSteps + 100);
    for (std::size_t step = 0; step < maximumSteps; ++step) {
        RayState next                              = advanceRay(previous, dt_m);
        const Vector_t<double, 3> nextDisplacement = next.position - exit.position;
        const double nextDistance                  = dot(nextDisplacement, exit.sAxis);
        if (crossed(previousDistance, nextDistance)) {
            double lower   = 0.0;
            double upper   = dt_m;
            RayState trial = next;
            for (int iteration = 0; iteration < 60; ++iteration) {
                const double middle                    = 0.5 * (lower + upper);
                trial                                  = advanceRay(previous, middle);
                const Vector_t<double, 3> displacement = trial.position - exit.position;
                const double distance                  = dot(displacement, exit.sAxis);
                if (crossed(previousDistance, distance)) {
                    upper = middle;
                } else {
                    lower = middle;
                }
                if (std::abs(upper - lower) <= boundaryTolerance * std::abs(dt_m)) break;
            }
            return trial;
        }
        previous         = next;
        previousDistance = nextDistance;
    }
    throw OpalException(
            "OrbitThreader::trackRayToExit",
            "A transfer-map shadow ray did not reach the recorded exit plane.");
}

LinearTransferMap OrbitThreader::makeLinearTransferMap(
        const LinearTransferMapReference& entrance, const LinearTransferMapReference& exit,
        const std::size_t pass) {
    matrix6x6_t inputDifferences(0.0);
    matrix6x6_t outputDifferences(0.0);
    for (int column = 0; column < 6; ++column) {
        std::array<double, 6> plus{};
        std::array<double, 6> minus{};
        plus[column]                 = transferMapSteps[column];
        minus[column]                = -transferMapSteps[column];
        const RayState plusEntrance  = rayFromCoordinates(plus, entrance);
        const RayState minusEntrance = rayFromCoordinates(minus, entrance);
        const RayState plusExit  = trackRayToExit(plusEntrance, exit, exit.time - entrance.time);
        const RayState minusExit = trackRayToExit(minusEntrance, exit, exit.time - entrance.time);
        const auto encodedPlusEntrance  = coordinates(plusEntrance, entrance);
        const auto encodedMinusEntrance = coordinates(minusEntrance, entrance);
        const auto encodedPlusExit      = coordinates(plusExit, exit);
        const auto encodedMinusExit     = coordinates(minusExit, exit);
        for (int row = 0; row < 6; ++row) {
            inputDifferences(row, column) =
                    0.5 * (encodedPlusEntrance[row] - encodedMinusEntrance[row]);
            outputDifferences(row, column) = 0.5 * (encodedPlusExit[row] - encodedMinusExit[row]);
        }
    }

    matrix6x6_t inverseInput(0.0);
    const double condition = invert6x6(inputDifferences, inverseInput);
    if (!std::isfinite(condition) || condition > 1.0e12) {
        throw OpalException(
                "OrbitThreader::makeLinearTransferMap",
                "The finite-difference input matrix is ill-conditioned (condition number "
                        + std::to_string(condition) + ").");
    }

    LinearTransferMap result;
    result.matrix                    = prod(outputDifferences, inverseInput);
    result.finiteDifferenceSteps     = transferMapSteps;
    result.entrance                  = entrance;
    result.exit                      = exit;
    result.pass                      = pass;
    result.inputConditionNumber      = condition;
    result.determinantResidual       = determinantResidual(result.matrix);
    result.symplecticResidual        = symplecticResidual(result.matrix);
    result.includesOverlappingFields = false;
    return result;
}

void OrbitThreader::calculateLinearTransferMaps() {
    if (referenceSamples_m.size() < 2) return;

    // OrbitThreader first tracks backwards by the bunch extent.  That path is needed to build
    // the IndexMap, but a LINE/one-turn transfer map must start at the requested reference s.
    const auto firstAtOrAfterStart = std::lower_bound(
            referenceSamples_m.begin(), referenceSamples_m.end(), transferMapStartPathLength_m,
            [](const ReferenceSample& sample, const double pathLength) {
                return sample.state.pathLength < pathLength;
            });
    if (firstAtOrAfterStart == referenceSamples_m.end()) return;
    if (firstAtOrAfterStart != referenceSamples_m.begin()
        && firstAtOrAfterStart->state.pathLength > transferMapStartPathLength_m) {
        const auto& before    = std::prev(firstAtOrAfterStart)->state;
        const auto& after     = firstAtOrAfterStart->state;
        const double fraction = (transferMapStartPathLength_m - before.pathLength)
                                / (after.pathLength - before.pathLength);
        const RayState clipped = advanceRay(
                {before.position, before.momentum, before.time},
                fraction * (after.time - before.time));
        LinearTransferMapReference start = transportFrame(before, clipped.momentum);
        start.position                   = clipped.position;
        start.momentum                   = clipped.momentum;
        start.time                       = clipped.time;
        start.pathLength                 = transferMapStartPathLength_m;

        std::vector<ReferenceSample> clippedSamples;
        clippedSamples.reserve(
                1 + static_cast<std::size_t>(referenceSamples_m.end() - firstAtOrAfterStart));
        clippedSamples.push_back({start});
        clippedSamples.insert(clippedSamples.end(), firstAtOrAfterStart, referenceSamples_m.end());
        referenceSamples_m = std::move(clippedSamples);
    } else {
        referenceSamples_m.erase(referenceSamples_m.begin(), firstAtOrAfterStart);
    }
    if (referenceSamples_m.size() < 2) return;

    const auto validateActiveSet = [](const IndexMap::value_t& active) {
        for (const auto& element : active) {
            if (element->getType() == ElementType::RFCAVITY
                || element->getType() == ElementType::TRAVELINGWAVE) {
                throw OpalException(
                        "OrbitThreader::calculateLinearTransferMaps",
                        "RF elements are not supported by the first linear-transfer-map "
                        "implementation (element '"
                                + element->getName() + "').");
            }
        }
    };

    using ElementPtr = std::shared_ptr<ElementBase>;
    struct BoundaryEvent {
        LinearTransferMapReference state;
        ElementPtr element;
        bool entering;
    };
    struct Segment {
        LinearTransferMapReference entrance;
        LinearTransferMapReference exit;
        IndexMap::value_t active;
    };

    std::vector<BoundaryEvent> events;
    auto previousSet = itsOpalBeamline_m.getElements(referenceSamples_m.front().state.position);
    validateActiveSet(previousSet);

    for (std::size_t sample = 1; sample < referenceSamples_m.size(); ++sample) {
        const auto& before = referenceSamples_m[sample - 1].state;
        const auto& after  = referenceSamples_m[sample].state;
        auto currentSet    = itsOpalBeamline_m.getElements(after.position);
        validateActiveSet(currentSet);

        IndexMap::value_t exited;
        std::set_difference(
                previousSet.begin(), previousSet.end(), currentSet.begin(), currentSet.end(),
                std::inserter(exited, exited.begin()));
        for (const auto& element : exited) {
            events.push_back({refineBoundary(element, before, after, false), element, false});
        }

        IndexMap::value_t entered;
        std::set_difference(
                currentSet.begin(), currentSet.end(), previousSet.begin(), previousSet.end(),
                std::inserter(entered, entered.begin()));
        for (const auto& element : entered) {
            events.push_back({refineBoundary(element, before, after, true), element, true});
        }
        previousSet = std::move(currentSet);
    }

    std::sort(events.begin(), events.end(), [](const auto& left, const auto& right) {
        if (left.state.pathLength != right.state.pathLength) {
            return left.state.pathLength < right.state.pathLength;
        }
        // At an exactly shared boundary, close the old set before opening the new one.
        if (left.entering != right.entering) return left.entering < right.entering;
        return left.element->getName() < right.element->getName();
    });

    std::vector<Segment> segments;
    IndexMap::value_t active =
            itsOpalBeamline_m.getElements(referenceSamples_m.front().state.position);
    LinearTransferMapReference segmentEntrance = referenceSamples_m.front().state;
    constexpr double minimumSegmentLength      = 1.0e-12;
    for (const auto& event : events) {
        if (event.state.pathLength - segmentEntrance.pathLength > minimumSegmentLength) {
            segments.push_back({segmentEntrance, event.state, active});
        }
        if (event.entering) {
            active.insert(event.element);
        } else {
            active.erase(event.element);
        }
        segmentEntrance = event.state;
    }
    const auto& finalState = referenceSamples_m.back().state;
    if (finalState.pathLength - segmentEntrance.pathLength > minimumSegmentLength) {
        segments.push_back({segmentEntrance, finalState, active});
    }

    std::map<ElementPtr, std::size_t, std::owner_less<ElementPtr>> passes;
    transferMapSegments_m.reserve(segments.size());
    for (std::size_t segment = 0; segment < segments.size(); ++segment) {
        const auto& interval = segments[segment];
        validateActiveSet(interval.active);
        LinearTransferMap map = makeLinearTransferMap(interval.entrance, interval.exit, segment);
        map.segment           = segment;
        map.includesOverlappingFields = interval.active.size() > 1;
        for (const auto& element : interval.active) {
            map.activeElements.push_back(element->getName());
        }
        std::sort(map.activeElements.begin(), map.activeElements.end());
        transferMapSegments_m.push_back(map);
        for (const auto& element : interval.active) {
            LinearTransferMap attached = map;
            attached.pass              = passes[element]++;
            element->addLinearTransferMap(std::move(attached));
        }
    }

    if (transferMapSegments_m.empty()) return;
    matrix6x6_t combined;
    for (const auto& segment : transferMapSegments_m) {
        combined = prod(segment.matrix, combined);
    }
    combinedLinearTransferMap_m   = combined;
    combinedDeterminantResidual_m = determinantResidual(combined);
    combinedSymplecticResidual_m  = symplecticResidual(combined);
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
