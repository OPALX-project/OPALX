//
// Class PlacementResolver
//   The single PLACE stage of the input -> computeExternalFields() pipeline.
//
// Copyright (c) 2024, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#include "Elements/PlacementResolver.h"

#include "AbsBeamline/ElementBase.h"
#include "Algorithms/Quaternion.hpp"
#include "Utilities/Options.h"

#include <cmath>
#include <vector>

void PlacementResolver::resolve(FieldList& elements, const CoordinateSystemTrafo& labFrame) {
    static unsigned int order     = 0;
    const FieldList::iterator end = elements.end();

    // Phase 1 — 6D-pose (Mode A) elements: compose the recorded global-to-local pose with the
    // lab frame and fix them in place. The reference-path walk below (Mode B / ELEMEDGE) then
    // skips them via isPositioned(). Both modes are resolved here, in one pass.
    for (FieldList::iterator it = elements.begin(); it != end; ++it) {
        const std::shared_ptr<ElementBase> element = (*it).getElement();
        if (element->isElementPositionSet()) {
            continue;  // Mode B (ELEMEDGE): placed by the reference-path walk below
        }
        CoordinateSystemTrafo toElement = element->getCSTrafoGlobal2Local();
        toElement *= labFrame;
        element->setCSTrafoGlobal2Local(toElement);
        element->fixPosition();
    }

    unsigned int minOrder = order;
    {
        double endPriorPathLength               = 0.0;
        CoordinateSystemTrafo currentCoordTrafo = labFrame;

        FieldList::iterator it = elements.begin();
        for (; it != end; ++it) {
            std::shared_ptr<ElementBase> element = (*it).getElement();
            if (element->isPositioned()) {
                continue;
            }
            (*it).order_m = minOrder;

            if (element->getType() != ElementType::SBEND && element->getType() != ElementType::RBEND
                && element->getType() != ElementType::RBEND3D) {
                continue;
            }

            double beginThisPathLength = element->getElementPosition();
            Vector_t<double, 3> beginThis3D(0, 0, beginThisPathLength - endPriorPathLength);
            double thisLength    = element->getChordLength();
            double bendAngle     = element->getBendAngle();
            double entranceAngle = element->getEntranceAngle();
            double arcLength     = element->getArcLength();

            double rotationAngleAboutZ = element->getRotationAboutZ();
            Quaternion_t rotationAboutZ(
                    cos(0.5 * rotationAngleAboutZ),
                    sin(-0.5 * rotationAngleAboutZ) * Vector_t<double, 3>(0, 0, 1));

            Vector_t<double, 3> effectiveRotationAxis =
                    rotationAboutZ.rotate(Vector_t<double, 3>(0, -1, 0));
            effectiveRotationAxis = effectiveRotationAxis / euclidean_norm(effectiveRotationAxis);

            Quaternion_t rotationAboutAxis(
                    cos(0.5 * bendAngle), sin(0.5 * bendAngle) * effectiveRotationAxis);
            Quaternion_t halfRotationAboutAxis(
                    cos(0.25 * bendAngle), sin(0.25 * bendAngle) * effectiveRotationAxis);
            Quaternion_t entryFaceRotation(
                    cos(0.5 * entranceAngle), sin(0.5 * entranceAngle) * effectiveRotationAxis);

            if (!Options::idealized) {
                std::vector<Vector_t<double, 3>> truePath = element->getDesignPath();
                Quaternion_t directionExitHardEdge(
                        cos(0.5 * (0.5 * bendAngle - entranceAngle)),
                        sin(0.5 * (0.5 * bendAngle - entranceAngle)) * effectiveRotationAxis);
                Vector_t<double, 3> exitHardEdge =
                        thisLength * directionExitHardEdge.rotate(Vector_t<double, 3>(0, 0, 1));
                double distanceEntryHETruePath = euclidean_norm(truePath.front());
                Vector_t<double, 3> exitDelta =
                        rotationAboutZ.rotate(truePath.back()) - exitHardEdge;
                double distanceExitHETruePath = euclidean_norm(exitDelta);
                // Bend path length = body length (end-start of the field interval, which the
                // BeamlineFieldElement stored as startField + getElementLength()).
                double pathLengthTruePath = element->getElementLength();
                arcLength = pathLengthTruePath - distanceEntryHETruePath - distanceExitHETruePath;
            }

            Vector_t<double, 3> chord =
                    thisLength * halfRotationAboutAxis.rotate(Vector_t<double, 3>(0, 0, 1));
            Vector_t<double, 3> endThis3D = beginThis3D + chord;
            double endThisPathLength      = beginThisPathLength + arcLength;

            CoordinateSystemTrafo fromEndLastToBeginThis(
                    beginThis3D, (entryFaceRotation * rotationAboutZ).conjugate());
            CoordinateSystemTrafo fromEndLastToEndThis(endThis3D, rotationAboutAxis.conjugate());

            element->setCSTrafoGlobal2Local(fromEndLastToBeginThis * currentCoordTrafo);

            currentCoordTrafo = (fromEndLastToEndThis * currentCoordTrafo);

            endPriorPathLength = endThisPathLength;
        }
    }

    double endPriorPathLength               = 0.0;
    CoordinateSystemTrafo currentCoordTrafo = labFrame;

    FieldList::iterator it = elements.begin();
    for (; it != end; ++it) {
        std::shared_ptr<ElementBase> element = (*it).getElement();
        if (element->isPositioned()) continue;

        (*it).order_m = order++;

        double beginThisPathLength = element->getElementPosition();
        double thisLength          = element->getElementLength();
        Vector_t<double, 3> beginThis3D(0, 0, beginThisPathLength - endPriorPathLength);

        if (element->getType() == ElementType::SOURCE) {
            beginThis3D(2) -= thisLength;
        }

        Vector_t<double, 3> endThis3D;
        if (element->getType() == ElementType::SBEND || element->getType() == ElementType::RBEND
            || element->getType() == ElementType::RBEND3D) {
            thisLength       = element->getChordLength();
            double bendAngle = element->getBendAngle();

            double rotationAngleAboutZ = element->getRotationAboutZ();
            Quaternion_t rotationAboutZ(
                    cos(0.5 * rotationAngleAboutZ),
                    sin(-0.5 * rotationAngleAboutZ) * Vector_t<double, 3>(0, 0, 1));

            Vector_t<double, 3> effectiveRotationAxis =
                    rotationAboutZ.rotate(Vector_t<double, 3>(0, -1, 0));
            effectiveRotationAxis = effectiveRotationAxis / euclidean_norm(effectiveRotationAxis);

            Quaternion_t rotationAboutAxis(
                    cos(0.5 * bendAngle), sin(0.5 * bendAngle) * effectiveRotationAxis);
            Quaternion halfRotationAboutAxis(
                    cos(0.25 * bendAngle), sin(0.25 * bendAngle) * effectiveRotationAxis);

            double arcLength = element->getArcLength();
            if (!Options::idealized) {
                std::vector<Vector_t<double, 3>> truePath = element->getDesignPath();
                double entranceAngle                      = element->getEntranceAngle();
                Quaternion_t directionExitHardEdge(
                        cos(0.5 * (0.5 * bendAngle - entranceAngle)),
                        sin(0.5 * (0.5 * bendAngle - entranceAngle)) * effectiveRotationAxis);
                Vector_t<double, 3> exitHardEdge =
                        thisLength * directionExitHardEdge.rotate(Vector_t<double, 3>(0, 0, 1));
                double distanceEntryHETruePath = euclidean_norm(truePath.front());
                Vector_t<double, 3> exitDelta =
                        rotationAboutZ.rotate(truePath.back()) - exitHardEdge;
                double distanceExitHETruePath = euclidean_norm(exitDelta);
                // Bend path length = body length (end-start of the field interval, which the
                // BeamlineFieldElement stored as startField + getElementLength()).
                double pathLengthTruePath = element->getElementLength();
                arcLength = pathLengthTruePath - distanceEntryHETruePath - distanceExitHETruePath;
            }

            endThis3D =
                    (beginThis3D
                     + halfRotationAboutAxis.rotate(Vector_t<double, 3>(0, 0, thisLength)));
            CoordinateSystemTrafo fromEndLastToEndThis(endThis3D, rotationAboutAxis.conjugate());
            currentCoordTrafo = fromEndLastToEndThis * currentCoordTrafo;

            endPriorPathLength = beginThisPathLength + arcLength;
        } else {
            double rotationAngleAboutZ = (*it).getElement()->getRotationAboutZ();
            Quaternion_t rotationAboutZ(
                    cos(0.5 * rotationAngleAboutZ),
                    sin(-0.5 * rotationAngleAboutZ) * Vector_t<double, 3>(0, 0, 1));

            CoordinateSystemTrafo fromLastToThis(beginThis3D, rotationAboutZ);

            element->setCSTrafoGlobal2Local(fromLastToThis * currentCoordTrafo);
        }

        element->fixPosition();
    }
}
