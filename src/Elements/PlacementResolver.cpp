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
#include "BeamlineGeometry/Geometry.h"
#include "Utilities/OpalException.h"

#include <cmath>

/**
 * @brief Give every element its global -> local (lab -> entrance) transform.
 *
 * Two placement modes are resolved here:
 *   - Mode A (explicit 6D pose): X/Y/Z/PHI/PSI/THETA already are the entrance
 *     frame, so just compose them with the lab frame (Phase 1).
 *   - Mode B (ELEMEDGE): the element sits at a path length along the reference
 *     orbit; its frame comes from walking that orbit (Phases 2 and 3).
 *
 * The walk only changes orientation at bends; a straight element keeps the
 * previous orientation and its spacing is carried by the ELEMEDGE path lengths.
 * Phase 2 sets each bend's entrance frame; Phase 3 places every remaining
 * element and fixes all positions.
 */
void PlacementResolver::resolve(ElementList& elements, const CoordinateSystemTrafo& labFrame) {
    const ElementList::iterator end = elements.end();

    // Guard: the whole beamline must use ONE placement convention. Either every element is placed
    // by ELEMEDGE (Mode B) or every element by a 6D lab-frame pose (Mode A). Mixing them is
    // rejected because a 6D-posed bend does not deflect the reference orbit that the ELEMEDGE walk
    // (Phases 2 and 3) follows, so every ELEMEDGE element downstream of a positioned bend would be
    // silently misplaced. The two modes are otherwise resolved independently below.
    {
        std::shared_ptr<ElementBase> elemEdgeCheck;  // an element placed by ELEMEDGE (Mode B)
        std::shared_ptr<ElementBase> poseCheck;      // an element placed by 6D pose (Mode A)
        for (ElementList::iterator it = elements.begin(); it != end; ++it) {
            if ((*it)->isElementPositionSet()) {
                if (!elemEdgeCheck) {
                    elemEdgeCheck = *it;
                }
            } else if (!poseCheck) {
                poseCheck = *it;
            }
        }
        if (elemEdgeCheck && poseCheck) {
            throw OpalException(
                    "PlacementResolver::resolve",
                    "beamline mixes placement conventions: \"" + elemEdgeCheck->getName()
                            + "\" is placed by ELEMEDGE while \"" + poseCheck->getName()
                            + "\" is placed by a 6D lab-frame pose (X/Y/Z/PHI/PSI/THETA). Use one "
                              "convention for the whole beamline: ELEMEDGE for every element or a "
                              "6D pose for every element.");
        }
    }

    // Phase 1 — 6D-pose (Mode A) elements: compose the recorded global-to-local pose with the
    // lab frame and fix them in place. The reference-path walk below (Mode B / ELEMEDGE) then
    // skips them via isPositioned(). Both modes are resolved here, in one pass.
    for (ElementList::iterator it = elements.begin(); it != end; ++it) {
        const std::shared_ptr<ElementBase> element = (*it);
        if (element->isElementPositionSet()) {
            continue;  // Mode B (ELEMEDGE): placed by the reference-path walk below
        }
        // Mode-A X/Y/Z/PHI/PSI/THETA IS the element's geometrical ENTRANCE frame (the frame the
        // tracker and field kernels use), exactly like ELEMEDGE placement. Just compose it with
        // the lab frame; no body-centre-to-entrance correction is needed.
        CoordinateSystemTrafo toElement = element->getCSTrafoGlobal2Local();
        toElement *= labFrame;
        element->setCSTrafoGlobal2Local(toElement);
        element->fixPosition();
    }

    // Phase 2 — set the entrance frame of each ELEMEDGE bend by walking the reference orbit.
    {
        // Running path length and frame at the exit of the previous bend.
        double endPriorPathLength               = 0.0;
        CoordinateSystemTrafo currentCoordTrafo = labFrame;

        ElementList::iterator it = elements.begin();
        for (; it != end; ++it) {
            std::shared_ptr<ElementBase> element = (*it);
            if (element->isPositioned()) {
                continue;  // already placed by Phase 1 (Mode A)
            }

            // This pass only frames bends; every other element is placed in Phase 3.
            if (element->getType() != ElementType::SBEND && element->getType() != ElementType::RBEND
                && element->getType() != ElementType::RBEND3D) {
                continue;
            }

            // Entrance offset from the previous bend's exit, along the straight run between them.
            double beginThisPathLength = element->getElementPosition();
            Vector_t<double, 3> beginThis3D(0, 0, beginThisPathLength - endPriorPathLength);
            double thisLength = element->getGeometry().getChordLength();
            double bendAngle  = element->getGeometry().getBendAngle();
            double arcLength  = element->getGeometry().getArcLength();

            // Roll about z: tilts the bend plane (e.g. a vertical bend rolls the axis out of x-z).
            double rotationAngleAboutZ = element->getRotationAboutZ();
            Quaternion_t rotationAboutZ(
                    cos(0.5 * rotationAngleAboutZ),
                    sin(-0.5 * rotationAngleAboutZ) * Vector_t<double, 3>(0, 0, 1));

            // Bend axis: the nominal -y axis, rolled about z. The orbit turns about this axis.
            Vector_t<double, 3> effectiveRotationAxis =
                    rotationAboutZ.rotate(Vector_t<double, 3>(0, -1, 0));
            effectiveRotationAxis = effectiveRotationAxis / euclidean_norm(effectiveRotationAxis);

            // Turn by the full bend angle (entrance -> exit tangent) and by half of it
            // (entrance tangent -> chord).
            Quaternion_t rotationAboutAxis(
                    cos(0.5 * bendAngle), sin(0.5 * bendAngle) * effectiveRotationAxis);
            Quaternion_t halfRotationAboutAxis(
                    cos(0.25 * bendAngle), sin(0.25 * bendAngle) * effectiveRotationAxis);

            // Exit point = entrance + chord. The chord has the chord length and points half the
            // bend angle off the entrance tangent.
            Vector_t<double, 3> chord =
                    thisLength * halfRotationAboutAxis.rotate(Vector_t<double, 3>(0, 0, 1));
            Vector_t<double, 3> endThis3D = beginThis3D + chord;
            // Advance the running path length by the ARC (not the chord), so a following ELEMEDGE
            // element sits flush against the bend's arc exit.
            double endThisPathLength = beginThisPathLength + arcLength;

            // This element's own frame is its GEOMETRICAL ENTRANCE frame (lab -> entrance face),
            // the same frame the field kernel evaluates in:
            //  - SBEND: a curved body, so its entrance face frame coincides with the design-orbit
            //    entrance tangent (identity rotation relative to the incoming orbit here).
            //    SBend::apply measures arc length from it directly, no pole-face de-tilt.
            //  - RBEND/RBEND3D: a straight box, so its entrance face frame is the box frame — +z
            //    along the chord, half the bend angle off the incoming orbit tangent. It does NOT
            //    coincide with the orbit tangent. RBend::apply gates on the box z, uniform dipole.
            // E1/E2 pole-face rotations are rejected at parse time, so there is no face tilt to add.
            // Only this element's own frame is set here; the running frame handed to the next
            // element (currentCoordTrafo) is advanced separately below.
            const bool isRectangular = element->getType() == ElementType::RBEND
                                       || element->getType() == ElementType::RBEND3D;
            const Quaternion_t entryFrameRotation =
                    isRectangular ? halfRotationAboutAxis
                                  : Quaternion_t(1.0, Vector_t<double, 3>(0, 0, 0));

            // Map prev-exit -> this entrance frame, and prev-exit -> this exit frame.
            CoordinateSystemTrafo fromEndLastToBeginThis(
                    beginThis3D, (entryFrameRotation * rotationAboutZ).conjugate());
            CoordinateSystemTrafo fromEndLastToEndThis(endThis3D, rotationAboutAxis.conjugate());

            // Global->local for this bend = (running lab frame) then (prev exit -> this entrance).
            element->setCSTrafoGlobal2Local(fromEndLastToBeginThis * currentCoordTrafo);

            // Advance the running frame and path length to this bend's exit.
            currentCoordTrafo = (fromEndLastToEndThis * currentCoordTrafo);

            endPriorPathLength = endThisPathLength;
        }
    }

    // Phase 3 — place every remaining (non-bend) element and fix all positions. Re-walk from the
    // lab frame; the running frame again only turns at bends (their own frames were set in Phase 2).
    double endPriorPathLength               = 0.0;
    CoordinateSystemTrafo currentCoordTrafo = labFrame;

    ElementList::iterator it = elements.begin();
    for (; it != end; ++it) {
        std::shared_ptr<ElementBase> element = (*it);
        if (element->isPositioned()) continue;  // already placed by Phase 1 (Mode A)

        // Entrance offset from the previous element along the straight run.
        double beginThisPathLength = element->getElementPosition();
        Vector_t<double, 3> beginThis3D(0, 0, beginThisPathLength - endPriorPathLength);

        Vector_t<double, 3> endThis3D;
        if (element->getType() == ElementType::SBEND || element->getType() == ElementType::RBEND
            || element->getType() == ElementType::RBEND3D) {
            // Bend: its own frame is already set (Phase 2); here we only advance the running frame
            // across it so the following elements are placed correctly.
            double thisLength = element->getGeometry().getChordLength();
            double bendAngle  = element->getGeometry().getBendAngle();

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

            const double arcLength = element->getGeometry().getArcLength();

            // Exit point via the chord, then advance the running frame and path length (by arc).
            endThis3D =
                    (beginThis3D
                     + halfRotationAboutAxis.rotate(Vector_t<double, 3>(0, 0, thisLength)));
            CoordinateSystemTrafo fromEndLastToEndThis(endThis3D, rotationAboutAxis.conjugate());
            currentCoordTrafo = fromEndLastToEndThis * currentCoordTrafo;

            endPriorPathLength = beginThisPathLength + arcLength;
        } else {
            // Straight element: entrance offset plus roll about z, composed with the running frame.
            double rotationAngleAboutZ = (*it)->getRotationAboutZ();
            Quaternion_t rotationAboutZ(
                    cos(0.5 * rotationAngleAboutZ),
                    sin(-0.5 * rotationAngleAboutZ) * Vector_t<double, 3>(0, 0, 1));

            CoordinateSystemTrafo fromLastToThis(beginThis3D, rotationAboutZ);

            element->setCSTrafoGlobal2Local(fromLastToThis * currentCoordTrafo);
        }

        element->fixPosition();
    }
}
