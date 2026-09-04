//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#include "Elements/OpalMultipoleT.h"

#include "Attributes/Attributes.h"
#include "BeamlineCore/MultipoleTRep.h"
#include "Utilities/OpalException.h"

#include <string>
#include <vector>

OpalMultipoleT::OpalMultipoleT()
    : OpalElement(
              SIZE, "MULTIPOLET",
              "The \"MULTIPOLET\" element defines a combined function multipole.") {
    itsAttr[TP] = Attributes::makeRealArray(
            "TP",
            "Mid-plane field profile: the coefficients b_k of B = b_k x^k in the flat top "
            "[T m^(-k)], dipole first, sign included");
    itsAttr[LFRINGE] = Attributes::makeReal(
            "LFRINGE", "The length of the left (entrance) end field [m]; 0 is a hard edge");
    itsAttr[RFRINGE] = Attributes::makeReal(
            "RFRINGE", "The length of the right (exit) end field [m]; 0 is a hard edge");
    itsAttr[MAXFORDER] = Attributes::makeReal(
            "MAXFORDER", "Number of terms used in the vertical field expansion",
            DefaultMAXFORDER);
    itsAttr[ANGLE] = Attributes::makeReal(
            "ANGLE", "The bend angle of a curved magnet [rad]; 0 makes it straight", 0.0);
    itsAttr[HAPERT] = Attributes::makeReal("HAPERT", "The full aperture width [m]");
    itsAttr[VAPERT] = Attributes::makeReal("VAPERT", "The full aperture height [m]");
    itsAttr[SCALING_MODEL] = Attributes::makeUpperCaseString(
            "SCALING_MODEL",
            "The name of the time dependence model, which should give a scaling factor.");

    registerOwnership();
    setElement(new MultipoleTRep("MULTIPOLET"));
}

OpalMultipoleT::OpalMultipoleT(const std::string& name, OpalMultipoleT* parent)
    : OpalElement(name, parent) {
    setElement(new MultipoleTRep(name));
}

OpalMultipoleT* OpalMultipoleT::clone(const std::string& name) {
    return new OpalMultipoleT(name, this);
}

void OpalMultipoleT::print(std::ostream& os) const {
    OpalElement::print(os);
}

void OpalMultipoleT::update() {
    OpalElement::update();

    auto* magnet = dynamic_cast<MultipoleTRep*>(getElement());

    // Define geometry. A bend angle makes the body a planar arc of that angle, with L the arc
    // length (the same convention as SBEND); otherwise the body is straight.
    const double length = Attributes::getReal(itsAttr[LENGTH]);
    double angle        = Attributes::getReal(itsAttr[ANGLE]);
    Geometry& geometry  = magnet->getGeometry();

    if (angle == 0.0) {
        geometry = Geometry::makeStraight(length);
    } else {
        if (length != 0.0) {
            geometry = Geometry::makeSBend(length, angle / length);
        } else {
            geometry = Geometry::makeSBend(0.0, 0.0);
            geometry.setBendAngle(angle);
        }

        // A negative bend angle is the positive one in a frame rolled by 180 degrees about z,
        // the same treatment SBEND gives it.
        if (magnet->isPositioned() && angle < 0.0) {
            angle = -angle;

            Quaternion rotAboutZ(0, 0, 0, 1);
            CoordinateSystemTrafo g2l = magnet->getCSTrafoGlobal2Local();
            magnet->releasePosition();
            magnet->setCSTrafoGlobal2Local(
                    CoordinateSystemTrafo(g2l.getOrigin(), rotAboutZ * g2l.getRotation()));
            magnet->fixPosition();

            geometry = (length != 0.0) ? Geometry::makeSBend(length, angle / length)
                                       : Geometry::makeSBend(0.0, 0.0);
            geometry.setBendAngle(angle);
        }
    }

    // Define field. TP holds the physical mid-plane field, so it is handed over unscaled; the
    // setters do the range checks and throw.
    magnet->setMaxFOrder(
            static_cast<unsigned int>(Attributes::getReal(itsAttr[MAXFORDER])));
    magnet->setTransverseProfile(Attributes::getRealArray(itsAttr[TP]));
    magnet->setFringeLengths(
            Attributes::getReal(itsAttr[LFRINGE]), Attributes::getReal(itsAttr[RFRINGE]));

    // Transverse aperture. HAPERT and VAPERT are the full width and height of a rectangular
    // chamber and go together; APERTURE is the generic alternative, already installed by
    // OpalElement::update() (full widths, halved by getApert()).
    //
    // "Given" means written in the deck. HAPERT and VAPERT have defaults, so itsAttr reads as
    // set even when the deck omits them; only defaultUsed() tells the two apart.
    const bool hasApert  = !itsAttr[APERT].defaultUsed();
    const bool hasHapert = !itsAttr[HAPERT].defaultUsed();
    const bool hasVapert = !itsAttr[VAPERT].defaultUsed();
    const double hapert  = Attributes::getReal(itsAttr[HAPERT]);
    const double vapert  = Attributes::getReal(itsAttr[VAPERT]);

    if (hasApert && (hasHapert || hasVapert)) {
        throw OpalException(
                "OpalMultipoleT::update",
                getOpalName()
                        + ": APERTURE and HAPERT/VAPERT cannot both be given. Use APERTURE for "
                          "any shape, or HAPERT and VAPERT for a rectangle.");
    }
    if (hasHapert != hasVapert) {
        throw OpalException(
                "OpalMultipoleT::update",
                getOpalName() + ": HAPERT and VAPERT must be given together.");
    }
    if (hasHapert && (hapert <= 0.0 || vapert <= 0.0)) {
        throw OpalException(
                "OpalMultipoleT::update",
                getOpalName() + ": HAPERT and VAPERT must be positive.");
    }
    if (hasHapert) {
        // The stored aperture holds half widths.
        magnet->setAperture(ApertureType::RECTANGULAR, {0.5 * hapert, 0.5 * vapert});
    }

    magnet->setScalingName(Attributes::getString(itsAttr[SCALING_MODEL]));

    if (itsAttr[WAKEF]) {
        throw OpalException(
                "OpalMultipoleT::update", "WAKEF is not supported yet for MULTIPOLET.");
    }

    if (itsAttr[PARTICLEMATTERINTERACTION]) {
        throw OpalException(
                "OpalMultipoleT::update",
                "PARTICLEMATTERINTERACTION is not supported yet for MULTIPOLET.");
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(magnet);
}
