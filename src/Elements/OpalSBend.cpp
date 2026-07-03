//
// Class OpalSBend
//   The SBEND element.
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "Elements/OpalSBend.h"
#include <cmath>
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/SBendRep.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"

OpalSBend::OpalSBend()
    : OpalBend("SBEND", "The \"SBEND\" element defines a sector bending magnet.") {
    registerOwnership();

    setElement((new SBendRep("SBEND")));
}

OpalSBend::OpalSBend(const std::string& name, OpalSBend* parent) : OpalBend(name, parent) {
    setElement((new SBendRep(name)));
}

OpalSBend::~OpalSBend() {}

OpalSBend* OpalSBend::clone(const std::string& name) { return new OpalSBend(name, this); }

void OpalSBend::update() {
    OpalElement::update();

    // Define geometry.
    SBendRep* bend     = dynamic_cast<SBendRep*>(getElement());
    double length      = Attributes::getReal(itsAttr[LENGTH]);
    double angle       = Attributes::getReal(itsAttr[ANGLE]);
    double e1          = Attributes::getReal(itsAttr[E1]);
    double e2          = Attributes::getReal(itsAttr[E2]);
    Geometry& geometry = static_cast<Geometry&>(bend->getGeometry());

    if (length) {
        geometry = Geometry::makeSBend(length, angle / length);
    } else {
        geometry = Geometry::makeSBend(0.0, 0.0);
        geometry.setBendAngle(angle);
    }
    // Define field. MAD length convention: L is the design-orbit arc length, so the
    // default dipole strength is the arc curvature k0 = angle / L = 1 / rho. This matches
    // the geometry above (makeSBend treats L as the arc), so the field bends the orbit on
    // exactly the placed arc at any angle. (The legacy chord form 2 sin(angle/2) / L
    // belonged to the OPAL convention where L was the chord.)
    double factor = OpalData::getInstance()->getP0() / Physics::c;
    double k0     = itsAttr[K0] ? Attributes::getReal(itsAttr[K0])
                    : length    ? angle / length
                                : angle;
    const std::vector<double> normal = {
            factor * k0, factor * Attributes::getReal(itsAttr[K1]),
            factor * Attributes::getReal(itsAttr[K2]) / 2.0,
            factor * Attributes::getReal(itsAttr[K3]) / 6.0};
    const std::vector<double> skew = {
            factor * Attributes::getReal(itsAttr[K0S]), factor * Attributes::getReal(itsAttr[K1S]),
            factor * Attributes::getReal(itsAttr[K2S]) / 2.0,
            factor * Attributes::getReal(itsAttr[K3S]) / 6.0};
    bend->setFieldComponents(normal, skew);

    // Set field amplitude or bend angle.
    if (itsAttr[ANGLE]) {
        if (bend->isPositioned() && angle < 0.0) {
            e1    = -e1;
            e2    = -e2;
            angle = -angle;

            Quaternion rotAboutZ(0, 0, 0, 1);
            CoordinateSystemTrafo g2l = bend->getCSTrafoGlobal2Local();
            bend->releasePosition();
            bend->setCSTrafoGlobal2Local(
                    CoordinateSystemTrafo(g2l.getOrigin(), rotAboutZ * g2l.getRotation()));
            bend->fixPosition();
        }
        geometry.setBendAngle(angle);
    }

    if (itsAttr[GREATERTHANPI])
        throw OpalException("OpalSBend::update", "GREATERTHANPI not supportet any more");

    geometry.setEntranceAngle(e1);
    geometry.setExitAngle(e2);

    // Units are eV.
    if (itsAttr[DESIGNENERGY]) {
        bend->setDesignEnergy(Attributes::getReal(itsAttr[DESIGNENERGY]), false);
    }

    if (itsAttr[HGAP])
        throw OpalException(
                "OpalSBend::update", "HGAP is not supported; specify the full GAP instead.");
    bend->setFullGap(Attributes::getReal(itsAttr[GAP]));
    bend->setFringeIntegral(Attributes::getReal(itsAttr[FINT]));

    if (itsAttr[APERT])
        throw OpalException(
                "OpalSBend::update", "APERTURE in SBEND not supported; use GAP and HAPERT instead");

    // Only override the aperture when a positive HAPERT is actually given. HAPERT has a
    // default of 0.0, and itsAttr[HAPERT] reads as "set" even when the deck omits it, so
    // gating on presence would install a zero-width rectangular aperture that makes
    // isInsideTransverse always false and silently drops the bend from the beamline.
    if (Attributes::getReal(itsAttr[HAPERT]) > 0.0) {
        double hapert = Attributes::getReal(itsAttr[HAPERT]);
        bend->setAperture(ApertureType::RECTANGULAR, std::vector<double>({hapert, hapert, 1.0}));
    }

    if (itsAttr[LENGTH])
        bend->getGeometry().setElementLength(Attributes::getReal(itsAttr[LENGTH]));
    else
        bend->getGeometry().setElementLength(0.0);

    if (itsAttr[WAKEF]) {
        throw OpalException(
                "OpalSBend::update", "WAKEF is not supported yet for the OPALX-native SBEND port.");
    }

    if (itsAttr[PARTICLEMATTERINTERACTION]) {
        throw OpalException(
                "OpalSBend::update",
                "PARTICLEMATTERINTERACTION is not supported yet for the OPALX-native SBEND port.");
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(bend);
}
