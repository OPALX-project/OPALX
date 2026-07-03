//
// Class OpalRBend
//   The RBEND element.
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
#include "Elements/OpalRBend.h"
#include <cmath>
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/RBendRep.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"

OpalRBend::OpalRBend()
    : OpalBend("RBEND", "The \"RBEND\" element defines a rectangular bending magnet.") {
    registerOwnership();

    setElement(new RBendRep("RBEND"));
}

OpalRBend::OpalRBend(const std::string& name, OpalRBend* parent) : OpalBend(name, parent) {
    setElement(new RBendRep(name));
}

OpalRBend::~OpalRBend() {}

OpalRBend* OpalRBend::clone(const std::string& name) { return new OpalRBend(name, this); }

void OpalRBend::update() {
    OpalElement::update();

    // Define geometry. L is the straight body (box) length — the placed hardware, and the
    // chord of the design orbit. The orbit arc through the box is longer,
    // arc = L (angle/2) / sin(angle/2), reported by Geometry::getArcLength().
    RBendRep* bend     = dynamic_cast<RBendRep*>(getElement());
    double length      = Attributes::getReal(itsAttr[LENGTH]);
    double angle       = Attributes::getReal(itsAttr[ANGLE]);
    double e1          = Attributes::getReal(itsAttr[E1]);
    double e2          = Attributes::getReal(itsAttr[E2]);
    Geometry& geometry = static_cast<Geometry&>(bend->getGeometry());
    geometry.setElementLength(length);
    if (angle < 0) {
    }
    geometry.setBendAngle(angle);

    // Define field. With L the box (chord) length, the design radius is
    // rho = L / (2 sin(angle/2)), so the default dipole strength is
    // k0 = 1 / rho = 2 sin(angle/2) / L (= angle / arc length).
    double factor = OpalData::getInstance()->getP0() / Physics::c;
    double k0     = itsAttr[K0] ? Attributes::getReal(itsAttr[K0])
                    : length    ? 2 * sin(angle / 2) / length
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

    // Set the bend angle and pole-face rotations.
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
    geometry.setEntranceAngle(e1);
    geometry.setExitAngle(e2);

    // Energy in eV.
    if (itsAttr[DESIGNENERGY]) {
        bend->setDesignEnergy(Attributes::getReal(itsAttr[DESIGNENERGY]), false);
    }

    if (itsAttr[HGAP])
        throw OpalException(
                "OpalRBend::update", "HGAP is not supported; specify the full GAP instead.");
    bend->setFullGap(Attributes::getReal(itsAttr[GAP]));
    bend->setFringeIntegral(Attributes::getReal(itsAttr[FINT]));

    if (itsAttr[HAPERT]) {
        double hapert = Attributes::getReal(itsAttr[HAPERT]);
        bend->setAperture(ApertureType::RECTANGULAR, std::vector<double>({hapert, hapert, 1.0}));
    }

    if (itsAttr[LENGTH])
        bend->getGeometry().setElementLength(Attributes::getReal(itsAttr[LENGTH]));
    else
        bend->getGeometry().setElementLength(0.0);

    if (itsAttr[WAKEF]) {
        throw OpalException(
                "OpalRBend::update", "WAKEF is not supported yet for the OPALX-native RBEND port.");
    }

    if (itsAttr[PARTICLEMATTERINTERACTION]) {
        throw OpalException(
                "OpalRBend::update",
                "PARTICLEMATTERINTERACTION is not supported yet for the OPALX-native RBEND port.");
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(bend);
}
