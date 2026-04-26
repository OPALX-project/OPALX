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
#include "Attributes/Attributes.h"
#include "BeamlineCore/RBendRep.h"
#include "Fields/BMultipoleField.h"
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

    // OpalData updates all registered objects, including the exemplar
    // "RBEND" definition. The exemplar carries parser metadata only and
    // must not be forced to satisfy analytic bend input constraints.
    if (getParent() == nullptr && getOpalName() == "RBEND") {
        return;
    }

    RBendRep* bend       = dynamic_cast<RBendRep*>(getElement());
    double length        = Attributes::getReal(itsAttr[LENGTH]);
    double angle         = Attributes::getReal(itsAttr[ANGLE]);
    double e1            = Attributes::getReal(itsAttr[E1]);
    const double fullGap = Attributes::getReal(itsAttr[GAP]);
    const double halfGap =
            itsAttr[HGAP] ? Attributes::getReal(itsAttr[HGAP]) : 0.5 * std::abs(fullGap);
    const double fringeIntegral = Attributes::getReal(itsAttr[FINT]);
    RBendGeometry& geometry     = bend->getGeometry();
    geometry.setElementLength(length);
    geometry.setBendAngle(angle);

    bend->setLength(length);
    bend->setBendAngle(angle);
    bend->setEntranceAngle(e1);
    bend->setFullGap(fullGap);
    bend->setFringeHalfGap(halfGap);
    bend->setFringeIntegral(fringeIntegral);

    validateAnalyticBendDefinition(
            "OpalRBend::update(" + getOpalName() + ")", !itsAttr[ANGLE].defaultUsed(),
            !itsAttr[K0].defaultUsed(), bend->getEffectiveFieldLength(), angle,
            Attributes::getReal(itsAttr[K0]));

    // Define number of slices for map tracking
    bend->setNSlices(Attributes::getReal(itsAttr[NSLICES]));

    // Define pole face angles.
    bend->setEntryFaceRotation(Attributes::getReal(itsAttr[E1]));
    bend->setExitFaceRotation(Attributes::getReal(itsAttr[E2]));
    bend->setEntryFaceCurvature(Attributes::getReal(itsAttr[H1]));
    bend->setExitFaceCurvature(Attributes::getReal(itsAttr[H2]));

    // Define integration parameters.
    bend->setSlices(Attributes::getReal(itsAttr[SLICES]));
    bend->setStepsize(Attributes::getReal(itsAttr[STEPSIZE]));

    // Define field.
    BMultipoleField field;
    double k0  = deriveAnalyticDipoleCoefficient(bend->getEffectiveFieldLength(), angle);
    double k0s = itsAttr[K0S] ? Attributes::getReal(itsAttr[K0S]) : 0.0;

    field.setNormalComponent(0, k0);
    field.setSkewComponent(0, Attributes::getReal(itsAttr[K0S]));
    field.setNormalComponent(1, Attributes::getReal(itsAttr[K1]));
    field.setSkewComponent(1, Attributes::getReal(itsAttr[K1S]));
    field.setNormalComponent(2, Attributes::getReal(itsAttr[K2]) / 2.0);
    field.setSkewComponent(2, Attributes::getReal(itsAttr[K2S]) / 2.0);
    field.setNormalComponent(3, Attributes::getReal(itsAttr[K3]) / 6.0);
    field.setSkewComponent(3, Attributes::getReal(itsAttr[K3S]) / 6.0);
    bend->setNormalizedField(field);

    if (bend->isPositioned() && angle < 0.0) {
        e1    = -e1;
        angle = -angle;

        Quaternion rotAboutZ(0, 0, 0, 1);
        CoordinateSystemTrafo g2l = bend->getCSTrafoGlobal2Local();
        bend->releasePosition();
        bend->setCSTrafoGlobal2Local(
                CoordinateSystemTrafo(g2l.getOrigin(), rotAboutZ * g2l.getRotation()));
        bend->fixPosition();
    }
    bend->setBendAngle(angle);
    bend->setFieldAmplitude(k0, k0s);
    bend->setEntranceAngle(e1);

    if (itsAttr[ROTATION])
        throw OpalException(
                "OpalRBend::update", "ROTATION not supported any more; use PSI instead");

    if (itsAttr[FMAPFN])
        bend->setFieldMapFN(Attributes::getString(itsAttr[FMAPFN]));
    else if (bend->getName() != "RBEND") {
        bend->setFieldMapFN("1DPROFILE1-DEFAULT");
    }

    // Energy in eV.
    if (itsAttr[DESIGNENERGY]) {
        bend->setDesignEnergy(Attributes::getReal(itsAttr[DESIGNENERGY]), false);
    }

    if (itsAttr[HAPERT]) {
        double hapert = Attributes::getReal(itsAttr[HAPERT]);
        if (hapert > 0.0) {
            bend->setAperture(
                    ApertureType::RECTANGULAR, std::vector<double>({hapert, hapert, 1.0}));
        }
    }

    if (itsAttr[WAKEF]) {
        throw OpalException(
                "OpalRBend::update", "WAKEF is not supported yet for the OPALX-native RBEND port.");
    }

    if (itsAttr[K1])
        bend->setK1(Attributes::getReal(itsAttr[K1]));
    else
        bend->setK1(0.0);

    if (itsAttr[PARTICLEMATTERINTERACTION]) {
        throw OpalException(
                "OpalRBend::update",
                "PARTICLEMATTERINTERACTION is not supported yet for the OPALX-native RBEND port.");
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(bend);
}
