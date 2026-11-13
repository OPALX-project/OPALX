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
#include <memory>
#include <string>
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/SBendRep.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"
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

std::shared_ptr<ElementBase> OpalSBend::build(const ElementValues& values) {
    auto element = std::make_shared<SBendRep>(values.getName());
    fill(*element, values);
    return element;
}

void OpalSBend::update() { fill(*getElement(), getValues()); }

void OpalSBend::fill(ElementBase& element, const ElementValues& values) {
    setCommonFields(element, values);

    // E1/E2 pole-face rotations are not wired into the OPALX-native bend geometry/field yet, so
    // reject them explicitly rather than silently ignoring them.
    if (!values.defaultUsed("E1") || !values.defaultUsed("E2")) {
        throw OpalException(
                "OpalSBend::update",
                values.getName()
                        + ": pole-face rotations (E1, E2) are not yet implemented in the "
                          "OPALX-native bend port. Remove E1/E2 from the element definition.");
    }

    // Define geometry.
    SBendRep* bend     = dynamic_cast<SBendRep*>(&element);
    double length      = values.getReal("L");
    double angle       = values.getReal("ANGLE");
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
    double k0     = values.hasValue("K0") ? values.getReal("K0") : length ? angle / length : angle;
    const std::vector<double> normal = {
            factor * k0, factor * values.getReal("K1"), factor * values.getReal("K2") / 2.0,
            factor * values.getReal("K3") / 6.0};
    const std::vector<double> skew = {
            factor * values.getReal("K0S"), factor * values.getReal("K1S"),
            factor * values.getReal("K2S") / 2.0, factor * values.getReal("K3S") / 6.0};
    bend->setFieldComponents(normal, skew);

    // Set field amplitude or bend angle.
    if (values.hasValue("ANGLE")) {
        if (bend->isPositioned() && angle < 0.0) {
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

    if (values.hasValue("GREATERTHANPI"))
        throw OpalException("OpalSBend::update", "GREATERTHANPI not supportet any more");

    // Units are eV.
    if (values.hasValue("DESIGNENERGY")) {
        // bend->setDesignEnergy(values.getReal("DESIGNENERGY") * Units::MeV2eV, false);
        throw OpalException(
                "OpalSBend::update",
                "DESIGNENERGY is not supported yet for the OPALX-native SBEND port.");
    }

    if (values.hasValue("GAP"))
        throw OpalException(
                "OpalSBend::update", "GAP is not supported; specify HGAP (the half gap) instead.");
    // gap_m stores the full gap; HGAP is the half gap, so double it.
    bend->setFullGap(2.0 * values.getReal("HGAP"));
    bend->setFringeIntegral(values.getReal("FINT"));

    // Transverse aperture. The vertical half-aperture is the pole gap HGAP; the horizontal
    // one comes from HAPERT or from APERTURE, never from both. setCommonFields() has
    // already installed either the APERTURE string (full widths, halved by getApert()) or
    // the 1e6 default, so this block only adds the gap and rejects the combinations that
    // have no physical meaning.
    //
    // "Given" means written in the input. HAPERT has a default, so values.hasValue("HAPERT") is
    // true even when the input omits it; only defaultUsed() tells the two apart.
    const bool hasApert  = !values.defaultUsed("APERTURE");
    const bool hasHapert = !values.defaultUsed("HAPERT");
    const double hgap    = values.getReal("HGAP");
    const double hapert  = values.getReal("HAPERT");

    if (hasApert && hasHapert) {
        throw OpalException(
                "OpalSBend::update",
                values.getName()
                        + ": APERTURE and HAPERT cannot both be given. APERTURE sets both "
                          "half-apertures, HAPERT only the horizontal one.");
    }
    if ((hasApert || hasHapert) && hgap <= 0.0) {
        throw OpalException(
                "OpalSBend::update",
                values.getName()
                        + ": HAPERT/APERTURE requires a positive HGAP; the vertical aperture "
                          "is the half gap.");
    }
    if (hasHapert && hapert <= 0.0) {
        throw OpalException("OpalSBend::update", values.getName() + ": HAPERT must be positive.");
    }

    auto aperture = bend->getAperture();
    if (hasApert) {
        // Keep the input's own shape and half-widths, but the chamber cannot be taller than
        // the poles it sits between. Equality is fine: flush with the pole faces.
        if (aperture.second[1] > hgap) {
            throw OpalException(
                    "OpalSBend::update",
                    values.getName() + ": vertical half-aperture "
                            + std::to_string(aperture.second[1])
                            + " m exceeds the half gap HGAP = " + std::to_string(hgap) + " m.");
        }
    } else if (hgap > 0.0) {
        // Flat poles above and below, a flat side wall left and right. Without HAPERT the
        // horizontal half-width stays at the 1e6 default, i.e. vertical scraping only.
        aperture.first     = ApertureType::RECTANGULAR;
        aperture.second[1] = hgap;
        if (hasHapert) {
            aperture.second[0] = hapert;
        }
        bend->setAperture(aperture.first, aperture.second);
    }

    if (values.hasValue("L"))
        bend->getGeometry().setElementLength(values.getReal("L"));
    else
        bend->getGeometry().setElementLength(0.0);

    if (values.hasValue("WAKEF")) {
        throw OpalException(
                "OpalSBend::update", "WAKEF is not supported yet for the OPALX-native SBEND port.");
    }

    if (values.hasValue("PARTICLEMATTERINTERACTION")) {
        throw OpalException(
                "OpalSBend::update",
                "PARTICLEMATTERINTERACTION is not supported yet for the OPALX-native SBEND port.");
    }
}
