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
#include <string>
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/RBendRep.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"
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

    // Explicit pole-face angles, design energy or a profile select OPAL-T compatibility. An
    // ordinary native RBEND keeps its existing box field and FINT correction.
    const bool legacy = itsAttr[FMAPFN] || itsAttr[DESIGNENERGY] || !itsAttr[E1].defaultUsed()
                        || !itsAttr[E2].defaultUsed();

    // Define geometry. L is the straight body (box) length — the placed hardware, and the
    // chord of the design orbit. The orbit arc through the box is longer,
    // arc = L (angle/2) / sin(angle/2), reported by Geometry::getArcLength().
    RBendRep* bend     = dynamic_cast<RBendRep*>(getElement());
    double length      = Attributes::getReal(itsAttr[LENGTH]);
    double angle       = Attributes::getReal(itsAttr[ANGLE]);
    Geometry& geometry = static_cast<Geometry&>(bend->getGeometry());
    geometry.setElementLength(length);
    geometry.setBendAngle(angle);
    geometry.setLegacyRBend(legacy);
    const double e1 = Attributes::getReal(itsAttr[E1]);
    if (legacy) {
        if (!std::isfinite(length) || length <= 0.0 || !std::isfinite(angle) || !itsAttr[ANGLE]
            || std::abs(angle) >= Physics::pi || !std::isfinite(e1)
            || std::abs(e1) >= 0.5 * Physics::pi || std::abs(angle - e1) >= 0.5 * Physics::pi)
            throw OpalException("OpalRBend::update", getOpalName()
                + ": legacy RBEND requires L > 0, ANGLE, and entrance/exit angles strictly between -pi/2 and pi/2.");
        if (!itsAttr[E2].defaultUsed()
            && (!std::isfinite(Attributes::getReal(itsAttr[E2]))
                || std::abs(Attributes::getReal(itsAttr[E2]) - (angle - e1)) > 1e-12))
            throw OpalException("OpalRBend::update", getOpalName()
                + ": rectangular pole faces require E2 = ANGLE - E1; omit E2 to derive it.");
        for (int attr : {K0, K0S, K1, K1S, K2, K2S, K3, K3S, FINT})
            if (!itsAttr[attr].defaultUsed())
                throw OpalException("OpalRBend::update", getOpalName()
                    + ": the legacy profile supports an ANGLE-calibrated dipole only; multipoles, K0/K0S and explicit FINT cannot be combined with it.");
        geometry.setEntranceAngle(e1);
        geometry.setExitAngle(angle - e1);
    }
    if (itsAttr[H1] || itsAttr[H2])
        throw OpalException("OpalRBend::update", "Curved pole faces H1/H2 are not supported.");

    // Define field. With L the box (chord) length, the design radius is
    // rho = L / (2 sin(angle/2)), so the default dipole strength is
    // k0 = 1 / rho = 2 sin(angle/2) / L (= angle / arc length).
    double factor                    = OpalData::getInstance()->getP0() / Physics::c;
    double k0                        = itsAttr[K0] ? Attributes::getReal(itsAttr[K0])
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

    // Set the bend angle.
    if (itsAttr[ANGLE]) {
        if (!legacy && bend->isPositioned() && angle < 0.0) {
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

    // Energy in eV.
    if (itsAttr[DESIGNENERGY]) {
        const double energy = Attributes::getReal(itsAttr[DESIGNENERGY]);
        if (!std::isfinite(energy) || energy <= 0.0)
            throw OpalException("OpalRBend::update", "DESIGNENERGY must be positive and finite.");
        bend->setDesignEnergy(energy * Units::MeV2eV, false);
    } else if (legacy) {
        throw OpalException("OpalRBend::update", "Legacy RBEND requires DESIGNENERGY in MeV.");
    }

    const bool hasGap  = !itsAttr[GAP].defaultUsed();
    const bool hasHgap = !itsAttr[HGAP].defaultUsed();
    if (hasGap && hasHgap)
        throw OpalException("OpalRBend::update", "Specify GAP or HGAP, not both.");
    double fullGap =
            hasGap ? Attributes::getReal(itsAttr[GAP]) : 2.0 * Attributes::getReal(itsAttr[HGAP]);
    if (!std::isfinite(fullGap) || fullGap < 0.0)
        throw OpalException("OpalRBend::update", "The pole gap must be finite and nonnegative.");
    // gap_m stores the full gap; HGAP is the half gap, so double it.
    bend->setFullGap(fullGap);
    bend->setFringeIntegral(Attributes::getReal(itsAttr[FINT]));
    bend->setLegacyFieldMap(
            legacy ? (itsAttr[FMAPFN] ? Attributes::getString(itsAttr[FMAPFN])
                                      : "1DPROFILE1-DEFAULT")
                   : "");
    if (legacy && !bend->usesLegacyFieldMap())
        throw OpalException("OpalRBend::update", "Legacy FMAPFN cannot be empty.");
    if (legacy && fullGap == 0.0) fullGap = bend->getFullGap();

    // Transverse aperture. The vertical half-aperture is the pole gap HGAP; the horizontal
    // one comes from HAPERT or from APERTURE, never from both. OpalElement::update() has
    // already installed either the APERTURE string (full widths, halved by getApert()) or
    // the 1e6 default, so this block only adds the gap and rejects the combinations that
    // have no physical meaning.
    //
    // "Given" means written in the deck. HAPERT has a default, so itsAttr[HAPERT] reads as
    // set even when the deck omits it; only defaultUsed() tells the two apart.
    const bool hasApert  = !itsAttr[APERT].defaultUsed();
    const bool hasHapert = !itsAttr[HAPERT].defaultUsed();
    const double hgap    = 0.5 * fullGap;
    const double hapert  = Attributes::getReal(itsAttr[HAPERT]);

    if (hasApert && hasHapert) {
        throw OpalException(
                "OpalRBend::update",
                getOpalName()
                        + ": APERTURE and HAPERT cannot both be given. APERTURE sets both "
                          "half-apertures, HAPERT only the horizontal one.");
    }
    if ((hasApert || hasHapert) && hgap <= 0.0) {
        throw OpalException(
                "OpalRBend::update",
                getOpalName()
                        + ": HAPERT/APERTURE requires a positive HGAP; the vertical aperture "
                          "is the half gap.");
    }
    if (hasHapert && hapert <= 0.0) {
        throw OpalException("OpalRBend::update", getOpalName() + ": HAPERT must be positive.");
    }

    auto aperture = bend->getAperture();
    if (hasApert) {
        // Keep the deck's own shape and half-widths, but the chamber cannot be taller than
        // the poles it sits between. Equality is fine: flush with the pole faces.
        if (aperture.second[1] > hgap) {
            throw OpalException(
                    "OpalRBend::update",
                    getOpalName() + ": vertical half-aperture " + std::to_string(aperture.second[1])
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
