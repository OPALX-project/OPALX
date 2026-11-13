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
#include <memory>
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

std::shared_ptr<ElementBase> OpalRBend::build(const ElementValues& values) {
    auto element = std::make_shared<RBendRep>(values.getName());
    fill(*element, values);
    return element;
}

void OpalRBend::update() { fill(*getElement(), getValues()); }

void OpalRBend::fill(ElementBase& element, const ElementValues& values) {
    setCommonFields(element, values);

    // Geometry, field profile and calibration are independent RBEND settings.
    const bool hasProfile = !values.defaultUsed("FMAPFN");
    const bool hasE1      = !values.defaultUsed("E1");
    const bool hasE2      = !values.defaultUsed("E2");
    RBendRep* bend        = dynamic_cast<RBendRep*>(&element);
    const double length   = values.getReal("L");
    double angle          = values.getReal("ANGLE");
    Geometry& geometry    = static_cast<Geometry&>(bend->getGeometry());
    geometry.setElementLength(length);
    geometry.setBendAngle(angle);
    geometry.setPoleFaceAnglesEnabled(hasE1 || hasE2);
    const double e1 = hasE1   ? values.getReal("E1")
                      : hasE2 ? angle - values.getReal("E2")
                              : 0.5 * angle;
    const double e2 = hasE2 ? values.getReal("E2") : angle - e1;
    if (hasE1 || hasE2 || hasProfile) {
        if (!std::isfinite(length) || length <= 0.0 || !std::isfinite(angle)
            || !values.hasValue("ANGLE") || std::abs(angle) >= Physics::pi || !std::isfinite(e1)
            || !std::isfinite(e2) || std::abs(e1) >= 0.5 * Physics::pi
            || std::abs(e2) >= 0.5 * Physics::pi)
            throw OpalException("OpalRBend::update", values.getName()
                + ": RBEND requires L > 0, ANGLE, and pole-face angles strictly between -pi/2 and pi/2.");
        if (std::abs(e1 + e2 - angle) > 1e-12)
            throw OpalException("OpalRBend::update", values.getName()
                + ": parallel rectangular pole faces require E1 + E2 = ANGLE; omit either angle to derive it.");
    }
    geometry.setEntranceAngle(e1);
    geometry.setExitAngle(e2);
    if (hasProfile) {
        for (const char* attr : {"K0", "K0S", "K1", "K1S", "K2", "K2S", "K3", "K3S", "FINT"})
            if (!values.defaultUsed(attr))
                throw OpalException("OpalRBend::update", values.getName()
                    + ": a 1DProfile1 field is an ANGLE-calibrated dipole; explicit multipoles and FINT are unsupported.");
    }
    if (values.hasValue("H1") || values.hasValue("H2"))
        throw OpalException("OpalRBend::update", "Curved pole faces H1/H2 are not supported.");

    // Field coefficients use the configured reference-orbit curvature. DESIGNENERGY
    // changes their rigidity at initialization, once the particle species is known.
    double factor                    = OpalData::getInstance()->getP0() / Physics::c;
    double k0                        = values.hasValue("K0") ? values.getReal("K0")
                                       : length ? (geometry.hasPoleFaceAngles() ? angle / geometry.getArcLength()
                                                                                : 2 * sin(angle / 2) / length)
                                                : angle;
    const std::vector<double> normal = {
            factor * k0, factor * values.getReal("K1"), factor * values.getReal("K2") / 2.0,
            factor * values.getReal("K3") / 6.0};
    const std::vector<double> skew = {
            factor * values.getReal("K0S"), factor * values.getReal("K1S"),
            factor * values.getReal("K2S") / 2.0, factor * values.getReal("K3S") / 6.0};
    bend->setFieldComponents(normal, skew);
    bend->setReferenceMomentum(OpalData::getInstance()->getP0());

    // Set the bend angle.
    if (values.hasValue("ANGLE")) {
        if (!geometry.hasPoleFaceAngles() && bend->isPositioned() && angle < 0.0) {
            angle = -angle;

            Quaternion rotAboutZ(0, 0, 0, 1);
            CoordinateSystemTrafo g2l = bend->getCSTrafoGlobal2Local();
            bend->releasePosition();
            bend->setCSTrafoGlobal2Local(
                    CoordinateSystemTrafo(g2l.getOrigin(), rotAboutZ * g2l.getRotation()));
            bend->fixPosition();
        }
        geometry.setBendAngle(angle);
        if (!geometry.hasPoleFaceAngles()) {
            geometry.setEntranceAngle(0.5 * angle);
            geometry.setExitAngle(0.5 * angle);
        }
    }

    // Energy in eV.
    if (values.hasValue("DESIGNENERGY")) {
        const double energy = values.getReal("DESIGNENERGY");
        if (!std::isfinite(energy) || energy <= 0.0)
            throw OpalException("OpalRBend::update", "DESIGNENERGY must be positive and finite.");
        bend->setDesignEnergy(energy * Units::MeV2eV, false);
    }

    const bool hasGap  = !values.defaultUsed("GAP");
    const bool hasHgap = !values.defaultUsed("HGAP");
    if (hasGap && hasHgap)
        throw OpalException("OpalRBend::update", "Specify GAP or HGAP, not both.");
    double fullGap = hasGap ? values.getReal("GAP") : 2.0 * values.getReal("HGAP");
    if (!std::isfinite(fullGap) || fullGap < 0.0)
        throw OpalException("OpalRBend::update", "The pole gap must be finite and nonnegative.");
    // gap_m stores the full gap; HGAP is the half gap, so double it.
    bend->setFullGap(fullGap);
    bend->setFringeIntegral(values.getReal("FINT"));
    const std::string filename = hasProfile ? values.getString("FMAPFN") : "";
    if (hasProfile && filename.empty())
        throw OpalException("OpalRBend::update", "FMAPFN cannot be empty.");
    bend->setFieldMap(filename);
    if (hasProfile && fullGap == 0.0) fullGap = bend->getFullGap();

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
    const double hgap    = 0.5 * fullGap;
    const double hapert  = values.getReal("HAPERT");

    if (hasApert && hasHapert) {
        throw OpalException(
                "OpalRBend::update",
                values.getName()
                        + ": APERTURE and HAPERT cannot both be given. APERTURE sets both "
                          "half-apertures, HAPERT only the horizontal one.");
    }
    if ((hasApert || hasHapert) && hgap <= 0.0) {
        throw OpalException(
                "OpalRBend::update",
                values.getName()
                        + ": HAPERT/APERTURE requires a positive HGAP; the vertical aperture "
                          "is the half gap.");
    }
    if (hasHapert && hapert <= 0.0) {
        throw OpalException("OpalRBend::update", values.getName() + ": HAPERT must be positive.");
    }

    auto aperture = bend->getAperture();
    if (hasApert) {
        // Keep the input's own shape and half-widths, but the chamber cannot be taller than
        // the poles it sits between. Equality is fine: flush with the pole faces.
        if (aperture.second[1] > hgap) {
            throw OpalException(
                    "OpalRBend::update",
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
                "OpalRBend::update", "WAKEF is not supported yet for the OPALX-native RBEND port.");
    }

    if (values.hasValue("PARTICLEMATTERINTERACTION")) {
        throw OpalException(
                "OpalRBend::update",
                "PARTICLEMATTERINTERACTION is not supported yet for the OPALX-native RBEND port.");
    }
}
