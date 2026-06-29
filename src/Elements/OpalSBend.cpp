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
#include "Attributes/Attributes.h"
#include "BeamlineCore/SBendRep.h"
#include "Fields/BMultipoleField.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"

namespace {

    /**
     * @brief Return the OPAL-compatible SBEND dipole normalization length.
     *
     * Historical OPAL interprets the `SBEND` `L` attribute as the magnet chord
     * length when it derives the design radius:
     * \f[
     *   R = \frac{L}{2\sin(\theta/2)} .
     * \f]
     * The dipole coefficient is therefore normalized with the corresponding arc
     * length
     * \f[
     *   L_\mathrm{norm} = R\theta
     *                  = L\frac{\theta}{2\sin(\theta/2)},
     *   \qquad k_0 = \frac{\theta}{L_\mathrm{norm}} .
     * \f]
     * The small-angle limit returns `L`, preserving the straight-element limit.
     */
    double getOpalChordToArcNormalizationLength(const double chordLength, const double angle) {
        const double halfAngle = 0.5 * angle;
        if (std::abs(halfAngle) <= 1.0e-12) {
            return chordLength;
        }

        const double sine = std::sin(halfAngle);
        if (std::abs(sine) <= 1.0e-15) {
            return chordLength;
        }

        return chordLength * halfAngle / sine;
    }

}  // namespace

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

    // OpalData updates all registered objects, including the exemplar
    // "SBEND" definition. The exemplar carries parser metadata only and
    // must not be forced to satisfy analytic bend input constraints.
    if (getParent() == nullptr && getOpalName() == "SBEND") {
        return;
    }

    SBendRep* bend       = dynamic_cast<SBendRep*>(getElement());
    double length        = Attributes::getReal(itsAttr[LENGTH]);
    double angle         = Attributes::getReal(itsAttr[ANGLE]);
    double e1            = Attributes::getReal(itsAttr[E1]);
    double e2            = Attributes::getReal(itsAttr[E2]);
    const double fullGap = Attributes::getReal(itsAttr[GAP]);
    const double halfGap =
            itsAttr[HGAP] ? Attributes::getReal(itsAttr[HGAP]) : 0.5 * std::abs(fullGap);
    const double fringeIntegral           = Attributes::getReal(itsAttr[FINT]);
    const double fieldNormalizationLength = getOpalChordToArcNormalizationLength(length, angle);
    PlanarArcGeometry& geometry           = bend->getGeometry();

    if (length) {
        geometry = PlanarArcGeometry(fieldNormalizationLength, angle / fieldNormalizationLength);
    } else {
        geometry = PlanarArcGeometry(angle);
    }

    bend->setBendAngle(angle);
    bend->setEntranceAngle(e1);
    bend->setExitAngle(e2);
    bend->setFullGap(fullGap);
    bend->setFringeHalfGap(halfGap);
    bend->setFringeIntegral(fringeIntegral);

    validateAnalyticBendDefinition(
            "OpalSBend::update(" + getOpalName() + ")", !itsAttr[ANGLE].defaultUsed(),
            !itsAttr[K0].defaultUsed(), fieldNormalizationLength, angle,
            Attributes::getReal(itsAttr[K0]));

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
    double k0  = deriveAnalyticDipoleCoefficient(fieldNormalizationLength, angle);
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

    bend->setBendAngle(angle);
    bend->setFieldAmplitude(k0, k0s);

    if (itsAttr[GREATERTHANPI])
        throw OpalException("OpalSBend::update", "GREATERTHANPI not supportet any more");

    if (itsAttr[ROTATION])
        throw OpalException(
                "OpalSBend::update", "ROTATION not supportet any more; use PSI instead");

    if (itsAttr[FMAPFN])
        bend->setFieldMapFN(Attributes::getString(itsAttr[FMAPFN]));
    else if (bend->getName() != "SBEND") {
        bend->setFieldMapFN("1DPROFILE1-DEFAULT");
    }

    bend->setEntranceAngle(e1);
    bend->setExitAngle(e2);

    // Units are eV.
    if (itsAttr[DESIGNENERGY]) {
        bend->setDesignEnergy(Attributes::getReal(itsAttr[DESIGNENERGY]), false);
    }

    if (itsAttr[APERT])
        throw OpalException(
                "OpalSBend::update", "APERTURE in SBEND not supported; use GAP and HAPERT instead");

    if (itsAttr[HAPERT]) {
        double hapert = Attributes::getReal(itsAttr[HAPERT]);
        if (hapert > 0.0) {
            bend->setAperture(
                    ApertureType::RECTANGULAR, std::vector<double>({hapert, hapert, 1.0}));
        }
    }

    if (itsAttr[WAKEF]) {
        throw OpalException(
                "OpalSBend::update", "WAKEF is not supported yet for the OPALX-native SBEND port.");
    }

    if (itsAttr[K1])
        bend->setK1(Attributes::getReal(itsAttr[K1]));
    else
        bend->setK1(0.0);

    if (itsAttr[PARTICLEMATTERINTERACTION]) {
        throw OpalException(
                "OpalSBend::update",
                "PARTICLEMATTERINTERACTION is not supported yet for the OPALX-native SBEND port.");
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(bend);
}
