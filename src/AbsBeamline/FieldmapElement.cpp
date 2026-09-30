//
// Class FieldmapElement
//   A beamline element whose only field source is a tabulated field map.
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPALX.
//
// OPALX is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPALX. If not, see <https://www.gnu.org/licenses/>.
//
#include "AbsBeamline/FieldmapElement.h"

#include "AbsBeamline/BeamlineVisitor.h"
#include "BeamlineGeometry/Geometry.h"
#include "Fields/Fieldmap.h"
#include "PartBunch/PartBunch.h"
#include "Utilities/GeneralOpalException.h"

#include <Kokkos_Core.hpp>

#include <string>
#include <vector>

extern Inform* gmsg;

namespace {
    /// Map types this element accepts: the ones that are static in time, whether they carry
    /// a magnetic field, an electric field, or both. Dynamic (RF) maps are excluded: their
    /// readers ignore the scale arguments of applyField(), and their phase lives on the
    /// element rather than in the map, so RFCAVITY and TRAVELINGWAVE remain the elements
    /// for those.
    bool isStatic(MapType type) {
        switch (type) {
            case T1DMagnetoStatic:
            case TAstraMagnetoStatic:
            case T2DMagnetoStatic:
            case T2DMagnetoStatic_cspline:
            case TG4BL2DMagnetoStatic:
            case TG4BL3DGrid:
                return true;
            default:
                return false;
        }
    }
}  // namespace

/* ============================== Constructors ============================== */
FieldmapElement::FieldmapElement() : FieldmapElement("") {}

FieldmapElement::FieldmapElement(const FieldmapElement& right)
    : ElementBase(right),
      filename_m(right.filename_m),
      fieldmap_m(right.fieldmap_m),
      bscale_m(right.bscale_m),
      escale_m(right.escale_m),
      isZReversed_m(right.isZReversed_m),
      startField_m(right.startField_m),
      endField_m(right.endField_m),
      hasTransverseExtent_m(right.hasTransverseExtent_m),
      halfWidthX_m(right.halfWidthX_m),
      halfWidthY_m(right.halfWidthY_m) {}

FieldmapElement::FieldmapElement(const std::string& name)
    : ElementBase(name),
      filename_m(""),
      fieldmap_m(nullptr),
      bscale_m(1.0),
      escale_m(1.0),
      isZReversed_m(false),
      startField_m(0.0),
      endField_m(0.0),
      hasTransverseExtent_m(false),
      halfWidthX_m(0.0),
      halfWidthY_m(0.0) {}

FieldmapElement::~FieldmapElement() {}

/* ============================== Apply functions =========================== */
void FieldmapElement::apply(const std::shared_ptr<ParticleContainer_t>& pc) {
    if (fieldmap_m == nullptr) {
        return;
    }
    fieldmap_m->applyField(pc, bscale_m, escale_m);
}

void FieldmapElement::apply(
        const Vector_t<double, 3>& R, const Vector_t<double, 3>& /*P*/, const double& /*t*/,
        Vector_t<double, 3>& E, Vector_t<double, 3>& B) {
    if (fieldmap_m == nullptr) {
        return;
    }

    // getFieldstrength() accumulates into its output arguments despite the getter name, so
    // the temporaries have to start at zero.
    Vector_t<double, 3> tmpE(0.0), tmpB(0.0);

    if (fieldmap_m->getFieldstrength(R, tmpE, tmpB)) {
        return;  // outside the map: no field, and not an error
    }

    // The two scales are separate because G4beamline scales the two fields separately. A
    // map carrying only one of them leaves the other temporary at zero.
    B += bscale_m * tmpB;
    E += escale_m * tmpE;
}

bool FieldmapElement::applyToReferenceParticle(
        const Vector_t<double, 3>& R, const Vector_t<double, 3>& /*P*/, const double& /*t*/,
        Vector_t<double, 3>& E, Vector_t<double, 3>& B) {
    if (fieldmap_m == nullptr) {
        return false;
    }

    // A genuine aperture hit is the only thing that counts as a loss here. Leaving the map
    // transversely is not one: the particle just stops feeling a field. Reporting it as a
    // loss would make OrbitThreader flag HITMATERIAL and add a metre to the path length.
    if (R(2) >= startField_m && R(2) < endField_m
        && !ApertureHelper::isInsideAperture(R, aperture_m)) {
        return true;
    }

    Vector_t<double, 3> tmpE(0.0), tmpB(0.0);
    if (fieldmap_m->getFieldstrength(R, tmpE, tmpB)) {
        return false;
    }

    B += bscale_m * tmpB;
    E += escale_m * tmpE;
    return false;
}

/* ============================== Functions ================================= */
void FieldmapElement::accept(BeamlineVisitor& visitor) const {
    visitor.visitFieldmapElement(*this);
}

void FieldmapElement::initialise(PartBunch_t* bunch) {
    Inform msg("FieldmapElement ", *gmsg);

    RefPartBunch_m = bunch;

    if (filename_m.empty()) {
        throw GeneralOpalException(
                "FieldmapElement::initialise",
                "\"" + getName() + "\" has no field map. FMAPFN is required.");
    }

    fieldmap_m = Fieldmap::getFieldmap(filename_m, false, isZReversed_m);

    if (fieldmap_m == nullptr) {
        throw GeneralOpalException(
                "FieldmapElement::initialise",
                "\"" + getName() + "\": could not load the field map \"" + filename_m + "\".");
    }

    if (!isStatic(fieldmap_m->getType())) {
        throw GeneralOpalException(
                "FieldmapElement::initialise",
                "\"" + getName() + "\": \"" + filename_m
                        + "\" is not a static field map. A FIELDMAP element carries a field that "
                          "does not change with time, electric or magnetic or both; a "
                          "time-dependent map needs the phase and frequency that live on "
                          "RFCAVITY or TRAVELINGWAVE, so use one of those instead.");
    }

    // The reader's constructor has already parsed the header, so both extents answer now,
    // before goOnline() reads the data.
    fieldmap_m->getFieldDimensions(startField_m, endField_m);
    getGeometry().setElementLength(endField_m - startField_m);
    // The local frame is the map's frame, so the body starts where the map's z range starts.
    getGeometry().setStartZ(startField_m);

    // Three of the five readers throw here: a one-dimensional map has no transverse extent
    // at all, because its off-axis field is an expansion about the axis rather than a
    // tabulated box.
    hasTransverseExtent_m = false;
    try {
        double xIni = 0.0, xFinal = 0.0, yIni = 0.0, yFinal = 0.0, zIni = 0.0, zFinal = 0.0;
        fieldmap_m->getFieldDimensions(xIni, xFinal, yIni, yFinal, zIni, zFinal);
        halfWidthX_m          = 0.5 * std::abs(xFinal - xIni);
        halfWidthY_m          = 0.5 * std::abs(yFinal - yIni);
        hasTransverseExtent_m = halfWidthX_m > 0.0 && halfWidthY_m > 0.0;
    } catch (...) {
        halfWidthX_m = 0.0;
        halfWidthY_m = 0.0;
    }

    msg << level2 << getName() << " using file ";
    fieldmap_m->getInfo(&msg);
    msg << level2 << getName() << ": field extends over z = " << startField_m << " .. "
        << endField_m << " m in the element frame";
    if (hasTransverseExtent_m) {
        msg << ", half-widths " << halfWidthX_m << " x " << halfWidthY_m << " m" << endl;
    } else {
        msg << ". The map declares no transverse extent, so give an explicit APERTURE if the "
               "beam can leave it."
            << endl;
    }
}

void FieldmapElement::finalise() { online_m = false; }

void FieldmapElement::goOnline(const double& /*kineticEnergy*/) {
    Fieldmap::readMap(filename_m);
    online_m = true;
}

void FieldmapElement::goOffline() {
    // freeMap() drops the reference count and deletes the map when it reaches zero, so the
    // pointer must not be kept. Every use of fieldmap_m guards against null for this reason.
    Fieldmap::freeMap(filename_m);
    fieldmap_m = nullptr;
    online_m   = false;
}

ElementType FieldmapElement::getType() const { return ElementType::FIELDMAP; }

void FieldmapElement::getFieldExtent(double& zBegin, double& zEnd) const {
    zBegin = startField_m;
    zEnd   = endField_m;
}

BoundingBox FieldmapElement::getBoundingBoxInLabCoords() const {
    // Overridden only to narrow the box to the map's transverse extent. The z range is the
    // field window, which is the same as the geometry's [startZ, startZ + L].
    double x = aperture_m.second[0];
    double y = aperture_m.second[1];
    if (hasTransverseExtent_m) {
        x = std::min(x, halfWidthX_m);
        y = std::min(y, halfWidthY_m);
    }

    std::vector<Vector_t<double, 3>> corners(8);
    for (int i = -1; i < 2; i += 2) {
        for (int j = -1; j < 2; j += 2) {
            const unsigned int idx = (i + 1) / 2 + (j + 1);
            corners[idx]           = csTrafoGlobal2Local_m.transformFrom(
                    Vector_t<double, 3>({i * x, j * y, startField_m}));
            corners[idx + 4] = csTrafoGlobal2Local_m.transformFrom(
                    Vector_t<double, 3>({i * x, j * y, endField_m}));
        }
    }

    return BoundingBox::getBoundingBox(corners);
}

bool FieldmapElement::getSupportEnvelope(double& horizontalRadius, double& verticalRadius) const {
    const auto aperture = getAperture();
    if (aperture.second.size() >= 2 && std::abs(aperture.second[0]) < 1e5
        && std::abs(aperture.second[1]) < 1e5) {
        horizontalRadius = std::abs(aperture.second[0]);
        verticalRadius   = std::abs(aperture.second[1]);
        return horizontalRadius > 0.0 && verticalRadius > 0.0;
    }

    if (!hasTransverseExtent_m) {
        return false;
    }

    horizontalRadius = halfWidthX_m;
    verticalRadius   = halfWidthY_m;
    return true;
}

/* ============================== Accessors ================================= */
void FieldmapElement::setFieldMapFN(const std::string& fn) { filename_m = fn; }

const std::string& FieldmapElement::getFieldMapFN() const { return filename_m; }

void FieldmapElement::setBScale(double bscale) { bscale_m = bscale; }

double FieldmapElement::getBScale() const { return bscale_m; }

void FieldmapElement::setEScale(double escale) { escale_m = escale; }

double FieldmapElement::getEScale() const { return escale_m; }

void FieldmapElement::setIsZReversed(bool zReverse) { isZReversed_m = zReverse; }

bool FieldmapElement::getIsZReversed() const { return isZReversed_m; }

bool FieldmapElement::hasTransverseExtent() const { return hasTransverseExtent_m; }
