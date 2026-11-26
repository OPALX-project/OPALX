// ------------------------------------------------------------------------
// $RCSfile: Solenoid.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Solenoid
//   Defines a solenoid magnet.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 20.1.2026 $
// $Author: PSI $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Solenoid.h"
#include "Fields/Fieldmap.h"
#include "PartBunch/PartBunch.h"
#include "Physics/Physics.h"

#include <cmath>
#include <fstream>
#include <iostream>

extern Inform* gmsg;

/* ============================== Constructors ============================== */
Solenoid::Solenoid() : Solenoid("") {}

Solenoid::Solenoid(const Solenoid& right)
    : ElementBase(right),
      filename_m(right.filename_m),
      fieldmap_m(right.fieldmap_m),
      scale_m(right.scale_m),
      scaleError_m(right.scaleError_m),
      startField_m(right.startField_m),
      endField_m(right.endField_m),
      fast_m(right.fast_m),
      isZReversed_m(right.isZReversed_m) {}

Solenoid::Solenoid(const std::string& name)
    : ElementBase(name),
      filename_m(""),
      fieldmap_m(nullptr),
      scale_m(1.0),
      scaleError_m(0.0),
      startField_m(0.0),
      endField_m(0.0),
      fast_m(true),
      isZReversed_m(false) {}

Solenoid::~Solenoid() {
    //    _Fieldmap::deleteFieldmap(filename_m);
}
/* ========================================================================== */
/* ============================== Apply Functions =========================== */
/**
 * @brief Apply the solenoid field to all particles in the bunch.
 *
 * The map gives B times KS + DKS, and its E unscaled if it has one.
 */
void Solenoid::apply(const std::shared_ptr<ParticleContainer_t>& pc, const StepTime& /*time*/) {
    if (fieldmap_m == nullptr) return;
    fieldmap_m->applyField(pc, scale_m + scaleError_m, 1.0, getRegions());
}

/**
 * @brief Add the solenoid field at local position R (host).
 *
 * Same gate, map and scales as apply(pc).
 */
void Solenoid::applyOnHost(
        const Vector_t<double, 3>& R, double /*t*/, Vector_t<double, 3>& E,
        Vector_t<double, 3>& B) {
    if (fieldmap_m == nullptr || !getRegions().inFieldRegion(R)) return;

    Vector_t<double, 3> tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);
    // Outside the map data the map gives no field.
    if (fieldmap_m->getFieldstrength(R, tmpE, tmpB)) return;

    B += (scale_m + scaleError_m) * tmpB;
    E += tmpE;
}
/* ========================================================================== */
/* ============================== Functions ================================= */
/// @brief Set the strength scaling factor ks
void Solenoid::setKS(double ks) { scale_m = ks; }

/// @brief Set the strength scaling error dks
void Solenoid::setDKS(double ks) { scaleError_m = ks; }

/**
 * @brief Initialises the solenoid elements by reading the header of the
 * fieldmap, saving the start and endpoints of the fieldmaps.
 *
 * @param bunch Particle bunch
 */
void Solenoid::initialise(PartBunch_t* bunch) {
    Inform msg("Solenoid ", *gmsg);

    RefPartBunch_m = bunch;

    fieldmap_m = Fieldmap::getFieldmap(filename_m, fast_m, isZReversed_m);

    if (fieldmap_m != nullptr) {
        switch (fieldmap_m->getType()) {
            case T1DDynamic:
            case TAstraDynamic:
            case T2DDynamic:
            case T2DDynamic_cspline:
            case T3DDynamic:
            case T3DDynamicH5Block:
                // A time-dependent map needs a frequency and a phase, which a SOLENOID does not
                // have; its readers would also ignore KS.
                throw GeneralOpalException(
                        "Solenoid::initialise",
                        "SOLENOID '" + getName() + "': the field map '" + filename_m
                                + "' is time dependent. Use RFCAVITY or FIELDMAP for it.");
            default:
                break;
        }
        msg << level2 << getName() << " using file ";
        fieldmap_m->getInfo(&msg);

        double zBegin = 0.0, zEnd = 0.0;
        fieldmap_m->getFieldDimensions(zBegin, zEnd);

        startField_m = zBegin;
        endField_m   = zEnd;
        Fieldmap::readMap(filename_m);
    } else {
        startField_m = 0.0;
        endField_m   = 0.0;
    }
}

/// @brief Free field map and go offline
void Solenoid::goOffline() {
    // freeMap() deletes the map when no element uses it any more: do not keep the pointer.
    Fieldmap::freeMap(filename_m);
    fieldmap_m = nullptr;
    online_m   = false;
}

/// @brief Assign the field filename
void Solenoid::setFieldMapFN(std::string fn) { filename_m = fn; }

/// @brief Set the fast flag
void Solenoid::setFast(bool fast) { fast_m = fast; }

/// @brief Set the flag that reads the field map back to front
void Solenoid::setIsZReversed(bool zReverse) { isZReversed_m = zReverse; }

/// @brief Get the flag that reads the field map back to front
bool Solenoid::getIsZReversed() const { return isZReversed_m; }

/// @brief Get the dimensions of the solenoid
/// @param zBegin Start position
/// @param zEnd End position
void Solenoid::getFieldExtent(double& zBegin, double& zEnd) const {
    zBegin = startField_m;
    zEnd   = endField_m;
}

/// @brief Get the element type
/// @returns ElementType::SOLENOID
ElementType Solenoid::getType() const { return ElementType::SOLENOID; }

/// @brief Check if a point is inside the solenoid
/// @param r Position
/// @returns true if inside, false otherwise
bool Solenoid::isInside(const Vector_t<double, 3>& r) const {
    return fieldmap_m != nullptr && ElementBase::isInside(r) && fieldmap_m->isInside(r);
}

std::optional<ApertureShape> Solenoid::getFieldMapRange() const {
    if (fieldmap_m == nullptr) return std::nullopt;
    return fieldmap_m->getTransverseRange();
}

/// @brief Get the dimensions of the solenoid
/// @param begin Start position
/// @param end End position
bool Solenoid::getSupportEnvelope(double& horizontalRadius, double& verticalRadius) const {
    const auto aperture = getAperture();
    if (aperture.second.size() >= 2 && std::abs(aperture.second[0]) < 1e5
        && std::abs(aperture.second[1]) < 1e5) {
        horizontalRadius = std::abs(aperture.second[0]);
        verticalRadius   = std::abs(aperture.second[1]);
        return horizontalRadius > 0.0 && verticalRadius > 0.0;
    }

    if (fieldmap_m == nullptr) {
        return false;
    }

    try {
        double xIni = 0.0, xFinal = 0.0, yIni = 0.0, yFinal = 0.0, zIni = 0.0, zFinal = 0.0;
        fieldmap_m->getFieldDimensions(xIni, xFinal, yIni, yFinal, zIni, zFinal);
        horizontalRadius = 0.5 * std::abs(xFinal - xIni);
        verticalRadius   = 0.5 * std::abs(yFinal - yIni);
        return horizontalRadius > 0.0 && verticalRadius > 0.0;
    } catch (...) {
        return false;
    }
}
