// ------------------------------------------------------------------------
// $RCSfile: Solenoid.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Solenoid
//   Defines the abstract interface for a solenoid magnet.
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
#include "AbsBeamline/BeamlineVisitor.h"
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
      zReverse_m(right.zReverse_m) {}

Solenoid::Solenoid(const std::string& name)
    : ElementBase(name),
      filename_m(""),
      fieldmap_m(nullptr),
      scale_m(1.0),
      scaleError_m(0.0),
      startField_m(0.0),
      endField_m(0.0),
      fast_m(true),
      zReverse_m(false) {}

Solenoid::~Solenoid() {
    //    _Fieldmap::deleteFieldmap(filename_m);
}
/* ========================================================================== */
/* ============================== Apply Functions =========================== */
/**
 * @brief apply the solenoid field to all particles in the bunch
 * @note currently not implemented
 * @returns true if at least one particle is lost, false otherwise
 * (not implemented, always returns false)
 */
void Solenoid::apply(const std::shared_ptr<ParticleContainer_t>& pc) {
    Inform m("Solenoid::apply");
    m << level5 << "Solenoid::apply() called." << endl;

    fieldmap_m->applyField(pc, scale_m + scaleError_m);
}

/**
 * @brief Apply to particle with position R and momentum P
 *
 * @param R Position
 * @param P Momentum
 * @param t Time
 * @param E Electric Field
 * @param B Magnetic Field
 */
void Solenoid::apply(
        const Vector_t<double, 3>& R, const Vector_t<double, 3>& /*P*/, const double& /*t*/,
        Vector_t<double, 3>& /*E*/, Vector_t<double, 3>& B) {
    if (R(2) >= startField_m && R(2) < endField_m) {
        Vector_t<double, 3> tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);

        const bool outOfBounds = fieldmap_m->getFieldstrength(R, tmpE, tmpB);
        if (outOfBounds) {
            return;
        }

        B += (scale_m + scaleError_m) * tmpB;
    }
}

/**
 * @brief Apply to reference particle with position R and momemtum P
 *
 * @param R Position
 * @param P Momentum
 * @param t Time
 * @param E Electric Field
 * @param B Magnetic Field
 *
 * @returns true if particle is lost, false otherwise
 */
bool Solenoid::applyToReferenceParticle(
        const Vector_t<double, 3>& R, const Vector_t<double, 3>& /*P*/, const double& /*t*/,
        Vector_t<double, 3>& /*E*/, Vector_t<double, 3>& B) {
    if (R(2) >= startField_m && R(2) < endField_m) {
        Vector_t<double, 3> tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);

        const bool outOfBounds = fieldmap_m->getFieldstrength(R, tmpE, tmpB);
        if (outOfBounds) return true;

        B += scale_m * tmpB;
    }

    return false;
}
/* ========================================================================== */
/* ============================== Functions ================================= */
/// @brief Apply visitor to Solenoid
void Solenoid::accept(BeamlineVisitor& visitor) const { visitor.visitSolenoid(*this); }

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

    fieldmap_m = Fieldmap::getFieldmap(filename_m, fast_m, zReverse_m);

    if (fieldmap_m != nullptr) {
        msg << level2 << getName() << " using file ";
        fieldmap_m->getInfo(&msg);

        double zBegin = 0.0, zEnd = 0.0;
        fieldmap_m->getFieldDimensions(zBegin, zEnd);

        startField_m = zBegin;
        endField_m   = zEnd;
    } else {
        startField_m = 0.0;
        endField_m   = 0.0;
    }
}

void Solenoid::finalise() {}

/// @brief Read field map and go online
void Solenoid::goOnline(const double&) {
    Fieldmap::readMap(filename_m);
    online_m = true;
}

/// @brief Free field map and go offline
void Solenoid::goOffline() {
    Fieldmap::freeMap(filename_m);
    online_m = false;
}

/// @brief Assign the field filename
void Solenoid::setFieldMapFN(std::string fn) { filename_m = fn; }

/// @brief Set the fast flag
void Solenoid::setFast(bool fast) { fast_m = fast; }

/// @brief Get the fast flag
bool Solenoid::getFast() const { return fast_m; }

/// @brief Set the flag that reads the field map back to front
void Solenoid::setZReverse(bool zReverse) { zReverse_m = zReverse; }

/// @brief Get the flag that reads the field map back to front
bool Solenoid::getZReverse() const { return zReverse_m; }

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
    return fieldmap_m != nullptr && ApertureHelper::isInsideAperture(r, aperture_m)
           && fieldmap_m->isInside(r);
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
