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

#include <fstream>
#include <iostream>

extern Inform* gmsg;

/* ============================== Constructors ============================== */
Solenoid::Solenoid() : Solenoid("") {}

Solenoid::Solenoid(const Solenoid& right)
    : Component(right),
      filename_m(right.filename_m),
      fieldmap_m(right.fieldmap_m),
      scale_m(right.scale_m),
      scaleError_m(right.scaleError_m),
      startField_m(right.startField_m),
      fast_m(right.fast_m) {}

Solenoid::Solenoid(const std::string& name)
    : Component(name),
      filename_m(""),
      fieldmap_m(nullptr),
      scale_m(1.0),
      scaleError_m(0.0),
      startField_m(0.0),
      fast_m(true) {}

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
bool Solenoid::apply() {
    std::cout << "Solenoid::apply() called" << std::endl;

    std::shared_ptr<ParticleContainer_t> pc = RefPartBunch_m->getParticleContainer();

    fieldmap_m->applyField(pc);

    return false;
}

/**
 * @brief apply the solenoid field to particle i
 *
 * @param i Particle index
 * @param t Time
 * @param E Electric Field
 * @param B Magnetic Field
 *
 * @returns true if particle is lost, false otherwise
 */
bool Solenoid::apply(
    const size_t& i, const double& t, Vector_t<double, 3>& E, Vector_t<double, 3>& B
) {
    std::shared_ptr<ParticleContainer_t> pc = RefPartBunch_m->getParticleContainer();
    auto Rview                              = pc->R.getView();
    auto Pview                              = pc->P.getView();

    const Vector_t<double, 3> R = Rview(i);
    const Vector_t<double, 3> P = Pview(i);
    return apply(R, P, t, E, B);
}

/**
 * @brief Apply to particle with position R and momentum P
 *
 * @param R Position
 * @param P Momentum
 * @param t Time
 * @param E Electric Field
 * @param B Magnetic Field
 *
 * @returns true if particle is lost, false otherwise
 */
bool Solenoid::apply(
    const Vector_t<double, 3>& R, const Vector_t<double, 3>& /*P*/, const double& /*t*/,
    Vector_t<double, 3>& /*E*/, Vector_t<double, 3>& B
) {
    if (R(2) >= startField_m && R(2) < startField_m + getElementLength()) {
        Vector_t<double, 3> tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);

        const bool outOfBounds = fieldmap_m->getFieldstrength(R, tmpE, tmpB);
        if (outOfBounds) {
            return getFlagDeleteOnTransverseExit();
        }

        B += (scale_m + scaleError_m) * tmpB;
    }

    return false;
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
    Vector_t<double, 3>& /*E*/, Vector_t<double, 3>& B
) {
    if (R(2) >= startField_m && R(2) < startField_m + getElementLength()) {
        Vector_t<double, 3> tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);

        const bool outOfBounds = fieldmap_m->getFieldstrength(R, tmpE, tmpB);
        if (outOfBounds)
            return true;

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
 * @param startField Starting position of the field
 * @param endField Ending position of the field
 */
void Solenoid::initialise(PartBunch_t* bunch, double& startField, double& endField) {
    Inform msg("Solenoid ", *gmsg);

    RefPartBunch_m = bunch;

    fieldmap_m = Fieldmap::getFieldmap(filename_m, fast_m);

    if (fieldmap_m != nullptr) {
        msg << level2 << getName() << " using file ";
        fieldmap_m->getInfo(&msg);

        double zBegin = 0.0, zEnd = 0.0;
        fieldmap_m->getFieldDimensions(zBegin, zEnd);

        startField_m = zBegin;
        setElementLength(zEnd - zBegin);
        endField = startField + getElementLength();
    } else {
        endField = startField;
    }
}

void Solenoid::finalise() {}

bool Solenoid::bends() const { return false; }

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

/// @brief Get the dimensions of the solenoid
/// @param zBegin Start position
/// @param zEnd End position
void Solenoid::getDimensions(double& zBegin, double& zEnd) const {
    zBegin = startField_m;
    zEnd   = startField_m + getElementLength();
}

/// @brief Get the element type
/// @returns ElementType::SOLENOID
ElementType Solenoid::getType() const { return ElementType::SOLENOID; }

/// @brief Check if a point is inside the solenoid
/// @param r Position
/// @returns true if inside, false otherwise
bool Solenoid::isInside(const Vector_t<double, 3>& r) const {
    return isInsideTransverse(r) && fieldmap_m->isInside(r);
}

/// @brief Get the dimensions of the solenoid
/// @param begin Start position
/// @param end End position
void Solenoid::getElementDimensions(double& begin, double& end) const {
    begin = startField_m;
    end   = begin + getElementLength();
}
