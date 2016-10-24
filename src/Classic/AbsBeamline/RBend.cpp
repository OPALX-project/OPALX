// ------------------------------------------------------------------------
// $RCSfile: RBend.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RBend
//   Defines the abstract interface for a rectangular bend magnet.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/RBend.h"
#include "Algorithms/PartBunch.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Utilities/Options.h"
#include "Fields/Fieldmap.h"
#include "AbstractObjects/OpalData.h"
#include <iostream>
#include <fstream>

extern Inform *gmsg;

// Class RBend
// ------------------------------------------------------------------------

RBend::RBend():
    Bend()
{
    setMessageHeader("RBend ");
}

RBend::RBend(const RBend &right):
    Bend(right)
{
    setMessageHeader("RBend ");
}

RBend::RBend(const std::string &name):
    Bend(name)
{
    setMessageHeader("RBend ");
}

RBend::~RBend() {
}

void RBend::accept(BeamlineVisitor &visitor) const {
    visitor.visitRBend(*this);
}

/*
 * OPAL-MAP methods
 * ================
 */
double RBend::getNormalComponent(int n) const {
    return getField().getNormalComponent(n);
}

double RBend::getSkewComponent(int n) const {
    return getField().getSkewComponent(n);
}

void RBend::setNormalComponent(int n, double v) {
    getField().setNormalComponent(n, v);
}

void RBend::setSkewComponent(int n, double v) {
    getField().setSkewComponent(n, v);
}

/*
 * BET methods.
 */
void RBend::addKR(int i, double t, Vector_t &K) {
    Inform msg("RBend::addK()");

    Vector_t tmpE(0.0, 0.0, 0.0);
    Vector_t tmpB(0.0, 0.0, 0.0);
    Vector_t tmpE_diff(0.0, 0.0, 0.0);
    Vector_t tmpB_diff(0.0, 0.0, 0.0);
    double pz = RefPartBunch_m->getZ(i) - getStartField() - ds_m;
    const Vector_t tmpA(RefPartBunch_m->getX(i) - dx_m, RefPartBunch_m->getY(i) - dy_m, pz);

    DiffDirection zdir(DZ);
    getFieldmap()->getFieldstrength(tmpA, tmpE, tmpB);
    getFieldmap()->getFieldDerivative(tmpA, tmpE_diff, tmpB_diff, zdir);

    double g = RefPartBunch_m->getGamma(i);

    if(std::abs(getFieldAmplitude() * tmpB_diff(2)) > 0.1) {
        double cf = Physics::q_e * tmpB(2) / (g * Physics::EMASS);
        K += Vector_t(-pow(cf * getFieldAmplitude() * tmpB(0), 2) / 3.0, -pow(cf * getFieldAmplitude() * tmpB(1), 2) / 3.0, 0.0);
    }
}

void RBend::addKT(int i, double t, Vector_t &K) {
    Inform msg("RBend::addK()");

    Vector_t tmpE(0.0, 0.0, 0.0);
    Vector_t tmpB(0.0, 0.0, 0.0);
    double pz = RefPartBunch_m->getZ(i) - getStartField() - ds_m;
    const Vector_t tmpA(RefPartBunch_m->getX(i) - dx_m, RefPartBunch_m->getY(i) - dy_m, pz);
    getFieldmap()->getFieldstrength(tmpA, tmpE, tmpB);

    double b = RefPartBunch_m->getBeta(i);
    double g = 1 / sqrt(1 - b * b);

    double cf = -Physics::q_e * Physics::c * b * tmpB(2) * getFieldAmplitude() / (g * Physics::EMASS);
    Vector_t temp(cf * tmpB(1), cf * tmpB(0), 0.0);

    //FIXME: K += ??
}


/*
 * OPAL-T Methods.
 * ===============
 */

/*
 *  This function merely repackages the field arrays as type Vector_t and calls
 *  the equivalent method but with the Vector_t data types.
 */

ElementBase::ElementType RBend::getType() const {
    return RBEND;
}

void RBend::SetBendAngle(double angle) {
    getAngle() = angle;
    getExitAngle() = std::abs(angle) - getEntranceAngle();
}

void RBend::SetEntranceAngle(double entranceAngle) {
    getEntranceAngle() = entranceAngle;
    getExitAngle() = std::abs(getAngle()) - entranceAngle;
}

bool RBend::FindChordLength(Inform &msg,
                            double &chordLength,
                            bool &chordLengthFromMap) {

    /*
     * Find bend chord length. If this was not set by the user using the
     * L (length) attribute, infer it from the field map.
     */
    chordLength = 2 * getLength() * sin(0.5 * std::abs(getAngle())) /
        (sin(getEntranceAngle()) + sin(std::abs(getAngle()) - getEntranceAngle()));
    if(chordLength > 0.0) {
        chordLengthFromMap = false;
        return true;
    } else {

        if(chordLength == 0.0) {
            double length = getMapLength();
            chordLength = 2 * length * sin(0.5 * getAngle()) /
                (sin(getEntranceAngle()) + sin(getAngle() - getEntranceAngle()));
        }
        chordLengthFromMap = true;

        if(chordLength <= 0.0) {
            ERRORMSG(level2 << "Magnet length inferred from field map is less than or equal"
                     << " to zero. Check your bend magnet input." << endl);
            return false;
        } else
            return true;

    }
}