// ------------------------------------------------------------------------
// $RCSfile: CyclotronValley.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: CyclotronValley
//   Defines the abstract interface for an accelerating structure.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/CyclotronValley.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Algorithms/PartBunch.h"
#include "Fields/Fieldmap.hh"
#include "Physics/Physics.h"

#include <iostream>
#include <fstream>

extern Inform *gmsg;

// Class CyclotronValley
// ------------------------------------------------------------------------

CyclotronValley::CyclotronValley():
    Component(),
    filename_m(""),
    fieldmap_m(NULL),
    scale_m(1.0),
    ElementEdge_m(0.0),
    startField_m(0.0),
    endField_m(0.0),
    fast_m(false) {
    setElType(isRF);
}


CyclotronValley::CyclotronValley(const CyclotronValley &right):
    Component(right),
    filename_m(right.filename_m),
    fieldmap_m(right.fieldmap_m),
    scale_m(right.scale_m),
    ElementEdge_m(right.ElementEdge_m),
    startField_m(right.startField_m),
    endField_m(right.endField_m),
    fast_m(right.fast_m) {
    setElType(isRF);
}


CyclotronValley::CyclotronValley(const std::string &name):
    Component(name),
    filename_m(""),
    fieldmap_m(NULL),
    scale_m(1.0),
    ElementEdge_m(0.0),
    startField_m(0.0),
    endField_m(0.0),
    fast_m(false) {
    setElType(isRF);
}


CyclotronValley::~CyclotronValley() {
  //    Fieldmap::deleteFieldmap(filename_m);
    /*if(RNormal_m) {
        delete[] RNormal_m;
        delete[] VrNormal_m;
        delete[] DvDr_m;
        }*/
}


void CyclotronValley::accept(BeamlineVisitor &visitor) const {
    visitor.visitCyclotronValley(*this);
}

void CyclotronValley::setFieldMapFN(std::string fn) {
    filename_m = fn;
}

std::string CyclotronValley::getFieldMapFN() const {
    return filename_m;
}

void CyclotronValley::setFast(bool fast) {
    fast_m = fast;
}


bool CyclotronValley::getFast() const {
    return fast_m;
}

bool CyclotronValley::apply(const size_t &i, const double &t, double E[], double B[]) {
    Vector_t Ev(0, 0, 0), Bv(0, 0, 0);
    Vector_t Rt(RefPartBunch_m->getX(i), RefPartBunch_m->getY(i), RefPartBunch_m->getZ(i));
    if(apply(Rt, Vector_t(0.0), t, Ev, Bv)) return true;

    E[0] = Ev(0);
    E[1] = Ev(1);
    E[2] = Ev(2);
    B[0] = Bv(0);
    B[1] = Bv(1);
    B[2] = Bv(2);

    return false;
}

bool CyclotronValley::apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) {


    const Vector_t tmpR(RefPartBunch_m->getX(i) - dx_m, RefPartBunch_m->getY(i) - dy_m , RefPartBunch_m->getZ(i) - startField_m - ds_m);
    Vector_t tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);
    if(!fieldmap_m->getFieldstrength(tmpR, tmpE, tmpB)) {
        // E += scale_m * cos(phase) * tmpE;
        //B -= scale_m * sin(phase) * tmpB;
	B +=scale_m * tmpB;
        return false;
    }

    return true;
}

bool CyclotronValley::apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B) {


    const Vector_t tmpR(R(0) - dx_m, R(1) - dy_m, R(2) - startField_m - ds_m);
    Vector_t  tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);

    if(!fieldmap_m->getFieldstrength(tmpR, tmpE, tmpB)) {
        // E += scale_m * cos(phase) * tmpE;
        //B -= scale_m * sin(phase) * tmpB;
	B +=scale_m * tmpB;
        return false;
    }
    return true;
}

void CyclotronValley::initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor) {// called by ParallelTTracker::visitCyclotronValley --> OpalBeamline::visit
    using Physics::two_pi;
    double zBegin = 0.0, zEnd = 0.0, rBegin = 0.0, rEnd = 0.0;
    Inform msg("CyclotronValley ");
    std::stringstream errormsg;
    RefPartBunch_m = bunch;

    fieldmap_m = Fieldmap::getFieldmap(filename_m, fast_m);
    if(fieldmap_m != NULL) {
        fieldmap_m->getFieldDimensions(zBegin, zEnd, rBegin, rEnd);
        if(zEnd > zBegin) {
            msg << getName() << " using file: " << "  ";
            fieldmap_m->getInfo(&msg);

            ElementEdge_m = startField;
            startField_m = startField = ElementEdge_m + zBegin;
            endField_m = endField = ElementEdge_m + zEnd;

        } else {
            endField = startField - 1e-3;
        }
    } else {
        endField = startField - 1e-3;
    }



}


void CyclotronValley::finalise()
{}

bool CyclotronValley::bends() const {
    return false;
}


void CyclotronValley::goOnline(const double &) {
    Fieldmap::readMap(filename_m);
    online_m = true;
}

void CyclotronValley::goOffline() {
    Fieldmap::freeMap(filename_m);
    online_m = false;
}





void CyclotronValley::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = startField_m;
    zEnd = endField_m;
}


const std::string &CyclotronValley::getType() const {
    static const std::string type("CyclotronValley");
    return type;
}
