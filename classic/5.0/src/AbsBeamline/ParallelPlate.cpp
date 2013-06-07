// ------------------------------------------------------------------------
// $RCSfile: ParallelPlate.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ParallelPlate
//   Defines the abstract interface for Parallel Plate. For benchmarking the
//   Secondary emisssion model
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2010/10/12 15:32:31 $
// $Author: Chuan Wang $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/ParallelPlate.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Algorithms/PartBunch.h"
#include "Fields/Fieldmap.hh"
#include "Physics/Physics.h"

#include <iostream>
#include <fstream>

extern Inform *gmsg;

// Class ParallelPlate
// ------------------------------------------------------------------------

ParallelPlate::ParallelPlate():
    Component(),
    filename_m(""),
    scale_m(1.0),
    phase_m(0.0),
    frequency_m(0.0),
    length_m(0.0),
    startField_m(0.0),
    endField_m(0.0),
    ElementEdge_m(0.0),
    ptime_m(0.0) {
    setElType(isRF);
}


ParallelPlate::ParallelPlate(const ParallelPlate &right):
    Component(right),
    filename_m(right.filename_m),
    scale_m(right.scale_m),
    phase_m(right.phase_m),
    frequency_m(right.frequency_m),
    length_m(right.length_m),
    startField_m(right.startField_m),
    endField_m(right.endField_m),
    ElementEdge_m(right.ElementEdge_m),
    ptime_m(0.0) {
    setElType(isRF);
}


ParallelPlate::ParallelPlate(const string &name):
    Component(name),
    filename_m(""),
    scale_m(1.0),
    phase_m(0.0),
    frequency_m(0.0),
    length_m(0.0),
    startField_m(0.0),
    endField_m(0.0),
    ElementEdge_m(0.0),
    ptime_m(0.0) {
    setElType(isRF);
}


ParallelPlate::~ParallelPlate() {


}


void ParallelPlate::accept(BeamlineVisitor &visitor) const {
    visitor.visitParallelPlate(*this);
}


string ParallelPlate::getFieldMapFN() const {
    return "";
}

void ParallelPlate::setAmplitude(double vPeak) {
    scale_m = vPeak;
}

double ParallelPlate::getAmplitude() const {
    return scale_m ;
}

void ParallelPlate::setFrequency(double freq) {
    frequency_m = freq;
}

double ParallelPlate::getFrequency()const {
    return  frequency_m;
}

void ParallelPlate::setPhase(double phase) {
    phase_m = phase;
}

double ParallelPlate::getPhase() const {
    return phase_m;
}

void ParallelPlate::setElementLength(double length) {
    length_m = length;
}

double ParallelPlate::getElementLength() const {
    return length_m;
}


bool ParallelPlate::apply(const size_t &i, const double &t, double E[], double B[]) {
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


bool ParallelPlate::apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) {
    const double phase = frequency_m * t + phase_m;
    Vector_t tmpE(0.0, 0.0, -1.0), tmpB(0.0, 0.0, 0.0);
    if (ptime_m!=t) {
        ptime_m=t;
        //cout<<"time: "<<ptime_m<<" phase: "<<phase<<" E: "<<scale_m/length_m * sin(phase) * tmpE<<endl;

    }
    const Vector_t tmpR(RefPartBunch_m->getX(i) - dx_m, RefPartBunch_m->getY(i) - dy_m , RefPartBunch_m->getZ(i) - startField_m - ds_m);

    //Here we don't check if the particle is outof bounds, just leave the matter to the particle Boundary collision model in BoundaryGeometry
    if( (tmpR(2)>=startField_m-0.00001)
        && (tmpR(2)<=endField_m+0.00001) ) {
        E += scale_m/length_m * sin(phase) * tmpE; //Here scale_m should be Voltage between the parallel plates(V).
        B = tmpB;       //B field is always zero for our parallel plate elements used for benchmarking.
        return false;
    }
    return true;
}

bool ParallelPlate::apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B) {
    const double phase = frequency_m * t + phase_m;
    Vector_t  tmpE(0.0, 0.0, -1.0), tmpB(0.0, 0.0, 0.0);
    if (ptime_m!=t) {
        ptime_m=t;
        //cout<<"time: "<<ptime_m<<" phase: "<<phase<<" E: "<<scale_m/length_m * sin(phase) * tmpE<<endl;

    }
    const Vector_t tmpR(R(0) - dx_m, R(1) - dy_m, R(2) - startField_m - ds_m);// Fixme: check this.


    //if(!fieldmap_m->getFieldstrength(tmpR, tmpE, tmpB)) {
    if( (tmpR(2)>=startField_m-0.00001)
        && (tmpR(2)<=endField_m+0.00001) ) {
        E += scale_m/length_m * sin(phase) * tmpE;
        B = tmpB;
        return false;
    }
    return true;
}

void ParallelPlate::initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor) {

    Inform msg("PARALLELPLATE");
    std::stringstream errormsg;
    RefPartBunch_m = bunch;
    startField_m = ElementEdge_m = startField;
    //cout<<"Lenth_m in initialise: "<<length_m<<endl;
    endField_m = endField = ElementEdge_m + length_m;//fixme:Has length_m already been initialised before?

}

// In current version ,not implemented yet.
void ParallelPlate::initialise(PartBunch *bunch, const double &scaleFactor) {
    using Physics::pi;

    Inform msg("ParallelPlate initialization for cyclotron tracker ");

    RefPartBunch_m = bunch;
    // ElementEdge_m = startField;
    msg << " Currently parallelplate initialization for cyclotron tracker is  empty! " << endl;

}



void ParallelPlate::finalise()
{}

bool ParallelPlate::bends() const {
    return false;
}


void ParallelPlate::goOnline() {

    online_m = true;
}

void ParallelPlate::goOffline() {

    online_m = false;
}

void ParallelPlate::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = startField_m;
    zEnd = endField_m;
}


const string &ParallelPlate::getType() const {
    static const string type("ParallelPlate");
    return type;
}



