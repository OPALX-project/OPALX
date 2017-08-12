// ------------------------------------------------------------------------
// $RCSfile: Degrader.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Degrader
//   Defines the abstract interface for a beam Degrader.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Degrader.h"
#include "Algorithms/PartBunch.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Fields/Fieldmap.h"
#include "Structure/LossDataSink.h"
#include "Utilities/Options.h"
#include "Physics/Physics.h"
#include <memory>

extern Inform *gmsg;

using namespace std;

// Class Degrader
// ------------------------------------------------------------------------

Degrader::Degrader():
    Component(),
    position_m(0.0),
    semiMinorAxis_m(1e6),
    semiMajorAxis_m(1e6)
{}

Degrader::Degrader(const Degrader &right):
    Component(right),
    position_m(right.position_m),
    semiMinorAxis_m(right.semiMinorAxis_m),
    semiMajorAxis_m(right.semiMajorAxis_m)
{}

Degrader::Degrader(const std::string &name):
    Component(name),
    position_m(0.0),
    semiMinorAxis_m(1e6),
    semiMajorAxis_m(1e6)
{}


Degrader::~Degrader() {

    if(online_m)
        goOffline();
}


void Degrader::accept(BeamlineVisitor &visitor) const {
    visitor.visitDegrader(*this);
}


bool Degrader::apply(const size_t &i, const double &t, double E[], double B[]) {
    Vector_t Ev(0, 0, 0), Bv(0, 0, 0);
    return apply(i, t, Ev, Bv);
}

bool Degrader::apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) {

    const Vector_t &R = RefPartBunch_m->R[i] - Vector_t(dx_m, dy_m, ds_m); // including the missaligment

    if (isInMaterial(R)) {
        //if particle was allready marked as -1 (it means it should have gone into degrader but didn't)
        //set the label to -2 (will not go into degrader and will be deleted when particles per core > 2)
        if (RefPartBunch_m->Bin[i] < 0)
            RefPartBunch_m->Bin[i] = -2;
        else
            RefPartBunch_m->Bin[i] = -1;
    }

    return false;
}

bool Degrader::apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B) {
    return false;
}


void Degrader::initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor) {
    if (!hasSurfacePhysics())
        throw GeneralClassicException("Degrader::initialize",
                                      "No surface physics attached");
    RefPartBunch_m = bunch;
    position_m = startField;
    endField = position_m + getElementLength();
}

void Degrader::initialise(PartBunch *bunch, const double &scaleFactor) {
    RefPartBunch_m = bunch;
}


void Degrader::finalise()
{
    *gmsg << "* Finalize Degrader" << endl;
}

void Degrader::goOnline(const double &) {
    Inform msg("Degrader::goOnline ");
    if(RefPartBunch_m == NULL) {
        throw GeneralClassicException("Degrader::goOnline",
                                      "'" + getName() + "' isn't initialized");
    }
    online_m = true;
}

void Degrader::goOffline() {
    Inform msg("Degrader::goOffline ");
    online_m = false;

    msg << " done..." << endl;
}

bool Degrader::bends() const {
    return false;
}

void Degrader::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = position_m;
    zEnd = position_m + getElementLength();

}

ElementBase::ElementType Degrader::getType() const {
    return DEGRADER;
}