// ------------------------------------------------------------------------
// $RCSfile: Component.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Component
//   An abstract base class which defines the common interface for all
//   CLASSIC components, i.e. beamline members which are not themselves
//   beamlines.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//c
// ------------------------------------------------------------------------

#include "AbsBeamline/Component.h"
#include "Utilities/LogicalError.h"
#include "Algorithms/PartBunch.h"

// Class Component
// ------------------------------------------------------------------------
//   Represents an arbitrary component in an accelerator.  A component is
//   the basic element in the accelerator model, and can be thought of as
//   acting as a leaf in the Composite pattern.  A Component is associated
//   with an electromagnetic field.

Component::Component():
    ElementBase(),
    exit_face_slope_m(0.0),
    RefPartBunch_m(NULL),
    online_m(false)
{}


Component::Component(const Component &right):
    ElementBase(right),
    exit_face_slope_m(right.exit_face_slope_m),
    RefPartBunch_m(right.RefPartBunch_m),
    online_m(right.online_m)
{}


Component::Component(const std::string &name):
    ElementBase(name),
    exit_face_slope_m(0.0),
    RefPartBunch_m(NULL),
    online_m(false)
{}


Component::~Component()
{}


const ElementBase &Component::getDesign() const {
    return *this;
}

void Component::trackBunch(PartBunch &, const PartData &, bool, bool) const {
    throw LogicalError("Component::trackBunch()",
                       "Called for component \"" + getName() + "\".");
}


void Component::
trackMap(FVps<double, 6> &, const PartData &, bool, bool) const {
    throw LogicalError("Component::trackMap()",
                       "Called for component \"" + getName() + "\".");
}

void Component::goOnline(const double &) {
    online_m = true;
}

void Component::goOffline() {
    online_m = false;
}

bool Component::Online() {
    return online_m;
}

ElementBase::ElementType Component::getType() const {
    return ElementBase::ANY;
}

bool Component::apply(const size_t &i,
                      const double &t,
                      Vector_t &E,
                      Vector_t &B) {
    const Vector_t &R = RefPartBunch_m->R[i];
    if (R(2) >= 0.0 && R(2) < getElementLength()) {
        if (!isInsideTransverse(R)) return true;
    }
    return false;
}

bool Component::apply(const Vector_t &R,
                      const Vector_t &P,
                      const double &t,
                      Vector_t &E,
                      Vector_t &B) {
    if (R(2) >= 0.0 && R(2) < getElementLength()) {
        if (!isInsideTransverse(R)) return true;
    }
    return false;
}

bool Component::applyToReferenceParticle(const Vector_t &R,
                                         const Vector_t &P,
                                         const double &t,
                                         Vector_t &E,
                                         Vector_t &B) {
    if (R(2) >= 0.0 && R(2) < getElementLength()) {
        if (!isInsideTransverse(R)) return true;
    }
    return false;
}
