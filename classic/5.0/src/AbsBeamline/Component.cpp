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
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Component.h"
#include "Utilities/LogicalError.h"


// Class Component
// ------------------------------------------------------------------------
//   Represents an arbitrary component in an accelerator.  A component is
//   the basic element in the accelerator model, and can be thought of as
//   acting as a leaf in the Composite pattern.  A Component is associated
//   with an electromagnetic field.

Component::Component():
    ElementBase(),
    dx_m(0.0),
    dy_m(0.0),
    ds_m(0.0),
    Orientation_m(Vector_t(0.0, 0.0, 0.0)),
    exit_face_slope_m(0.0),
    RefPartBunch_m(NULL),
    online_m(false)
{}


Component::Component(const Component &right):
    ElementBase(right),
    dx_m(right.dx_m),
    dy_m(right.dy_m),
    ds_m(right.ds_m),
    Orientation_m(right.Orientation_m),
    exit_face_slope_m(right.exit_face_slope_m),
    RefPartBunch_m(right.RefPartBunch_m),
    online_m(right.online_m)
{}


Component::Component(const string &name):
    ElementBase(name),
    dx_m(0.0),
    dy_m(0.0),
    ds_m(0.0),
    Orientation_m(Vector_t(0.0, 0.0, 0.0)),
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

void Component::initialise(PartBunch *bunch, const double &scaleFactor)
{}

void Component::initialise(PartBunch *bunch, const int &fieldflag, const double &scaleFactor)
{}


void Component::goOnline() {
    online_m = true;
}

void Component::goOffline() {
    online_m = false;
}

bool Component::Online() {
    return online_m;
}

void Component::setMisalignment(double x, double y, double s) {
    dx_m = x;
    dy_m = y;
    ds_m = s;
}

void Component::getMisalignment(double &x, double &y, double &s) const {
    x = dx_m;
    y = dy_m;
    s = ds_m;
}

const string &Component::getType() const {
    static const string type("");
    return type;
}
