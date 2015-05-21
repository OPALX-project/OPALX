// ------------------------------------------------------------------------
// $RCSfile: ElementBase.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ElementBase
//   The very base class for beamline representation objects.  A beamline
//   is modelled as a composite structure having a single root object
//   (the top level beamline), which contains both "single" leaf-type
//   elements (Components), as well as sub-lines (composites).
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/12/16 16:26:43 $
// $Author: mad $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/AlignWrapper.h"
#include "AbsBeamline/ElementImage.h"
#include "Channels/Channel.h"
#include <string>

#include "Structure/BoundaryGeometry.h"    // OPAL file
#include "Solvers/WakeFunction.hh"
#include "Solvers/SurfacePhysicsHandler.hh"

using namespace std;

// Class ElementBase
// ------------------------------------------------------------------------

ElementBase::ElementBase():
    RCObject(),
    shareFlag(true),
    elementID(""),
    userAttribs(),
    wake_m(NULL),
    bgeometry_m(NULL),
    sphys_m(NULL)
{}


ElementBase::ElementBase(const ElementBase &right):
    RCObject(),
    shareFlag(true),
    elementID(right.elementID),
    userAttribs(right.userAttribs),
    wake_m(right.wake_m),
    bgeometry_m(right.bgeometry_m),
    sphys_m(right.sphys_m) {

    if(sphys_m) {
        sphys_m->updateElement(this);
    }
    if(bgeometry_m)
        bgeometry_m->updateElement(this);
}


ElementBase::ElementBase(const std::string &name):
    RCObject(),
    shareFlag(true),
    elementID(name),
    userAttribs(),
    wake_m(NULL),
    bgeometry_m(NULL),
    sphys_m(NULL)
{}


ElementBase::~ElementBase()

{}


const std::string &ElementBase::getName() const

{
    return elementID;
}


void ElementBase::setName(const std::string &name) {
    elementID = name;
}


double ElementBase::getAttribute(const std::string &aKey) const {
    const ConstChannel *aChannel = getConstChannel(aKey);

    if(aChannel != NULL) {
        double val = *aChannel;
        delete aChannel;
        return val;
    } else {
        return 0.0;
    }
}


bool ElementBase::hasAttribute(const std::string &aKey) const {
    const ConstChannel *aChannel = getConstChannel(aKey);

    if(aChannel != NULL) {
        delete aChannel;
        return true;
    } else {
        return false;
    }
}


void ElementBase::removeAttribute(const std::string &aKey) {
    userAttribs.removeAttribute(aKey);
}


void ElementBase::setAttribute(const std::string &aKey, double val) {
    Channel *aChannel = getChannel(aKey, true);

    if(aChannel != NULL  &&  aChannel->isSettable()) {
        *aChannel = val;
        delete aChannel;
    } else
        cout << "Channel NULL or not Settable" << endl;
}


Channel *ElementBase::getChannel(const std::string &aKey, bool create) {
    return userAttribs.getChannel(aKey, create);
}


const ConstChannel *ElementBase::getConstChannel(const std::string &aKey) const {
    // Use const_cast to allow calling the non-const method GetChannel().
    // The const return value of this method will nevertheless inhibit set().
    return const_cast<ElementBase *>(this)->getChannel(aKey);
}


std::string ElementBase::getTypeString() const {
    switch (getType()) {
    case ALIGNWRAPPER:
        return "AlignWrapper";
    case BEAMBEAM:
        return "BeamBeam";
    case BEAMLINE:
        return "Beamline";
    case COLLIMATOR:
        return "Collimator";
    case CORRECTOR:
        return "Corrector";
    case CYCLOTRON:
        return "Cyclotron";
    case CYCLOTRONVALLEY:
        return "CyclotronValley";
    case DEGRADER:
        return "Degrader";
    case DIAGNOSTIC:
        return "Diagnostic";
    case DRIFT:
        return "Drift";
    case INTEGRATOR:
        return "Integrator";
    case LAMBERTSON:
        return "Lambertson";
    case MARKER:
        return "Marker";
    case MONITOR:
        return "Monitor";
    case MULTIPOLE:
        return "Multipole";
    case OFFSET:
        return "Offset";
    case PARALLELPLATE:
        return "ParallelPlate";
    case PATCH:
        return "Patch";
    case PROBE:
        return "Probe";
    case RBEND:
        return "RBend";
    case RFCAVITY:
        return "RFCavity";
    case RFQUADRUPOLE:
        return "RFQuadrupole";
    case RING:
        return "Ring";
    case SBEND3D:
        return "SBend3D";
    case SBEND:
        return "SBend";
    case SEPARATOR:
        return "Separator";
    case SEPTUM:
        return "Septum";
    case SOLENOID:
        return "Solenoid";
    case STRIPPER:
        return "Stripper";
    case TRAVELINGWAVE:
        return "TravelingWave";
    case VARIABLERFCAVITY:
        return "VariableRFCavity";
    case ANY:
    default:
        return "";
    }
}

ElementImage *ElementBase::getImage() const {
    std::string type = getTypeString();
    return new ElementImage(getName(), type, userAttribs);
}


ElementBase *ElementBase::copyStructure() {
    if(isSharable()) {
        return this;
    } else {
        return clone();
    }
}


void ElementBase::makeSharable() {
    shareFlag = true;
}


ElementBase *ElementBase::makeAlignWrapper() {
    ElementBase *wrap = new AlignWrapper(this);
    wrap->setName(getName());
    return wrap;
}


ElementBase *ElementBase::makeFieldWrapper() {
    return this;
}


ElementBase *ElementBase::makeWrappers() {
    return makeFieldWrapper()->makeAlignWrapper();
}


ElementBase *ElementBase::removeAlignWrapper() {
    return this;
}


const ElementBase *ElementBase::removeAlignWrapper() const {
    return this;
}


ElementBase *ElementBase::removeFieldWrapper() {
    return this;
}


const ElementBase *ElementBase::removeFieldWrapper() const {
    return this;
}


ElementBase *ElementBase::removeWrappers() {
    return this;
}


const ElementBase *ElementBase::removeWrappers() const {
    return this;
}


bool ElementBase::update(const AttributeSet &set) {
    for(AttributeSet::const_iterator i = set.begin(); i != set.end(); ++i) {
        setAttribute(i->first, i->second);
    }

    return true;
}

void ElementBase::setWake(WakeFunction *wk) {
    wake_m = wk;//->clone(getName() + std::string("_wake")); }
}

void ElementBase::setBoundaryGeometry(BoundaryGeometry *geo) {
    bgeometry_m = geo;//->clone(getName() + std::string("_wake")); }
}

void ElementBase::setSurfacePhysics(SurfacePhysicsHandler *sphys) {
    sphys_m = sphys;
}