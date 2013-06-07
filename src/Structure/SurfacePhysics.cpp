// ------------------------------------------------------------------------
// $RCSfile: SurfacePhysics.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: SurfacePhysics
//   The class for the OPAL SURFACEPHYSICS command.
//
// $Date: 2009/07/14 22:09:00 $
// $Author: C. Kraus $
//
// ------------------------------------------------------------------------

#include "Structure/SurfacePhysics.h"
#include "Solvers/CollimatorPhysics.hh"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"
#include "AbsBeamline/ElementBase.h"

extern Inform *gmsg;

using namespace Physics;


// Class SurfacePhysics
// ------------------------------------------------------------------------

// The attributes of class SurfacePhysics.
namespace {
    enum {
        // DESCRIPTION OF SINGLE PARTICLE:
        TYPE,       // The type of the wake
        MATERIAL,   // From of the tube
        RADIUS, // Radius of the tube
        SIGMA,
        TAU,
        SIZE
    };
}

SurfacePhysics::SurfacePhysics():
    Definition(SIZE, "SURFACEPHYSICS",
               "The \"SURFACE_PHYSICS\" statement defines data for the surface physics handler "
               "on an element."),
    handler_m(0) {
    itsAttr[TYPE] = Attributes::makeString
                    ("TYPE", "Specifies the surface physics handler: Collimator");

    itsAttr[MATERIAL] = Attributes::makeString
                        ("MATERIAL", "The material of the surface");

    itsAttr[RADIUS] = Attributes::makeReal
                      ("RADIUS", "The radius of the beam pipe [m]");

    itsAttr[SIGMA] = Attributes::makeReal
                     ("SIGMA", "Material constant dependant on the  beam pipe material");

    itsAttr[TAU] = Attributes::makeReal
                   ("TAU", "Material constant dependant on the  beam pipe material");

    SurfacePhysics *defSurfacePhysics = clone("UNNAMED_SURFACEPHYSICS");
    defSurfacePhysics->builtin = true;

    try {
        defSurfacePhysics->update();
        OpalData::getInstance()->define(defSurfacePhysics);
    } catch(...) {
        delete defSurfacePhysics;
    }
}


SurfacePhysics::SurfacePhysics(const string &name, SurfacePhysics *parent):
    Definition(name, parent),
    handler_m(parent->handler_m)
{}


SurfacePhysics::~SurfacePhysics() {
    if(handler_m)
        delete handler_m;
}


bool SurfacePhysics::canReplaceBy(Object *object) {
    // Can replace only by another SURFACEPHYSICS.
    return dynamic_cast<SurfacePhysics *>(object) != 0;
}


SurfacePhysics *SurfacePhysics::clone(const string &name) {
    return new SurfacePhysics(name, this);
}


void SurfacePhysics::execute() {
    update();
}


SurfacePhysics *SurfacePhysics::find(const string &name) {
    SurfacePhysics *sphys = dynamic_cast<SurfacePhysics *>(OpalData::getInstance()->find(name));

    if(sphys == 0) {
        throw OpalException("SurfacePhysics::find()", "SurfacePhysics \"" + name + "\" not found.");
    }
    return sphys;
}


void SurfacePhysics::update() {
    // Set default name.
    if(getOpalName().empty()) setOpalName("UNNAMED_SURFACEPHYSICS");
}


void SurfacePhysics::initSurfacePhysicsHandler(ElementBase &element, const double &major, const double &minor) {
    *gmsg << "* ************* S U R F A C E P H Y S I C S **************************************** " << endl;
    *gmsg << "* SurfacePhysics::initSurfacePhysicsHandler " << endl;
    *gmsg << "* ********************************************************************************** " << endl;

    itsElement_m = &element;
    material_m = Attributes::getString(itsAttr[MATERIAL]);

    if(Attributes::getString(itsAttr[TYPE]) == "COLLIMATOR" || Attributes::getString(itsAttr[TYPE]) == "DEGRADER") {
        Inform m("sph ");
        handler_m = new CollimatorPhysics(getOpalName(), itsElement_m, major, minor, material_m);
        m << *this << endl;
    } else {
        handler_m = 0;
        INFOMSG("no surface physics handler attached, TYPE == " << Attributes::getString(itsAttr[TYPE]) << endl);
    }

}

void SurfacePhysics::updateElement(ElementBase *element) {
    handler_m->updateElement(element);
}

void SurfacePhysics::print(std::ostream &os) const {
    os << "* ************* S U R F A C E P H Y S I C S **************************************** " << std::endl;
    os << "* SURFACEPHYSICS " << getOpalName() << '\n'
       << "* MATERIAL       " << Attributes::getString(itsAttr[MATERIAL]) << '\n'
       << "* RADIUS         " << Attributes::getReal(itsAttr[RADIUS]) << '\n'
       << "* SIGMA          " << Attributes::getReal(itsAttr[SIGMA]) << '\n'
       << "* TAU            " << Attributes::getReal(itsAttr[TAU]) << '\n';
    os << "* ********************************************************************************** " << std::endl;
}
