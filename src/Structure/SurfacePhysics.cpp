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
	ENABLERUTHERFORD,
	LOWENERGYTHR,
        SIGMA,
        TAU,
	NPART,
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

    itsAttr[NPART] = Attributes::makeReal("NPART", "Number of particles in bunch");

    itsAttr[ENABLERUTHERFORD] = Attributes::makeBool("ENABLERUTHERFORD", "Enable large angle scattering",true);

    itsAttr[LOWENERGYTHR] = Attributes::makeReal("LOWENERGYTHR", "Lower Energy threshold Bethe-Block, manual Ch. 18 ",1E-4);

    SurfacePhysics *defSurfacePhysics = clone("UNNAMED_SURFACEPHYSICS");
    defSurfacePhysics->builtin = true;

    try {
        defSurfacePhysics->update();
        OpalData::getInstance()->define(defSurfacePhysics);
    } catch(...) {
        delete defSurfacePhysics;
    }

    registerOwnership(AttributeHandler::STATEMENT);
}


SurfacePhysics::SurfacePhysics(const std::string &name, SurfacePhysics *parent):
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


SurfacePhysics *SurfacePhysics::clone(const std::string &name) {
    return new SurfacePhysics(name, this);
}


void SurfacePhysics::execute() {
    update();
}


SurfacePhysics *SurfacePhysics::find(const std::string &name) {
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


void SurfacePhysics::initSurfacePhysicsHandler(ElementBase &element) {
    *gmsg << "* ************* S U R F A C E P H Y S I C S **************************************** " << endl;
    *gmsg << "* SurfacePhysics::initSurfacePhysicsHandler " << endl;
    *gmsg << "* ********************************************************************************** " << endl;

    itsElement_m = &element;
    material_m = Attributes::getString(itsAttr[MATERIAL]);
    enableRutherfordScattering_m = Attributes::getBool(itsAttr[ENABLERUTHERFORD]);
    lowEnergyThr_m = Attributes::getReal(itsAttr[LOWENERGYTHR]);

    if(Attributes::getString(itsAttr[TYPE]) == "CCOLLIMATOR" || Attributes::getString(itsAttr[TYPE]) == "COLLIMATOR" || Attributes::getString(itsAttr[TYPE]) == "DEGRADER") {
      handler_m = new CollimatorPhysics(getOpalName(), itsElement_m, material_m, enableRutherfordScattering_m, lowEnergyThr_m);
        *gmsg << *this << endl;
    } else {
        handler_m = 0;
        INFOMSG("no surface physics handler attached, TYPE == " << Attributes::getString(itsAttr[TYPE]) << endl);
    }

    if (Ippl::getNodes()>1) {
      *gmsg << "\033[01;35m";
      *gmsg << "* ************* A D V I S O R Y  *************************************************** " << endl;
      *gmsg << "* It is known that surface physics in combination with a parallel run can have       " << endl;
      *gmsg << "* side-effects such as segmentation faults or the simulation just stops. In such a   " << endl;
      *gmsg << "* case, either play with the number of simulation particles or use 1 core only.      " << endl;
      *gmsg << "* We are working on the issue, progress on that matter can be obtained at            " << endl;
      *gmsg << "* https://gitlab.psi.ch/OPAL/src/issues/137.                                         " << endl;
      *gmsg << "* ********************************************************************************** " << endl;
      *gmsg << "\u001b[0m" << endl;
    }
}

void SurfacePhysics::updateElement(ElementBase *element) {
    handler_m->updateElement(element);
}

void SurfacePhysics::print(std::ostream &os) const {
    os << "* ************* S U R F A C E P H Y S I C S **************************************** " << std::endl;
    os << "* SURFACEPHYSICS  \t " << getOpalName() << '\n'
       << "* MATERIAL        \t " << Attributes::getString(itsAttr[MATERIAL]) << '\n'
       << "* RADIUS          \t " << Attributes::getReal(itsAttr[RADIUS]) << '\n'
       << "* SIGMA           \t " << Attributes::getReal(itsAttr[SIGMA]) << '\n'
       << "* TAU             \t " << Attributes::getReal(itsAttr[TAU]) << '\n'
       << "* LOWENERGYTHR    \t " << Attributes::getReal(itsAttr[LOWENERGYTHR]) << " (GeV) \n";

    if (Attributes::getBool(itsAttr[ENABLERUTHERFORD]))
      os << "* RUTHERFORD SCAT \t ENABLED" << '\n';
    else
      os << "* RUTHERFORD SCAT \t DISABLED" << '\n';
    os << "* ********************************************************************************** " << std::endl;
}
