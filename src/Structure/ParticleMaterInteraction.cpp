// ------------------------------------------------------------------------
// $RCSfile: ParticleMaterInteraction.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ParticleMaterInteraction
//   The class for the OPAL PARTICLEMATERINTERACTION command.
//
// $Date: 2009/07/14 22:09:00 $
// $Author: C. Kraus $
//
// ------------------------------------------------------------------------

#include "Structure/ParticleMaterInteraction.h"
#include "Solvers/CollimatorPhysics.hh"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"
#include "AbsBeamline/ElementBase.h"
#include "Utilities/Util.h"

extern Inform *gmsg;

using namespace Physics;


// Class ParticleMaterInteraction
// ------------------------------------------------------------------------

// The attributes of class ParticleMaterInteraction.
namespace {
    enum {
        // DESCRIPTION OF SINGLE PARTICLE:
        TYPE,       // The type of the wake
        MATERIAL,   // From of the tube
        RADIUS, // Radius of the tube
        SIGMA,
        TAU,
	NPART,
        SIZE
    };
}

ParticleMaterInteraction::ParticleMaterInteraction():
    Definition(SIZE, "PARTICLEMATERINTERACTION",
               "The \"SURFACE_PHYSICS\" statement defines data for the particle mater interaction handler "
               "on an element."),
    handler_m(0) {
    itsAttr[TYPE] = Attributes::makeString
                    ("TYPE", "Specifies the particle mater interaction handler: Collimator");

    itsAttr[MATERIAL] = Attributes::makeString
                        ("MATERIAL", "The material of the surface");

    itsAttr[RADIUS] = Attributes::makeReal
                      ("RADIUS", "The radius of the beam pipe [m]");

    itsAttr[SIGMA] = Attributes::makeReal
                     ("SIGMA", "Material constant dependant on the  beam pipe material");

    itsAttr[TAU] = Attributes::makeReal
                   ("TAU", "Material constant dependant on the  beam pipe material");

    itsAttr[NPART] = Attributes::makeReal("NPART", "Number of particles in bunch");

    ParticleMaterInteraction *defParticleMaterInteraction = clone("UNNAMED_PARTICLEMATERINTERACTION");
    defParticleMaterInteraction->builtin = true;

    try {
        defParticleMaterInteraction->update();
        OpalData::getInstance()->define(defParticleMaterInteraction);
    } catch(...) {
        delete defParticleMaterInteraction;
    }

    registerOwnership(AttributeHandler::STATEMENT);
}


ParticleMaterInteraction::ParticleMaterInteraction(const std::string &name, ParticleMaterInteraction *parent):
    Definition(name, parent),
    handler_m(parent->handler_m)
{}


ParticleMaterInteraction::~ParticleMaterInteraction() {
    if(handler_m)
        delete handler_m;
}


bool ParticleMaterInteraction::canReplaceBy(Object *object) {
    // Can replace only by another PARTICLEMATERINTERACTION.
    return dynamic_cast<ParticleMaterInteraction *>(object) != 0;
}


ParticleMaterInteraction *ParticleMaterInteraction::clone(const std::string &name) {
    return new ParticleMaterInteraction(name, this);
}


void ParticleMaterInteraction::execute() {
    update();
}


ParticleMaterInteraction *ParticleMaterInteraction::find(const std::string &name) {
    ParticleMaterInteraction *parmatint = dynamic_cast<ParticleMaterInteraction *>(OpalData::getInstance()->find(name));

    if(parmatint == 0) {
        throw OpalException("ParticleMaterInteraction::find()", "ParticleMaterInteraction \"" + name + "\" not found.");
    }
    return parmatint;
}


void ParticleMaterInteraction::update() {
    // Set default name.
    if(getOpalName().empty()) setOpalName("UNNAMED_PARTICLEMATERINTERACTION");
}


void ParticleMaterInteraction::initParticleMaterInteractionHandler(ElementBase &element) {
    *gmsg << "* ************* P A R T I C L E  M A T E R  I N T E R A C T I O N ****************** " << endl;
    *gmsg << "* ParticleMaterInteraction::initParticleMaterInteractionHandler " << endl;
    *gmsg << "* ********************************************************************************** " << endl;

    itsElement_m = &element;
    material_m = Util::toUpper(Attributes::getString(itsAttr[MATERIAL]));

    const std::string type = Util::toUpper(Attributes::getString(itsAttr[TYPE]));
    if(type == "CCOLLIMATOR" ||
       type == "COLLIMATOR" ||
       type == "DEGRADER") {

        handler_m = new CollimatorPhysics(getOpalName(), itsElement_m, material_m);
        *gmsg << *this << endl;
    } else {
        handler_m = 0;
        INFOMSG(getOpalName() + ": no particle mater interaction handler attached, TYPE == " << Attributes::getString(itsAttr[TYPE]) << endl);
    }

}

void ParticleMaterInteraction::updateElement(ElementBase *element) {
    handler_m->updateElement(element);
}

void ParticleMaterInteraction::print(std::ostream &os) const {
    os << "* ************* P A R T I C L E  M A T E R  I N T E R A C T I O N ****************** " << std::endl;
    os << "* PARTICLEMATERINTERACTION " << getOpalName() << '\n'
       << "* MATERIAL       " << Attributes::getString(itsAttr[MATERIAL]) << '\n';
    os << "* ********************************************************************************** " << std::endl;
}
