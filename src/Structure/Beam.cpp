// ------------------------------------------------------------------------
// $RCSfile: Beam.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Beam
//   The class for the OPAL BEAM command.
//
//   Note by JMJ 4/4/2000:
//   At present it is possible to define unnormalised and normalised
//   emittances independently.  Similarly for GAMMA, ENERGY and PC.
//   We need to specify a basic set of independent class data members
//   and then define methods for setting the basic ones via their values.
//   The OPAL9 parser can probably be relied upon to cause later definitions
//   on a command line to override earlier ones.
//   The basic set should be, e.g.,
//      {PARTICLE,MASS,CHARGE,PC,EX,EY,ET,
//          KBUNCH,NPART,BUNCHED,RADIATE,DAMP,QUANTUM}
//   This will need some care in transforming things.   The appropriate
//   setXXX methods would just involve solving an equations defining the
//   XXX quantity in terms of the base ones and deciding which one it is
//   to override.  E.g., setting EXN would reset EX.
//   Maybe print a warning message when this happens.
//   It would, of course, be preferable not to allow the user to specify
//   GAMMA when he can specify PC, or EXN when he can specify EX but
//   some people will want these conveniences.
// ------------------------------------------------------------------------
//
// $Date: 2003/08/11 22:09:00 $
// $Author: dbruhwil $
//
// ------------------------------------------------------------------------

#include "Structure/Beam.h"
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Expressions/SAutomatic.h"
#include "Expressions/SRefExpr.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"

#include <cmath>
#include <iterator>

using namespace Expressions;

// Class Beam
// ------------------------------------------------------------------------

// The attributes of class Beam.
namespace {
    enum {
        // DESCRIPTION OF SINGLE PARTICLE:
        PARTICLE,   // The particle name
        MASS,       // The particle rest mass in GeV
        CHARGE,     // The particle charge in proton charges
        ENERGY,     // The particle energy in GeV
        PC,         // The particle momentum in GeV/c
        GAMMA,      // ENERGY / MASS

        // BEAM CURRENT AND EMITTANCES:
        BCURRENT,   // Beam current in A
        EX,         // Horizontal emittance
        EXN,        // Normalised horizontal emittance
        EY,         // Vertical emittance
        EYN,        // Normalises vertical emittance
        ET,         // Longitudinal emittance

        // BEAM FREQUENCY
        BFREQ,  // Beam frequency in MHz

        // DESCRIPTION OF BUNCHES:
        KBUNCH,     // Number of bunches
        NPART,      // Number of particles per bunch
        NSLICE,     // Number of slices per bunch
        BUNCHED,    // True, if the beam is bunched
        RADIATE,    // True, if the particles radiate energy
        DAMP,       // True, if damping is considered
        QUANTUM,    // True, if quantum excitation is considered
        SPACECHARGE,// True, if space charge computation is desired

        // DESCIPTION OF FIELDSOLVER will be eventually a separate opject
        FIELDSOLVER,
        // THE INTEGRATION TIMESTEP IN SEC
        SIZE
    };
}


const double Beam::energy_scale = 1.0e9;


Beam::Beam():
    Definition(SIZE, "BEAM",
               "The \"BEAM\" statement defines data for a the particles "
               "in a beam."),
    reference(1.0, Physics::m_p *energy_scale, 1.0 * energy_scale) {

    itsAttr[FIELDSOLVER] = Attributes::makeString
                           ("FIELDSOLVER", "Name of the attached field solver ");

    // DESCRIPTION OF SINGLE PARTICLE:
    itsAttr[PARTICLE] = Attributes::makeString
                        ("PARTICLE", "Name of particle to be used (Ch 9 of the user manual) ");
    itsAttr[MASS] = Attributes::makeReal
                    ("MASS", "Particle rest mass in GeV");
    itsAttr[CHARGE] = Attributes::makeReal
                      ("CHARGE", "Particle charge in proton charges");
    itsAttr[ENERGY] = Attributes::makeReal
                      ("ENERGY", "Particle energy in GeV");
    itsAttr[PC] = Attributes::makeReal
                  ("PC", "Particle momentum in GeV/c");
    PtrToScalar<double> expr = new SRefExpr<double>("P0", "");
    itsAttr[PC].set(new SAutomatic<double>(expr));
    itsAttr[GAMMA] = Attributes::makeReal
                     ("GAMMA", "ENERGY / MASS");

    // BEAM CURRENT AND EMITTANCES:
    itsAttr[BCURRENT] = Attributes::makeReal
                        ("BCURRENT", "Beam current in A (all bunches)");
    itsAttr[EX] = Attributes::makeReal
                  ("EX", "Horizontal emittance");
    itsAttr[EXN] = Attributes::makeReal
                   ("EXN", "Normalised horizontal emittance");
    itsAttr[EY] = Attributes::makeReal
                  ("EY", "Vertical emittance");
    itsAttr[EYN] = Attributes::makeReal
                   ("EYN", "Normalised vertical emittance");
    itsAttr[ET] = Attributes::makeReal
                  ("ET", "Longitudinal emittance");

    // BEAM FREQUENCY
    itsAttr[BFREQ] = Attributes::makeReal
                     ("BFREQ", "Beam frequency in MHz (all bunches)");

    // DESCRIPTION OF BUNCHES:
    itsAttr[KBUNCH] = Attributes::makeReal
                      ("KBUNCH", "Number of bunches in beam", 1.0);
    itsAttr[NPART] = Attributes::makeReal
                     ("NPART", "Number of particles in bunch");
    itsAttr[NSLICE] = Attributes::makeReal
                      ("NSLICE", "Number of slices in bunch");
    itsAttr[BUNCHED] = Attributes::makeBool
                       ("BUNCHED", "True, if the beam is bunched");
    itsAttr[RADIATE] = Attributes::makeBool
                       ("RADIATE", "True, if the particles radiate energy");
    itsAttr[DAMP] = Attributes::makeBool
                    ("DAMP", "True, if damping is considered");
    itsAttr[QUANTUM] = Attributes::makeBool
                       ("QUANTUM", "True, if quantum excitation is considered");
    itsAttr[SPACECHARGE] = Attributes::makeBool
                           ("SPACECHARGE", "True, if space charge computation is desired");

    // Set up default beam.
    Beam *defBeam = clone("UNNAMED_BEAM");
    defBeam->builtin = true;

    try {
        defBeam->update();
        OpalData::getInstance()->define(defBeam);
    } catch(...) {
        delete defBeam;
    }

    registerOwnership(AttributeHandler::STATEMENT);
}


Beam::Beam(const std::string &name, Beam *parent):
    Definition(name, parent),
    reference(parent->reference)
{}


Beam::~Beam()
{}


bool Beam::canReplaceBy(Object *object) {
    // Can replace only by another BEAM.
    return dynamic_cast<Beam *>(object) != 0;
}


Beam *Beam::clone(const std::string &name) {
    return new Beam(name, this);
}


void Beam::execute() {
    update();
}


Beam *Beam::find(const std::string &name) {
    Beam *beam = dynamic_cast<Beam *>(OpalData::getInstance()->find(name));

    if(beam == 0) {
        throw OpalException("Beam::find()", "Beam \"" + name + "\" not found.");
    }

    return beam;
}

size_t Beam::getNumberOfParticles() {
    return (size_t)Attributes::getReal(itsAttr[NPART]);
}

size_t Beam::getNumberOfSlices() {
    return (size_t)Attributes::getReal(itsAttr[NSLICE]);
}

double Beam::getEX() const {
    return Attributes::getReal(itsAttr[EX]);
}


double Beam::getEY() const {
    return Attributes::getReal(itsAttr[EY]);
}


double Beam::getET() const {
    return Attributes::getReal(itsAttr[ET]);
}


const PartData &Beam::getReference() const {
    // Cast away const, to allow logically constant Beam to update.
    const_cast<Beam *>(this)->update();
    return reference;
}

double Beam::getCurrent() const {
    return Attributes::getReal(itsAttr[BCURRENT]);
}

double Beam::getCharge() const {
    return Attributes::getReal(itsAttr[CHARGE]);
}

double Beam::getMass() const {
    return Attributes::getReal(itsAttr[MASS]);
}

std::string Beam::getParticleName() const {
    return Attributes::getString(itsAttr[PARTICLE]);
}

double Beam::getFrequency() const {
    return Attributes::getReal(itsAttr[BFREQ]);
}

void Beam::setEX(double value) {
    Attributes::setReal(itsAttr[EX], value);
}


void Beam::setEY(double value) {
    Attributes::setReal(itsAttr[EY], value);
}


void Beam::setET(double value) {
    Attributes::setReal(itsAttr[ET], value);
}


void Beam::update() {
    // Find the particle name.
    if(itsAttr[PARTICLE]) {
        static const char *names[] = {
            "ELECTRON", "PROTON", "POSITRON", "ANTIPROTON", "CARBON", "HMINUS", "URANIUM", "MUON", "DEUTERON", "XENON", "H2P"
        };

        static const double masses[] = {
            Physics::m_e,
            Physics::m_p,
            Physics::m_e,
            Physics::m_p,
            Physics::m_c,
            Physics::m_hm,
            Physics::m_u,
            Physics::m_mu,
            Physics::m_d,
            Physics::m_xe,
            2 * Physics::m_p
        };

        static const double charges[] = {
            -1.0, 1.0, 1.0, -1.0, 12.0, -1.0, 35.0, -1.0, 1.0, 20.0, 1.0
        };
        const unsigned int numParticleNames = std::end(names) - std::begin(names);

        std::string pName  = Attributes::getString(itsAttr[PARTICLE]);
        for(unsigned int i = 0; i < numParticleNames; ++ i) {
            if(pName == names[i]) {
                Attributes::setReal(itsAttr[MASS], masses[i]);
                Attributes::setReal(itsAttr[CHARGE], charges[i]);
                break;
            }
        }
    }

    // Set up particle reference; convert all to eV for CLASSIC.
    double mass =
        (itsAttr[MASS] ? Attributes::getReal(itsAttr[MASS]) : Physics::m_p) * energy_scale;
    double charge = itsAttr[CHARGE] ? Attributes::getReal(itsAttr[CHARGE]) : 1.0;
    reference = PartData(charge, mass, 1.0);

    if(itsAttr[GAMMA]) {
        double gamma = Attributes::getReal(itsAttr[GAMMA]);
        if(gamma > 1.0) {
            reference.setGamma(gamma);
        } else {
            throw OpalException("Beam::execute()",
                                "\"GAMMA\" should be greater than 1.");
        }
    } else if(itsAttr[ENERGY]) {
        double energy = Attributes::getReal(itsAttr[ENERGY]) * energy_scale;
        if(energy > reference.getM()) {
            reference.setE(energy);
        } else {
            throw OpalException("Beam::execute()",
                                "\"ENERGY\" should be greater than \"MASS\".");
        }
    } else if(itsAttr[PC]) {
        double pc = Attributes::getReal(itsAttr[PC]) * energy_scale;
        if(pc > 0.0) {
            reference.setP(pc);
        } else {
            throw OpalException("Beam::execute()",
                                "\"PC\" should be greater than 0.");
        }
    };

    // Emittances added above by JMJ 4/4/2000 so that UNNAMED_BEAM has
    // proper defaults (problem found by Julien Pancin, see email today).
    // This did not seem to work for him to deleted them again.  Confused.
    // Slightly worried about this preventing anyone giving
    // normalised emittances EXN, EYN, ETN.

    // MISSING: Other data for BEAM.

    // Set default name.
    if(getOpalName().empty()) setOpalName("UNNAMED_BEAM");
}


//ff
double Beam::getGamma() const { //obtain value for gamma
    return Attributes::getReal(itsAttr[GAMMA]);
}

//ff
double Beam::getPC() const { //obtain value for PC
    return Attributes::getReal(itsAttr[PC]);
}


void Beam::print(std::ostream &os) const {
    double charge = Attributes::getReal(itsAttr[CHARGE]);
    os << "* ************* B E A M ************************************************************ " << std::endl;
    os << "* BEAM        " << getOpalName() << '\n'
       << "* PARTICLE    " << Attributes::getString(itsAttr[PARTICLE]) << '\n'
       << "* CURRENT     " << Attributes::getReal(itsAttr[BCURRENT]) << " A\n"
       << "* FREQUENCY   " << Attributes::getReal(itsAttr[BFREQ]) << " MHz\n"
       << "* CHARGE      " << (charge > 0 ? '+' : '-') << "e * " << std::abs(charge) << " \n"
       << "* REST MASS   " << Attributes::getReal(itsAttr[MASS]) << " GeV\n"
       << "* MOMENTUM    " << Attributes::getReal(itsAttr[PC])   << '\n'
       << "* NPART       " << Attributes::getReal(itsAttr[NPART])   << '\n';
    os << "* ********************************************************************************** " << std::endl;
}
