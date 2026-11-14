//
// Class Beam
//   The class for the OPAL BEAM command.
//   A BEAM definition is used by most physics commands to define the
//   particle charge and the reference momentum, together with some other data.
//
// Copyright (c) 200x - 2021, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#include "Structure/Beam.h"

#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/OpalData.h"
#include "Algorithms/PartData.h"
#include "Attributes/Attributes.h"
#include "Physics/ParticleProperties.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"
#include "Utilities/OpalException.h"

#include <cmath>

using namespace Expressions;

// The attributes of class Beam.
namespace {
    constexpr const char* photonParticleName = "PHOTON";

    enum {
        PARTICLE,         // The particle name
        MASS,             // The particle rest mass in GeV
        CHARGE,           // The particle charge in proton charges
        ENERGY,           // The particle energy in GeV
        PC,               // The particle momentum in GeV/c
        GAMMA,            // ENERGY / MASS
        BCURRENT,         // Legacy, unused in OPALX (holdover from OPALCycl)
        BFREQ,            // Legacy, unused in OPALX (holdover from OPALCycl)
        BCHARGE,          // Bunch charge in C
        NALLOC,           // Allocation size (macroparticles) for this beam
        SOURCES,          // Name of EMISSIONSOURCELIST
        GLOBALPROCESSES,  // Global physics processes active for this beam
        DAUGHTERBEAM,     // Name of the beam that receives decay daughter particles
        POLARIZATION,     // Initial polarization vector P (rest-frame, lab-frame axes)
        TUNEINITIAL,      // Spectral launch triples, independent of bunch distribution.
        ZSTOP,            // Path length at which this beam stops being tracked
        UNTILMAXSTEPS,    // Keep the run going for this beam until MAXSTEPS
        REFDIRECTION,     // Direction of the reference momentum for file distributions
        SIZE
    };
}  // namespace

Beam::Beam()
    : Definition(
              SIZE, "BEAM",
              "The \"BEAM\" statement defines data for the particles "
              "in a beam."),
      reference(1.0, Physics::m_p * Units::GeV2eV, 1.0 * Units::GeV2eV) {
    itsAttr[PARTICLE] = Attributes::makePredefinedString(
            "PARTICLE", "Name of particle to be used",
            {"PHOTON", "ELECTRON", "POSITRON", "MUON", "PION", "PROTON", "ANTIPROTON", "DEUTERON",
             "HMINUS", "H2P", "ALPHA", "CARBON", "XENON", "URANIUM"});

    itsAttr[MASS]        = Attributes::makeReal("MASS", "Particle rest mass [GeV]");
    itsAttr[TUNEINITIAL] = Attributes::makeRealArray(
            "TUNEINITIAL",
            "Repeated spectral launches: kinetic energy [GeV], radius [m], radial momentum [mc].");

    itsAttr[CHARGE] = Attributes::makeReal("CHARGE", "Particle charge in proton charges");

    itsAttr[ENERGY] = Attributes::makeReal(
            "ENERGY", "Particle energy [GeV]. Not allowed on another beam's DAUGHTERBEAM.");

    itsAttr[PC] = Attributes::makeReal(
            "PC", "Particle momentum [GeV/c]. Not allowed on another beam's DAUGHTERBEAM.");

    itsAttr[GAMMA] = Attributes::makeReal(
            "GAMMA", "ENERGY / MASS. Not allowed on another beam's DAUGHTERBEAM.");

    itsAttr[BCURRENT] =
            Attributes::makeReal("BCURRENT", "Legacy, unused in OPALX. Use BCHARGE instead.");

    itsAttr[BFREQ] = Attributes::makeReal("BFREQ", "Legacy, unused in OPALX. Use BCHARGE instead.");

    itsAttr[BCHARGE] = Attributes::makeReal("BCHARGE", "Bunch charge [C]");

    itsAttr[NALLOC] =
            Attributes::makeReal("NALLOC", "Allocation size (macroparticles) for this beam");

    itsAttr[SOURCES] = Attributes::makeString(
            "SOURCES",
            "Name of the emission sources list (EMISSIONSOURCELIST). Required, except on "
            "another beam's DAUGHTERBEAM, where it is not allowed.");

    itsAttr[GLOBALPROCESSES] = Attributes::makeUpperCaseStringArray(
            "GLOBALPROCESSES", "Global physics processes active for this beam.");

    itsAttr[DAUGHTERBEAM] = Attributes::makeString(
            "DAUGHTERBEAM",
            "Name of the BEAM that receives this beam's decay daughters (needs "
            "GLOBALPROCESSES = {\"DECAY\"}). That BEAM gets particles only from the decay and "
            "uses this beam's reference particle, so it must not set GAMMA, ENERGY, PC or "
            "SOURCES.");

    itsAttr[POLARIZATION] = Attributes::makeRealArray(
            "POLARIZATION",
            "Initial polarization vector P = {Px, Py, Pz} for this beam "
            "(rest-frame Pauli expectation values along lab-frame axes). "
            "Must be length-3 with |P| in [0, 1]. Only valid for PARTICLE=MUON; "
            "setting POLARIZATION on any other species is rejected. Setting "
            "POLARIZATION (even {0,0,0}) enables per-particle spin tracking for this "
            "beam (Thomas-BMT integration); leaving it unset disables spin tracking. "
            "For muons produced by an upstream decay (e.g. PionDecay) the per-particle "
            "value is overwritten by the decay, but POLARIZATION must still be set to "
            "enable the spin storage that receives it.");

    itsAttr[ZSTOP] = Attributes::makeReal(
            "ZSTOP",
            "Path length [m] at which this beam stops: once its reference reaches it, the beam "
            "is tracked until its last particle passes the plane through the reference normal to "
            "its momentum. Only for beams listed after the first in TRACK BEAMS that start at the "
            "start of the line; the first beam uses TRACK ZSTOP.");

    itsAttr[UNTILMAXSTEPS] = Attributes::makeBool(
            "UNTILMAXSTEPS",
            "If true, the run continues for this beam until MAXSTEPS is used up, after the beams "
            "with a ZSTOP are done. Not together with ZSTOP.",
            false);

    itsAttr[REFDIRECTION] = Attributes::makeRealArray(
            "REFDIRECTION",
            "Direction {dx, dy, dz} of the reference momentum in the beamline frame; OPALX "
            "normalises it. Only for FROMFILE and EMITTEDFROMFILE distributions, together with "
            "PC, ENERGY or GAMMA, which give its length. Default {0, 0, 1}.");

    registerOwnership(AttributeHandler::STATEMENT);
}

void Beam::defineBuiltins() {
    Beam* defBeam = new Beam();
    defBeam->setOpalName("UNNAMED_BEAM");
    defBeam->builtin = true;

    try {
        defBeam->update();
        OpalData::getInstance()->define(defBeam);
    } catch (...) {
        delete defBeam;
    }
}

Beam::~Beam() {}

bool Beam::canReplaceBy(Object* object) {
    // Can replace only by another BEAM.
    return dynamic_cast<Beam*>(object) != 0;
}

void Beam::execute() {
    const bool photon = itsAttr[PARTICLE] && getParticleName() == photonParticleName;

    if (photon) {
        if (!itsAttr[ENERGY]) {
            throw OpalException("Beam::execute()", "\"ENERGY\" must be set for PARTICLE=PHOTON.");
        }
        if (itsAttr[MASS]) {
            throw OpalException(
                    "Beam::execute()",
                    "\"MASS\" is not allowed for PARTICLE=PHOTON. Use \"ENERGY\".");
        }
        if (itsAttr[CHARGE]) {
            throw OpalException(
                    "Beam::execute()", "\"CHARGE\" is not allowed for PARTICLE=PHOTON.");
        }
        if (itsAttr[PC]) {
            throw OpalException(
                    "Beam::execute()",
                    "\"PC\" is not allowed for PARTICLE=PHOTON. Use \"ENERGY\".");
        }
        if (itsAttr[GAMMA]) {
            throw OpalException(
                    "Beam::execute()",
                    "\"GAMMA\" is not allowed for PARTICLE=PHOTON. Use \"ENERGY\".");
        }
        if (itsAttr[SOURCES]) {
            throw OpalException(
                    "Beam::execute()", "\"SOURCES\" is not allowed for PARTICLE=PHOTON.");
        }
    }

    if (itsAttr[BCURRENT] || itsAttr[BFREQ]) {
        throw OpalException(
                "Beam::execute()",
                "\"BCURRENT\" and \"BFREQ\" are no longer used in OPALX. "
                "Use \"BCHARGE\" [C] to specify the bunch charge directly.");
    }

    update();

    if (!(itsAttr[PARTICLE]) && (!itsAttr[MASS] || !(itsAttr[CHARGE]))) {
        throw OpalException(
                "Beam::execute()",
                "The beam particle hasn't been set. "
                "Set either \"PARTICLE\" or \"MASS\" and \"CHARGE\".");
    }

    if (!(itsAttr[NALLOC])) {
        throw OpalException("Beam::execute()", "\"NALLOC\" must be set.");
    }

    if (photon) {
        const double energy = Attributes::getReal(itsAttr[ENERGY]);
        if (energy <= 0.0) {
            throw OpalException(
                    "Beam::execute()", "\"ENERGY\" should be greater than 0 for PARTICLE=PHOTON.");
        }
        return;
    }

    // SOURCES is checked by TrackRun: only there is it known whether this beam is another
    // beam's DAUGHTERBEAM (which must not set SOURCES).

    // Currently supported global process names (extend as new processes are implemented).
    for (const std::string& name : getGlobalProcessNames()) {
        if (name != "DECAY") {
            throw OpalException(
                    "Beam::execute()", "Unsupported entry in \"GLOBALPROCESSES\": \"" + name
                                               + "\". Supported values: DECAY.");
        }
    }

    validatePolarization();
    validateReferenceDirection();
}

void Beam::validateReferenceDirection() const {
    if (!hasReferenceDirection()) {
        return;
    }
    const std::vector<double> dir = Attributes::getRealArray(itsAttr[REFDIRECTION]);
    if (dir.size() != 3) {
        throw OpalException(
                "Beam::execute()", "\"REFDIRECTION\" must be a length-3 vector {dx, dy, dz}.");
    }
    if (!(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2] > 0.0)) {
        throw OpalException("Beam::execute()", "\"REFDIRECTION\" must not be zero.");
    }
    if (!hasExplicitEnergy()) {
        throw OpalException(
                "Beam::execute()",
                "\"REFDIRECTION\" needs \"PC\", \"ENERGY\" or \"GAMMA\" for the length of "
                "the reference momentum.");
    }
}

void Beam::validatePolarization() const {
    if (!itsAttr[POLARIZATION]) {
        return;
    }
    const std::vector<double> pol = Attributes::getRealArray(itsAttr[POLARIZATION]);
    if (pol.empty()) {
        // Attribute exists in the registry but the user did not specify a value;
        // treat as "not set" and skip all polarization-specific validation.
        return;
    }
    if (pol.size() != 3) {
        throw OpalException(
                "Beam::execute()", "\"POLARIZATION\" must be a length-3 vector {Px, Py, Pz}.");
    }
    // Only muons are currently spin-tracked. Reject POLARIZATION on any other
    // species regardless of magnitude — setting it (even to zero) is meaningless
    // for spin-0 particles and for species whose g-factor isn't wired into T-BMT.
    if (!itsAttr[PARTICLE]) {
        throw OpalException(
                "Beam::execute()",
                "\"POLARIZATION\" can only be set when \"PARTICLE\" is also set.");
    }
    const ParticleType pType = ParticleProperties::getParticleType(getParticleName());
    if (pType != ParticleType::MUON) {
        throw OpalException(
                "Beam::execute()", "\"POLARIZATION\" is only supported for PARTICLE=MUON; got \""
                                           + getParticleName() + "\".");
    }
    const double pMag2 = pol[0] * pol[0] + pol[1] * pol[1] + pol[2] * pol[2];
    if (pMag2 > 1.0 + 1.0e-12) {
        throw OpalException(
                "Beam::execute()", "\"POLARIZATION\" magnitude must be in [0, 1]; got |P| = "
                                           + std::to_string(std::sqrt(pMag2)) + ".");
    }
}

std::string Beam::getEmissionSourceListName() const {
    if (!itsAttr[SOURCES]) {
        throw OpalException(
                "Beam::getEmissionSourceListName()",
                "\"SOURCES\" must be set for a beam (name of EMISSIONSOURCELIST).");
    }

    const std::string name = Attributes::getString(itsAttr[SOURCES]);
    if (name.empty()) {
        throw OpalException(
                "Beam::getEmissionSourceListName()",
                "\"SOURCES\" must not be empty for a beam (name of EMISSIONSOURCELIST).");
    }
    return name;
}

std::vector<std::string> Beam::getGlobalProcessNames() const {
    return Attributes::getStringArray(itsAttr[GLOBALPROCESSES]);
}

std::string Beam::getDaughterBeamName() const {
    return Attributes::getString(itsAttr[DAUGHTERBEAM]);
}

std::vector<double> Beam::getPolarization() const {
    std::vector<double> pol = Attributes::getRealArray(itsAttr[POLARIZATION]);
    if (pol.empty()) {
        return {0.0, 0.0, 0.0};
    }
    return pol;
}

std::optional<double> Beam::getZStop() const {
    if (!itsAttr[ZSTOP]) return std::nullopt;
    return Attributes::getReal(itsAttr[ZSTOP]);
}

bool Beam::getUntilMaxSteps() const { return Attributes::getBool(itsAttr[UNTILMAXSTEPS]); }

bool Beam::hasPolarization() const {
    return !Attributes::getRealArray(itsAttr[POLARIZATION]).empty();
}

std::vector<double> Beam::getReferenceDirection() const {
    if (!hasReferenceDirection()) {
        return {0.0, 0.0, 1.0};
    }
    std::vector<double> dir = Attributes::getRealArray(itsAttr[REFDIRECTION]);
    const double norm       = std::sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
    for (double& d : dir) {
        d /= norm;
    }
    return dir;
}

bool Beam::hasReferenceDirection() const {
    return !Attributes::getRealArray(itsAttr[REFDIRECTION]).empty();
}

Beam* Beam::find(const std::string& name) {
    Beam* beam = dynamic_cast<Beam*>(OpalData::getInstance()->find(name));

    if (beam == 0) {
        throw OpalException("Beam::find()", "Beam \"" + name + "\" not found.");
    }

    return beam;
}

size_t Beam::getNumAlloc() const {
    if (Attributes::getReal(itsAttr[NALLOC]) > 0) {
        return (size_t)Attributes::getReal(itsAttr[NALLOC]);
    } else {
        throw OpalException(
                "Beam::getNumAlloc()",
                "Wrong allocation size for beam! \"NALLOC\" must be positive");
    }
}

const PartData& Beam::getReference() const {
    // Cast away const, to allow logically constant Beam to update.
    const_cast<Beam*>(this)->update();
    return reference;
}

double Beam::getCurrent() const { return Attributes::getReal(itsAttr[BCURRENT]); }

double Beam::getBunchCharge() const { return Attributes::getReal(itsAttr[BCHARGE]); }

double Beam::getCharge() const { return Attributes::getReal(itsAttr[CHARGE]); }

double Beam::getMass() const { return Attributes::getReal(itsAttr[MASS]); }

double Beam::getMomentum() const { return reference.getP() * Units::eV2GeV; }

std::string Beam::getParticleName() const { return Attributes::getString(itsAttr[PARTICLE]); }

bool Beam::isPhoton() const { return itsAttr[PARTICLE] && getParticleName() == photonParticleName; }

double Beam::getFrequency() const { return Attributes::getReal(itsAttr[BFREQ]); }

bool Beam::hasExplicitEnergy() const { return itsAttr[GAMMA] || itsAttr[ENERGY] || itsAttr[PC]; }

bool Beam::hasEmissionSources() const { return itsAttr[SOURCES]; }

std::vector<double> Beam::getTuneInitial() const {
    return Attributes::getRealArray(itsAttr[TUNEINITIAL]);
}

double Beam::getChargePerParticle() const {
    return std::copysign(1.0, getCharge()) * getBunchCharge() / getNumAlloc();
}

double Beam::getMassPerParticle() const {
    return getMass() * getChargePerParticle() / (getCharge() * Physics::q_e);
}

void Beam::update() {
    if (itsAttr[PARTICLE]) {
        std::string pName  = getParticleName();
        ParticleType pType = ParticleProperties::getParticleType(pName);
        if (!itsAttr[MASS]) {
            Attributes::setReal(itsAttr[MASS], ParticleProperties::getParticleMass(pType));
        }
        if (!itsAttr[CHARGE]) {
            Attributes::setReal(itsAttr[CHARGE], ParticleProperties::getParticleCharge(pType));
        }
    }

    if (isPhoton()) {
        reference = PartData();
        reference.setQ(0.0);
        reference.setM(0.0);
        return;
    }

    // Set up particle reference; convert all to eV for OPALX.
    double mass   = (itsAttr[MASS] ? getMass() : Physics::m_p) * Units::GeV2eV;
    double charge = itsAttr[CHARGE] ? getCharge() : 1.0;

    reference = PartData(charge, mass, 1.0);

    // Magnetic moment anomaly G = (g-2)/2 for the species (used by T-BMT).
    if (itsAttr[PARTICLE]) {
        const ParticleType pType = ParticleProperties::getParticleType(getParticleName());
        reference.setAnomaly(ParticleProperties::getParticleAnomaly(pType));
    }

    if (itsAttr[GAMMA]) {
        double gamma = Attributes::getReal(itsAttr[GAMMA]);
        if (gamma > 1.0) {
            reference.setGamma(gamma);
        } else {
            throw OpalException("Beam::update()", "\"GAMMA\" should be greater than 1.");
        }
    } else if (itsAttr[ENERGY]) {
        double energy = Attributes::getReal(itsAttr[ENERGY]) * Units::GeV2eV;
        if (energy > reference.getM()) {
            reference.setE(energy);
        } else {
            throw OpalException("Beam::update()", "\"ENERGY\" should be greater than \"MASS\".");
        }
    } else if (itsAttr[PC]) {
        double pc = Attributes::getReal(itsAttr[PC]) * Units::GeV2eV;
        if (pc > 0.0) {
            reference.setP(pc);
        } else {
            throw OpalException("Beam::update()", "\"PC\" should be greater than 0.");
        }
    }

    // Set default name.
    if (getOpalName().empty()) setOpalName("UNNAMED_BEAM");
}

void Beam::print(std::ostream& os) const {
    double charge = Attributes::getReal(itsAttr[CHARGE]);
    os << "* ************* B E A M ************************************************************ "
       << std::endl;
    os << "* BEAM        " << getOpalName() << '\n'
       << "* PARTICLE    " << Attributes::getString(itsAttr[PARTICLE]) << '\n'
       << "* REST MASS   " << Attributes::getReal(itsAttr[MASS]) << " [GeV]\n"
       << "* CHARGE      " << (charge > 0 ? '+' : '-') << "e * " << std::abs(charge) << " \n"
       << "* MOMENTUM    " << getMomentum() << " [GeV/c]\n"
       << "* BCHARGE     " << Attributes::getReal(itsAttr[BCHARGE]) << " [C]\n"
       << "* NALLOC      " << Attributes::getReal(itsAttr[NALLOC]) << '\n';
    os << "* ********************************************************************************** "
       << std::endl;
}
