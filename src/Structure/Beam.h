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
#ifndef OPAL_Beam_HH
#define OPAL_Beam_HH

#include "AbstractObjects/Definition.h"
#include "Algorithms/PartData.h"

#include <ostream>
#include <string>
#include <vector>

class Inform;

class Beam: public Definition {

public:
    /// Exemplar constructor.
    Beam();

    virtual ~Beam();

    /// Test if replacement is allowed.
    //  Can replace only by another BEAM.
    virtual bool canReplaceBy(Object* object);

    /// Make clone.
    virtual Beam* clone(const std::string& name);

    /// Check the BEAM data.
    virtual void execute();

    /// Find named BEAM.
    static Beam* find(const std::string& name);

    /// Return the number of (macro)particles
    size_t getNumberOfParticles() const;

    /// Return the embedded CLASSIC PartData.
    const PartData& getReference() const;

    /// Return the beam current in A
    double getCurrent() const;

    /// Return the charge number in elementary charge
    double getCharge() const;

    /// Return the beam frequency in MHz
    double getFrequency() const;

    /// Return Particle's name
    std::string getParticleName() const;

    /// True if this beam is configured as a photon beam.
    bool isPhoton() const;

    /// Return Particle's rest mass in GeV
    double getMass() const;

    // Returns Momentum in GeV/c
    double getMomentum() const;

    /// Charge per macro particle in C
    double getChargePerParticle() const;

    /// Mass per macro particle in GeV/c^2
    double getMassPerParticle() const;

    /// Return the name of the EMISSIONSOURCELIST linked to this beam.
    /// Throws if `SOURCES` is not set.
    std::string getEmissionSourceListName() const;

    /// Return the configured global process names for this beam.
    std::vector<std::string> getGlobalProcessNames() const;

    /// Return the name of the daughter beam (for decay products), or empty if not set.
    std::string getDaughterBeamName() const;

    /// True if PC, ENERGY, or GAMMA was explicitly provided by the user.
    bool hasExplicitEnergy() const;

    /// Update the BEAM data.
    virtual void update();

    void print(std::ostream& os) const;

private:
    // Not implemented.
    Beam(const Beam&);
    void operator=(const Beam&);

    // Clone constructor.
    Beam(const std::string& name, Beam* parent);

    // The particle reference data.
    PartData reference;

    // The converstion from GeV to eV.
    static const double energy_scale;
};

inline std::ostream &operator<<(std::ostream& os, const Beam& b) {
    b.print(os);
    return os;
}

#endif // OPAL_Beam_HH
