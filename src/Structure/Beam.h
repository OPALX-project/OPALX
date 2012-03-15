#ifndef OPAL_Beam_HH
#define OPAL_Beam_HH

// ------------------------------------------------------------------------
// $RCSfile: Beam.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Beam
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:44 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Definition.h"
#include "Algorithms/PartData.h"


// Class Beam
// ------------------------------------------------------------------------
/// The BEAM definition.
//  A BEAM definition is used by most physics commands to define the
//  particle charge and the reference momentum, together with some other
//  data.

class Beam: public Definition {

public:

    /// Exemplar constructor.
    Beam();

    virtual ~Beam();

    /// Test if replacement is allowed.
    //  Can replace only by another BEAM.
    virtual bool canReplaceBy(Object *object);

    /// Make clone.
    virtual Beam *clone(const string &name);

    /// Check the BEAM data.
    virtual void execute();

    /// Find named BEAM.
    static Beam *find(const string &name);

    /// Return emittance for mode 1.
    double getEX() const;

    /// Return emittance for mode 2.
    double getEY() const;

    /// Return emittance for mode 3.
    double getET() const;

    //ff => get gamma value
    double getGamma() const;

    //ff => get PC value
    double getPC() const;

    /// Return the number of (macro)particles
    size_t getNumberOfParticles();

    /// Return the embedded CLASSIC PartData.
    const PartData &getReference() const;

    /// Return the beam current in A
    double getCurrent() const;

    /// Return the charge in unit protons
    double getCharge() const;

    /// Return the beam frequency in Hz
    double getFrequency() const;

    /// Return Particle's name
    string getParticleName() const;

    /// Return Particle's rest mass in unit protons
    double getMass() const;

    /// Store emittance for mode 1.
    void setEX(double);

    /// Store emittance for mode 2.
    void setEY(double);

    /// Store emittance for mode 3.
    void setET(double);

    /// Update the BEAM data.
    virtual void update();

    /// Print the TFS descriptors for the beam.
    void tfsDescriptors(std::ostream &os) const;

    Inform &print(Inform &os) const;

private:

    // Not implemented.
    Beam(const Beam &);
    void operator=(const Beam &);

    // Clone constructor.
    Beam(const string &name, Beam *parent);

    // The particle reference data.
    PartData reference;

    // The converstion from GeV to eV.
    static const double energy_scale;
};

inline Inform &operator<<(Inform &os, const Beam &b) {
    return b.print(os);
}


#endif // OPAL_Beam_HH
