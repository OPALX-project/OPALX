#ifndef OPAL_PARTICLEMATERINTERACTION_HH
#define OPAL_PARTICLEMATERINTERACTION_HH

// ------------------------------------------------------------------------
// $RCSfile: Wake.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Wake
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:44 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Definition.h"
#include "Algorithms/PartData.h"
#include "Solvers/ParticleMaterInteractionHandler.hh"
class ElementBase;
class Inform;

// Class Wake
// ------------------------------------------------------------------------
/// The WAKE definition.
//  A WAKE definition is used by most physics commands to define the
//  particle charge and the reference momentum, together with some other
//  data.

class ParticleMaterInteraction: public Definition {

public:

    /// Exemplar constructor.
    ParticleMaterInteraction();

    virtual ~ParticleMaterInteraction();

    /// Test if replacement is allowed.
    //  Can replace only by another WAKE.
    virtual bool canReplaceBy(Object *object);

    /// Make clone.
    virtual ParticleMaterInteraction *clone(const std::string &name);

    /// Check the PARTICLEMATERINTERACTION data.
    virtual void execute();

    /// Find named PARTICLEMATERINTERACTION.
    static ParticleMaterInteraction *find(const std::string &name);

    /// Update the PARTICLEMATERINTERACTION data.
    virtual void update();

    void print(std::ostream &os) const;

    void initParticleMaterInteractionHandler(ElementBase &element);

    void updateElement(ElementBase *element);

    ParticleMaterInteractionHandler *handler_m;

private:

    // Not implemented.
    ParticleMaterInteraction(const ParticleMaterInteraction &);
    void operator=(const ParticleMaterInteraction &);

    // Clone constructor.
    ParticleMaterInteraction(const std::string &name, ParticleMaterInteraction *parent);

    // The particle reference data.
    PartData reference;

    // The conversion from GeV to eV.
    static const double energy_scale;

    // the element the particle mater interaction is attached to
    ElementBase *itsElement_m;
    std::string material_m;

};

inline std::ostream &operator<<(std::ostream &os, const ParticleMaterInteraction &b) {
    b.print(os);
    return os;
}

#endif // OPAL_PARTICLEMATERINTERACTION_HH
