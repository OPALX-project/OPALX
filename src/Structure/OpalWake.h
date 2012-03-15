#ifndef OPAL_Wake_HH
#define OPAL_Wake_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalWake.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalWake
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:44 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Definition.h"
#include "Algorithms/PartData.h"

class ElementBase;
class WakeFunction;
// Class OpalWake
// ------------------------------------------------------------------------
/// The WAKE definition.
//  A WAKE definition is used by most physics commands to define the
//  particle charge and the reference momentum, together with some other
//  data.

class OpalWake: public Definition {

public:

    /// Exemplar constructor.
    OpalWake();

    virtual ~OpalWake();

    /// Test if replacement is allowed.
    //  Can replace only by another WAKE.
    virtual bool canReplaceBy(Object *object);

    /// Make clone.
    virtual OpalWake *clone(const string &name);

    /// Check the WAKE data.
    virtual void execute();

    /// Find named WAKE.
    static OpalWake *find(const string &name);

    /// Update the WAKE data.
    virtual void update();

    /// Print the TFS descriptors for the wake
    void tfsDescriptors(std::ostream &os) const;

    Inform &print(Inform &os) const;

    int getNumberOfBins();

    void initWakefunction(ElementBase &element);

    void updateElement(ElementBase *element);
    WakeFunction *wf_m;

private:

    // Not implemented.
    OpalWake(const OpalWake &);
    void operator=(const OpalWake &);

    // Clone constructor.
    OpalWake(const string &name, OpalWake *parent);

    // The particle reference data.
    PartData reference;

    // The conversion from GeV to eV.
    static const double energy_scale;

    // the element the wake is attached to
    ElementBase *itsElement_m;

};

inline Inform &operator<<(Inform &os, const OpalWake &b) {
    return b.print(os);
}

#endif // OPAL_Wake_HH
