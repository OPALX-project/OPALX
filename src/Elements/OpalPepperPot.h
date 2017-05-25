#ifndef OPAL_OpalPeperPot_HH
#define OPAL_OpalPeperPot_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalPeperPot.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalPepperPot
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"


// Class OpalPepperPot
// ------------------------------------------------------------------------
/// The PEPPERPOT element.

class ParticleMatterInteraction;

class OpalPepperPot: public OpalElement {

public:

    /// The attributes of class OpalPepperPot.
    enum {
        R = COMMON,  // The horizontal half-size of a hole
        PITCH,        // The separation of the pepperpot holes
        NHOLX,
        NHOLY,
        XSIZE,
        YSIZE,
        OUTFN,
        SIZE
    };

    /// Exemplar constructor.
    OpalPepperPot();

    virtual ~OpalPepperPot();

    /// Make clone.
    virtual OpalPepperPot *clone(const std::string &name);

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC collimator.
    virtual void update();

private:

    // Not implemented.
    OpalPepperPot(const OpalPepperPot &);
    void operator=(const OpalPepperPot &);

    // Clone constructor.
    OpalPepperPot(const std::string &name, OpalPepperPot *parent);

    ParticleMatterInteraction *parmatint_m;
};

#endif // OPAL_OpalPepperPot_HH
