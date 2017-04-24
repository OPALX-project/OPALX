#ifndef OPAL_OpalECollimator_HH
#define OPAL_OpalECollimator_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalECollimator.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalECollimator
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"

class ParticleMaterInteraction;

// Class OpalECollimator
// ------------------------------------------------------------------------
/// The ECOLLIMATOR element.

class OpalECollimator: public OpalElement {

public:

    /// The attributes of class OpalECollimator.
    enum {
        XSIZE = COMMON,  // The horizontal half-size.
        YSIZE,           // The vertical half-size.
        OUTFN,
        SIZE
    };

    /// Exemplar constructor.
    OpalECollimator();

    virtual ~OpalECollimator();

    /// Make clone.
    virtual OpalECollimator *clone(const std::string &name);

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC collimator.
    virtual void update();

private:

    // Not implemented.
    OpalECollimator(const OpalECollimator &);
    void operator=(const OpalECollimator &);

    // Clone constructor.
    OpalECollimator(const std::string &name, OpalECollimator *parent);

    ParticleMaterInteraction *parmatint_m;
};

#endif // OPAL_OpalECollimator_HH
