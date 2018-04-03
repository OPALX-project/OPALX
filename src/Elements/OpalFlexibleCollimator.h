#ifndef OPAL_OpalFlexibleCollimator_HH
#define OPAL_OpalFlexibleCollimator_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalFlexibleCollimator.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalFlexibleCollimator
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"

class ParticleMatterInteraction;

// Class OpalFlexibleCollimator
// ------------------------------------------------------------------------
/// The ECOLLIMATOR element.

class OpalFlexibleCollimator: public OpalElement {

public:

    /// The attributes of class OpalFlexibleCollimator.
    enum {
        FNAME = COMMON,  // The horizontal half-size.
        OUTFN,
        SIZE
    };

    /// Exemplar constructor.
    OpalFlexibleCollimator();

    virtual ~OpalFlexibleCollimator();

    /// Make clone.
    virtual OpalFlexibleCollimator *clone(const std::string &name);

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC collimator.
    virtual void update();

private:

    // Not implemented.
    OpalFlexibleCollimator(const OpalFlexibleCollimator &);
    void operator=(const OpalFlexibleCollimator &);

    // Clone constructor.
    OpalFlexibleCollimator(const std::string &name, OpalFlexibleCollimator *parent);

    ParticleMatterInteraction *partMatInt_m;
};

#endif // OPAL_OpalFlexibleCollimator_HH