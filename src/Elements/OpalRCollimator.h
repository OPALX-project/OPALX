#ifndef OPAL_OpalRCollimator_HH
#define OPAL_OpalRCollimator_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalRCollimator.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalRCollimator
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:40 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"

class SurfacePhysics;

// Class OpalRCollimator
// ------------------------------------------------------------------------
/// The RCOLLIMATOR element.

class OpalRCollimator: public OpalElement {

public:

    /// The attributes of class OpalRCollimator.
    enum {
        XSIZE = COMMON,  // The horizontal half-size.
        YSIZE,           // The vertical half-size.
        OUTFN,
        SIZE
    };

    /// Exemplar constructor.
    OpalRCollimator();

    virtual ~OpalRCollimator();

    /// Make clone.
    virtual OpalRCollimator *clone(const string &name);

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC collimator.
    virtual void update();

private:

    // Not implemented.
    OpalRCollimator(const OpalRCollimator &);
    void operator=(const OpalRCollimator &);

    // Clone constructor.
    OpalRCollimator(const string &name, OpalRCollimator *parent);

    SurfacePhysics *sphys_m;
};

#endif // OPAL_OpalRCollimator_HH
