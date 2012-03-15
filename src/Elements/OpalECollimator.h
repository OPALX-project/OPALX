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

class SurfacePhysics;

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
        DX,             // Misalignment: translation in x direction
        DY,             // Misalignment: translation in y direction
        DZ,             // Misalignment: translation in z direction
        SIZE
    };

    /// Exemplar constructor.
    OpalECollimator();

    virtual ~OpalECollimator();

    /// Make clone.
    virtual OpalECollimator *clone(const string &name);

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC collimator.
    virtual void update();

private:

    // Not implemented.
    OpalECollimator(const OpalECollimator &);
    void operator=(const OpalECollimator &);

    // Clone constructor.
    OpalECollimator(const string &name, OpalECollimator *parent);

    SurfacePhysics *sphys_m;
};

#endif // OPAL_OpalECollimator_HH
