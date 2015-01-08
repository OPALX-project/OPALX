#ifndef OPAL_OpalDegrader_HH
#define OPAL_OpalDegrader_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalDegrader.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalDegrader
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"

class SurfacePhysics;

// Class OpalDegrader
// ------------------------------------------------------------------------
/// The DEGRADER element.

class OpalDegrader: public OpalElement {

public:

    /// The attributes of class OpalDegrader.
    enum {
        XSIZE = COMMON,  // not used
        YSIZE,           // not used
//        ZSIZE,           // the thickness of the degrader
        OUTFN,
        DX,             // Misalignment: translation in x direction
        DY,             // Misalignment: translation in y direction
        SIZE
    };

    /// Exemplar constructor.
    OpalDegrader();

    virtual ~OpalDegrader();

    /// Make clone.
    virtual OpalDegrader *clone(const std::string &name);

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC collimator.
    virtual void update();

private:

    // Not implemented.
    OpalDegrader(const OpalDegrader &);
    void operator=(const OpalDegrader &);

    // Clone constructor.
    OpalDegrader(const std::string &name, OpalDegrader *parent);

    SurfacePhysics *sphys_m;
};

#endif // OPAL_OpalDegrader_HH
