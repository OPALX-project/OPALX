#ifndef OPAL_OpalCCollimator_HH
#define OPAL_OpalCCollimator_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalSlit.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalCCollimator
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"

class SurfacePhysics;
// Class OpalCCollimator
// ------------------------------------------------------------------------
/// The CCOLLIMATOR element.

class OpalCCollimator: public OpalElement {

public:

    /// The attributes of class OpalCCollimator.
    enum {
        ANGSTART = COMMON,  // Start of angle in rad.
        ANGEND,       //End of angle in rad.
        RSTART,  //Start of radius in mm.
        REND,           // End of radius in mm.
        WIDTH, //The width of collimator
        OUTFN,
        SIZE
    };

    /// Exemplar constructor.
    OpalCCollimator();

    virtual ~OpalCCollimator();

    /// Make clone.
    virtual OpalCCollimator *clone(const string &name);

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC collimator.
    virtual void update();

private:

    // Not implemented.
    OpalCCollimator(const OpalCCollimator &);
    void operator=(const OpalCCollimator &);

    // Clone constructor.
    OpalCCollimator(const string &name, OpalCCollimator *parent);
    SurfacePhysics *sphys_m;
};

#endif // OPAL_OpalCCollimator_HH
