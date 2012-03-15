#ifndef OPAL_OpalQuadrupole_HH
#define OPAL_OpalQuadrupole_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalQuadrupole.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalQuadrupole
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"

class SurfacePhysics;

// Class OpalQuadrupole
// ------------------------------------------------------------------------
/// The QUADRUPOLE element.

class OpalQuadrupole: public OpalElement {

public:

    /// The attributes of class OpalQuadrupole.
    enum {
        K1 = COMMON,  // The normal quadrupole coefficient.
        K1S,          // The skew quadrupole coefficient.
        SIZE
    };

    /// Exemplar constructor.
    OpalQuadrupole();

    virtual ~OpalQuadrupole();

    /// Make clone.
    virtual OpalQuadrupole *clone(const string &name);

    /// Print the quadrupole.
    //  Handle printing in OPAL-8 format.
    virtual void print(std::ostream &) const;

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC multipole.
    virtual void update();

private:

    // Not implemented.
    OpalQuadrupole(const OpalQuadrupole &);
    void operator=(const OpalQuadrupole &);

    // Clone constructor.
    OpalQuadrupole(const string &name, OpalQuadrupole *parent);

    SurfacePhysics *sphys_m;
};

#endif // OPAL_OpalQuadrupole_HH
