#ifndef OPAL_OpalPatch_HH
#define OPAL_OpalPatch_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalPatch.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalPatch
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:40 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"


// Class OpalPatch
// ------------------------------------------------------------------------
/// The PATCH element.

class OpalPatch: public OpalElement {

public:

    /// The attributes of class OpalPatch.
    enum {
        DX = COMMON,  // The horizontal orbit displacement.
        DY,           // The vertical orbit displacement.
        DS,           // The longitudinal orbit displacement.
        VX,           // The rotation around the x-axis.
        VY,           // The rotation around the y-axis.
        VS,           // The rotation around the s-axis.
        SIZE
    };

    /// Exemplar constructor.
    OpalPatch();

    virtual ~OpalPatch();

    /// Make clone.
    virtual OpalPatch *clone(const string &name);

    /// Test for patch.
    //  Return true.
    virtual bool isPatch() const;

    /// Update the embedded CLASSIC patch.
    virtual void update();

private:

    // Not implemented.
    OpalPatch(const OpalPatch &);
    void operator=(const OpalPatch &);

    // Clone constructor.
    OpalPatch(const string &name, OpalPatch *parent);
};

#endif // OPAL_OpalPatch_HH
