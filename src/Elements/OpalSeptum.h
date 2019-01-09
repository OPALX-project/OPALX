#ifndef OPAL_OpalSeptum_HH
#define OPAL_OpalSeptum_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalSeptum.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalSeptum
//
// ------------------------------------------------------------------------
//
// $Date: 2009/09/21 10:08:06 $
// $Author: bi $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"

class OpalWake;

// Class OpalSeptum
// ------------------------------------------------------------------------
/// The Septum element.

class OpalSeptum: public OpalElement {

public:

    /// The attributes of class OpalSeptum.
    enum {
        XSTART = COMMON, // Start of x coordinate
        XEND,            // End of x coordinate
        YSTART,          // Start of y coordinate
        YEND,            // End of y coordinate
        WIDTH,           // Width of the septum
        OUTFN,           // Output filename
        SIZE
    };
    /// Exemplar constructor.
    OpalSeptum();

    virtual ~OpalSeptum();

    /// Make clone.
    virtual OpalSeptum *clone(const std::string &name);

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC septum.
    virtual void update();

private:

    // Not implemented.
    OpalSeptum(const OpalSeptum &);
    void operator=(const OpalSeptum &);

    // Clone constructor.
    OpalSeptum(const std::string &name, OpalSeptum *parent);

    OpalWake *owk_m;
};

#endif // OPAL_OpalSeptum_HH
