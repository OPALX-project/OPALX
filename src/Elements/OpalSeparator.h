#ifndef OPAL_OpalSeparator_HH
#define OPAL_OpalSeparator_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalSeparator.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalSeparator
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:40 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"


// Class OpalSeparator
// ------------------------------------------------------------------------
/// The ELSEPARATOR element.

class OpalSeparator: public OpalElement {

public:

    /// The attributes of class OpalSeparator.
    enum {
        EX = COMMON,  // The horizontal field.
        EY,           // The vertical field.
        SIZE
    };

    /// Exemplar constructor.
    OpalSeparator();

    virtual ~OpalSeparator();

    /// Make clone.
    virtual OpalSeparator *clone(const string &name);

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC separator.
    virtual void update();

private:

    // Not implemented.
    OpalSeparator(const OpalSeparator &);
    void operator=(const OpalSeparator &);

    // Clone constructor.
    OpalSeparator(const string &name, OpalSeparator *parent);
};

#endif // OPAL_OpalSeparator_HH
