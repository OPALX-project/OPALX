#ifndef OPAL_OpalOctupole_HH
#define OPAL_OpalOctupole_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalOctupole.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalOctupole
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"


// Class OpalOctupole
// ------------------------------------------------------------------------
/// The OCTUPOLE element.

class OpalOctupole: public OpalElement {

public:

    /// The attributes of class OpalOctupole.
    enum {
        K3 = COMMON,  // The normal octupole coefficient.
        DK3,          // The normal octupole coefficient error.
        K3S,          // The skew octupole coefficient.
        DK3S,          // The skew octupole coefficient error.
        SIZE
    };

    /// Exemplar constructor.
    OpalOctupole();

    virtual ~OpalOctupole();

    /// Make clone.
    virtual OpalOctupole *clone(const std::string &name);

    /// Print the element.
    //  Handle printing in OPAL-8 format.
    virtual void print(std::ostream &) const;

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC multipole.
    virtual void update();

private:

    // Not implemented.
    OpalOctupole(const OpalOctupole &);
    void operator=(const OpalOctupole &);

    // Clone constructor.
    OpalOctupole(const std::string &name, OpalOctupole *parent);
};

#endif // OPAL_OpalOctupole_HH
