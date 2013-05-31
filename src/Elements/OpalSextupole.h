#ifndef OPAL_OpalSextupole_HH
#define OPAL_OpalSextupole_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalSextupole.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalSextupole
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:40 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"


// Class OpalSextupole
// ------------------------------------------------------------------------
/// The SEXTUPOLE element.

class OpalSextupole: public OpalElement {

public:

    /// The attributes of class OpalSextupole.
    enum {
        K2 = COMMON,  // The normal sextupole strength.
        K2S,          // The skew sextupole strength.
        SIZE
    };

    /// Exemplar constructor.
    OpalSextupole();

    virtual ~OpalSextupole();

    /// Make clone.
    virtual OpalSextupole *clone(const string &name);

    /// Print the sextupole.
    //  Handle printing in OPAL-8 format.
    virtual void print(std::ostream &) const;

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC multipole.
    virtual void update();

private:

    // Not implemented.
    OpalSextupole(const OpalSextupole &);
    void operator=(const OpalSextupole &);

    // Clone constructor.
    OpalSextupole(const string &name, OpalSextupole *parent);
};

#endif // OPAL_OpalSextupole_HH
