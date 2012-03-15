#ifndef OPAL_OpalMultipole_HH
#define OPAL_OpalMultipole_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalMultipole.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalMultipole
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"


// Class OpalMultipole
// ------------------------------------------------------------------------
/// The MULTIPOLE element.

class OpalMultipole: public OpalElement {

public:

    /// The attributes of class OpalMultipole.
    enum {
        KN = COMMON,  // The normal field components.
        KS,           // The skewed field components.
        SIZE
    };

    /// Exemplar constructor.
    OpalMultipole();

    virtual ~OpalMultipole();

    /// Make clone.
    virtual OpalMultipole *clone(const string &name);

    /// Read element from the DOOM data base.
    virtual void doomGet(const DoomReader &);

    /// Write element to the DOOM data base.
    virtual void doomPut(DoomWriter &) const;

    /// Print the object.
    //  Handle printing in OPAL-8 format.
    virtual void print(std::ostream &) const;

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC multipole.
    virtual void update();

private:

    // Not implemented.
    OpalMultipole(const OpalMultipole &);
    void operator=(const OpalMultipole &);

    // Clone constructor.
    OpalMultipole(const string &name, OpalMultipole *parent);
};

#endif // OPAL_OpalMultipole_HH
