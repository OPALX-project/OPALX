#ifndef OPAL_OpalHKicker_HH
#define OPAL_OpalHKicker_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalHKicker.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalHKicker
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"


// Class OpalHKicker
// ------------------------------------------------------------------------
/// The HKICKER element.
//  Note the sign convention:  A positive kick bend particles to positive x.


class OpalHKicker: public OpalElement {

public:

    /// The attributes of class OpalHKicker.
    enum {
        KICK = COMMON,  // The kicker strength.
        DESIGNENERGY,   // The mean kinetic energy at exit
        K0,             // The magnetic field
        SIZE
    };

    /// Exemplar constructor.
    OpalHKicker();

    virtual ~OpalHKicker();

    /// Make clone.
    virtual OpalHKicker *clone(const std::string &name);

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC corrector.
    virtual void update();

private:

    // Not implemented.
    OpalHKicker(const OpalHKicker &);
    void operator=(const OpalHKicker &);

    // Clone constructor.
    OpalHKicker(const std::string &name, OpalHKicker *parent);
};

#endif // OPAL_OpalHKicker_HH
