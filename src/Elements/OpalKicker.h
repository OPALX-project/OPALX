#ifndef OPAL_OpalKicker_HH
#define OPAL_OpalKicker_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalKicker.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalKicker
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:32:23 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"


// Class OpalKicker
// ------------------------------------------------------------------------
/// The KICKER element.
//  Note the sign convention:  Positive kicks bend particles to positive x or
//  y respectively.

class OpalKicker: public OpalElement {

public:

    /// The attributes of class OpalKicker.
    enum {
        HKICK = COMMON,  // The horizontal kicker strength.
        VKICK,           // The vertical kicker strength.
        DESIGNENERGY,    // The mean kinetic energy at exit
        K0,              // The normal dipole field
        K0S,             // The skew dipole field
        SIZE
    };

    /// Exemplar constructor.
    OpalKicker();

    virtual ~OpalKicker();

    /// Make clone.
    virtual OpalKicker *clone(const std::string &name);


    // JMJ 18/12/2000 Following method not needed, commented out, delete after next CVS commit.
    //BEGIN JMJ 15/12/2000, adding missing print method
    // Print the kicker
    //  Handle printing in OPAL-8 format.
    //  virtual void print(std::ostream &) const;
    //END   JMJ 15/12/2000, adding missing print method

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC corrector.
    virtual void update();

private:

    // Not implemented.
    OpalKicker(const OpalKicker &);
    void operator=(const OpalKicker &);

    // Clone constructor.
    OpalKicker(const std::string &name, OpalKicker *parent);
};

#endif // OPAL_OpalKicker_HH
