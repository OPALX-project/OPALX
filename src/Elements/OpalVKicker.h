#ifndef OPAL_OpalVKicker_HH
#define OPAL_OpalVKicker_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalVKicker.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalVKicker
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:32:24 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"


// Class OpalVKicker
// ------------------------------------------------------------------------
/// The VKICKER element.
//  Note the sign convention:  A positive kick bend particles to positive y.

class OpalVKicker: public OpalElement {

public:

    /// The attributes of class OpalVKicker.
    enum {
        KICK = COMMON,  // The kicker strength.
        SIZE
    };

    /// Exemplar constructor.
    OpalVKicker();

    virtual ~OpalVKicker();

    /// Make clone.
    virtual OpalVKicker *clone(const string &name);


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
    OpalVKicker(const OpalVKicker &);
    void operator=(const OpalVKicker &);

    // Clone constructor.
    OpalVKicker(const string &name, OpalVKicker *parent);
};

#endif // OPAL_OpalVKicker_HH
