#ifndef OPAL_OpalStripper_HH
#define OPAL_OpalStripper_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalStripper.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalStripper
//
// ------------------------------------------------------------------------
//
// $Date: 2011/07/08 08:21:00 $
// $Author: Jianjun Yang $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"

// Class OpalStripper
// ------------------------------------------------------------------------
/// The Stripper element.


class OpalStripper: public OpalElement {

public:

    /// The attributes of class OpalStripper.
    enum {
        XSTART,           // Start of x coordinate
        XEND,             // End of x coordinate
        YSTART,           // Start of y coordinate
        YEND,             // End of y coordinate
        WIDTH,            // Width of the probe
	OPCHARGE,           // Charge number of the outcome particle
	OPMASS,             // Mass of the outcome particle
	STOP,
        SIZE
    };

    OpalStripper();

    virtual ~OpalStripper();

    /// Make clone.
    virtual OpalStripper *clone(const string &name);

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC septum.
    virtual void update();

private:

    // Not implemented.
    OpalStripper(const OpalStripper &);
    void operator=(const OpalStripper &);

    // Clone constructor.
    OpalStripper(const string &name, OpalStripper *parent);

};

#endif // OPAL_OpalStripper_HH
