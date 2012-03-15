#ifndef OPAL_OpalCyclotron_HH
#define OPAL_OpalCyclotron_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalCyclotron.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalCyclotron
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"


// Class OpalCyclotron
// ------------------------------------------------------------------------
/// The OpalCyclotron element.

class OpalCyclotron: public OpalElement {

public:

    /// The attributes of class OpalCyclotron.
    enum {
        TYPE,
        CYHARMON,         // The harmonic number of the cyclotron
        SYMMETRY,         // The symetry of the field     
        RINIT,             // The initial radius [m]
        PRINIT,             // The initial radial momenta [pr/p0] []
        PHIINIT,               // The initial phase [deg]
        RFFREQ,                // First hamonic of the RF system
        FMAPFN,                // The filename of the fieldmap
        BSCALE,                // A scalar to scale the B-field
        TCR1,    //trim coil r1 (mm)
        TCR2,    //trim coil r2 (mm)
        MBTC,    //max bfield of trim coil (kG)
	SLPTC,    //slope of the rising edge
        SIZE

    };



    /// Exemplar constructor.
    OpalCyclotron();

    virtual ~OpalCyclotron();

    /// Make clone.
    virtual OpalCyclotron *clone(const string &name);

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);
  
    /// Update the embedded CLASSIC cavity.
    virtual void update();

private:

    // Not implemented.
    OpalCyclotron(const OpalCyclotron &);
    void operator=(const OpalCyclotron &);

    // Clone constructor.
    OpalCyclotron(const string &name, OpalCyclotron *parent);
};

#endif // OPAL_OpalCyclotron_HH
