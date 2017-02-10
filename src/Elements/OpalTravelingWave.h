#ifndef OPAL_OpalTravelingWave_HH
#define OPAL_OpalTravelingWave_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalTravelingWave.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalTravelingWave
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"

class OpalWake;

// Class OpalTravelingWave
// ------------------------------------------------------------------------
/// The RFCAVITY element.

class OpalTravelingWave: public OpalElement {

public:

    /// The attributes of class OpalTravelingWave.
    enum {
        VOLT = COMMON,  // The peak voltage.
        DVOLT,          // The peak voltage error
        FREQ,           // The RF frequency.
        LAG,            // The phase lag.
        DLAG,           // The phase lag error
        HARMON,         // The harmonic number.
        BETARF,         // The beta_RF.
        PG,             // The RF power.
        ZSHUNT,         // The shunt impedance.
        TFILL,          // The filling time.
        FMAPFN,         // The filename of the fieldmap
        APVETO,         // Do not use this cavity in the Autophase procedure
        FAST,           // Faster but less accurate
        CAVITYTYPE,     // STANDING or TRAVELING wave structure
        NUMCELLS,       // Number of cells in a TW structure
        DESIGNENERGY,   // The mean kinetic energy at exit
        SIZE
    };

    /// Exemplar constructor.
    OpalTravelingWave();

    virtual ~OpalTravelingWave();

    /// Make clone.
    virtual OpalTravelingWave *clone(const std::string &name);

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC cavity.
    virtual void update();

private:

    // Not implemented.
    OpalTravelingWave(const OpalTravelingWave &);
    void operator=(const OpalTravelingWave &);

    // Clone constructor.
    OpalTravelingWave(const std::string &name, OpalTravelingWave *parent);

    OpalWake *owk_m;
};

#endif // OPAL_OpalTravelingWave_HH