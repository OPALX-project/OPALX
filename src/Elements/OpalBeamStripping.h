#ifndef OPAL_OpalBeamStripping_HH
#define OPAL_OpalBeamStripping_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalSlit.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalBeamStripping
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann, Jianjun Yang $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"

class ParticleMatterInteraction;
// Class OpalBeamStripping
// ------------------------------------------------------------------------
/// The BEAMSTRIPPING element.

class OpalBeamStripping: public OpalElement {

public:

    /// The attributes of class OpalBeamStripping.
    enum {
        PRESSURE = COMMON,
		TEMPERATURE,
		CROSSSECTION,
		ENERGYCS,
        MINZ,      // minimal vertical extend of the machine
        MAXZ,      // maximal vertical extend of the machine
        MINR,      // minimal radial extend of the machine
        MAXR,      // maximal radial extend of the machine
		OUTFN,
		SIZE
    };

    /// Exemplar constructor.
    OpalBeamStripping();

    virtual ~OpalBeamStripping();

    /// Make clone.
    virtual OpalBeamStripping *clone(const std::string &name);

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC beam stripping.
    virtual void update();

private:

    // Not implemented.
    OpalBeamStripping(const OpalBeamStripping &);
    void operator=(const OpalBeamStripping &);

    // Clone constructor.
    OpalBeamStripping(const std::string &name, OpalBeamStripping *parent);
    ParticleMatterInteraction *parmatint_m;
};

#endif // OPAL_OpalBeamStripping_HH
