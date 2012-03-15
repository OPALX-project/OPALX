#ifndef OPAL_OpalBeamBeam_HH
#define OPAL_OpalBeamBeam_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalBeamBeam.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalBeamBeam
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:38 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"


// Class OpalBeamBeam
// ------------------------------------------------------------------------
/// The BEAMBEAM element.

class OpalBeamBeam: public OpalElement {

public:

    /// The attributes of class OpalBeamBeam.
    enum {
        DX = COMMON, // Horizontal displacement of opposite beam.
        DY,          // Vertical displacement of opposite beam.
        SIGX,        // Horizontal extent of opposite beam.
        SIGY,        // Vertical extent of opposite beam.
        CHARGE,      // Particle charge in opposite beam.
        NPART,       // Number of particles in opposite beam.
        SIZE
    };

    /// Exemplar constructor.
    OpalBeamBeam();

    virtual ~OpalBeamBeam();

    /// Make clone.
    virtual OpalBeamBeam *clone(const string &name);

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC beam-beam element.
    virtual void update();

private:

    // Not implemented.
    OpalBeamBeam(const OpalBeamBeam &);
    void operator=(const OpalBeamBeam &);

    // Clone constructor.
    OpalBeamBeam(const string &name, OpalBeamBeam *parent);
};

#endif // OPAL_OpalBeamBeam_H
