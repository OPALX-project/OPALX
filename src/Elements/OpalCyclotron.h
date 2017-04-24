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

    /// Need to remove common, old TYPE = COMMON, prevents the identification of the fieldmap type in the Cyclotron 
    /// element issue #84
    enum {
        TYPE,      
        GEOMETRY,  // geometry of boundary
        CYHARMON,  // The harmonic number of the cyclotron
        SYMMETRY,  // The symetry of the field
        RINIT,     // The initial radius [m]
        PRINIT,    // The initial radial momentum [pr/p0] []
        PHIINIT,   // The initial phase [deg]
        ZINIT,     // The initial z coordinate [m]
        PZINIT,    // The initial vertical momentum [pz/p0] []
        RFFREQ,    // First hamonic of the RF system [MHz]
        FMAPFN,    // The filename of the mid-plane fieldmap
        RFMAPFN,   // The filename(s) of the RF fieldmap
	RFFCFN,    // The filename(s) of coefficients for RF frequency function f(t)
	RFVCFN,    // The filename(s) of coefficients for RF voltage function v(t)
        BSCALE,    // A scalar to scale the B-field
        ESCALE,    // A scalar to scale the RF field
        TCR1,      // trim coil r1 (mm)
        TCR2,      // trim coil r2 (mm)
        MBTC,      // max bfield of trim coil (kG)
        SLPTC,     // slope of the rising edge
        TCR1V,      // trim coil r1 (mm)
        TCR2V,      // trim coil r2 (mm)
        MBTCV,      // max bfield of trim coil (kG)
        SLPTCV,     // slope of the rising edge
	RFPHI,     // the initial phase of RF field
	SUPERPOSE, // whether the electric field map are superposed or not
	MINZ,      // minimal vertical extend of the machine
	MAXZ,      // maximal vertical extend of the machine
	MINR,      // minimal radial extend of the machine
	MAXR,      // maximal radial extend of the machine
	FMLOWE,    // minimal energy of the field map
	FMHIGHE,   // maximal energy of the field map
	SPIRAL,    // flag whether or not this is a spiral inflector simulation
        SIZE
    };



    /// Exemplar constructor.
    OpalCyclotron();

    virtual ~OpalCyclotron();

    /// Make clone.
    virtual OpalCyclotron *clone(const std::string &name);

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC cavity.
    virtual void update();

private:

    // Not implemented.
    OpalCyclotron(const OpalCyclotron &);
    void operator=(const OpalCyclotron &);

    // Clone constructor.
    OpalCyclotron(const std::string &name, OpalCyclotron *parent);

    BoundaryGeometry *obgeo_m;

};

#endif // OPAL_OpalCyclotron_HH
