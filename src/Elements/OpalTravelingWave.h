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


// Class OpalTravelingWave
// ------------------------------------------------------------------------
/// The RFCAVITY element.

class OpalTravelingWave: public OpalElement {

public:

  /// The attributes of class OpalTravelingWave.
  enum {
    VOLT = COMMON,  // The peak voltage.
    FREQ,           // The RF frequency.
    LAG,            // The phase lag.
    HARMON,         // The harmonic number.
    BETARF,         // The beta_RF.
    PG,             // The RF power.
    ZSHUNT,         // The shunt impedance.
    TFILL,          // The filling time.
    FMAPFN,         // The filename of the fieldmap
    FAST,           // Faster but less accurate
    CAVITYTYPE,     // STANDING or TRAVELING wave structure
    NUMCELLS,       // Number of cells in a TW structure
    DX,             // Misalignment: translation in x direction
    DY,             // Misalignment: translation in y direction
    DZ,             // Misalignment: translation in z direction
    SIZE
  };

  /// Exemplar constructor.
  OpalTravelingWave();

  virtual ~OpalTravelingWave();

  /// Make clone.
  virtual OpalTravelingWave *clone(const string &name);

  /// Fill in all registered attributes.
  virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);
  
  /// Update the embedded CLASSIC cavity.
  virtual void update();

private:

  // Not implemented.
  OpalTravelingWave(const OpalTravelingWave &);
  void operator=(const OpalTravelingWave &);

  // Clone constructor.
  OpalTravelingWave(const string &name, OpalTravelingWave *parent);
};

#endif // OPAL_OpalTravelingWave_HH
