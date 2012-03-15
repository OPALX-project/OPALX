#ifndef OPAL_OpalBeamBeam3D_HH
#define OPAL_OpalBeamBeam3D_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalBeamBeam3D.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalBeamBeam3D
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"


// Class OpalBeamBeam3D
// ------------------------------------------------------------------------
/// The BEAMINT element.

class OpalBeamBeam3D: public OpalElement {

public:

  /// The attributes of class OpalBeamBeam3D.
  enum {
    SIGX = COMMON,  // Horizontal standard deviation of strong beam.
    SIGY,           // Vertical standard deviation of strong beam.
    XMA,            // Horizontal displacement of strong beam.
    YMA,            // Vertical displacement of strong beam.
    ZMA,            // Longitudinal displacement of strong beam.
	            
    SLICES,         // Number of slices in strong beam.
    ANGLE,          // Horizontal crossing angle.
    XIYN,           // Tune shift factor.
    FAST,           // If true, use tables for error function.
    ALFXS,          // Alpha* for the strong beam.
    ALFYS,          
    DXS,            // Dispersion functions for strong beam.
    DPXS,           
    DYS,            
    DPYS,           
    EXS,            // Emittances.
    EYS,            
    SIGTS,          // Bunch length of strong beam.
    SIGES,          // Energy spread of strong beam.
    SIZE
  };

  /// Exemplar constructor.
  OpalBeamBeam3D();

  virtual ~OpalBeamBeam3D();

  /// Make clone.
  virtual OpalBeamBeam3D *clone(const string &name);

  /// Fill in all registered attributes.
  virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

  /// Update the embedded CLASSIC element.
  virtual void update();

private:

  // Not implemented.
  OpalBeamBeam3D(const OpalBeamBeam3D &);
  void operator=(const OpalBeamBeam3D &);

  // Clone constructor.
  OpalBeamBeam3D(const string &name, OpalBeamBeam3D *parent);
};

#endif // OPAL_OpalBeamBeam3D_H
