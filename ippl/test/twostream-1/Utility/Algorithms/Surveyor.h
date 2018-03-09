#ifndef MAD_Surveyor_HH
#define MAD_Surveyor_HH

// ------------------------------------------------------------------------
// $RCSfile: Surveyor.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Surveyor
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Algorithms/DefaultVisitor.h"
#include "BeamlineGeometry/Euclid3D.h"


// Class Surveyor
// ------------------------------------------------------------------------
//: Survey algorithm.
//  This visitor class computes the survey of a beam line.

class Surveyor: public DefaultVisitor {

public:

  //: Constructor.
  //  Assume zero initial conditions.
  //  The beam line to be tracked is [b]bl[/b].
  //  If [b]revTrack[/b] is true, track from s = C to s = 0.
  Surveyor(Beamline &bl, bool revTrack);

  //: Constructor.
  //  Use given initial conditions.
  //  The beam line to be tracked is [b]bl[/b].
  //  If [b]revTrack[/b] is true, track from s = C to s = 0.
  Surveyor(Beamline &bl, double x0, double y0, double z0,
	   double theta0, double phi0, double psi0,
	   bool revTrack = false);

  //: Constructor.
  //  Use given initial conditions in terms of an Euclid3D object.
  //  The beam line to be tracked is [b]bl[/b].
  //  If [b]revTrack[/b] is true, track from s = C to s = 0.
  Surveyor(Beamline &bl, const Euclid3D &init, bool revTrack);

  virtual ~Surveyor();

  //: Return accumulated map.
  void getMap(Euclid3D &) const;

  //: Reset accumulated map for restart.
  void setMap(const Euclid3D &);


  //: Apply the algorithm to a patch.
  virtual void visitPatch(const Patch &pat);

private:

  // Not implemented.
  Surveyor();
  Surveyor(const Surveyor &);
  void operator=(const Surveyor &);

  //: Default action.
  //  Apply the default to all element (advance the position and direction).
  //  All visitXXX() methods except visitPatch() call applyDefault() which
  //  is overridden here to propagate the survey through the element.
  virtual void applyDefault(const ElementBase &element);

  //: The accumulated survey map.
  Euclid3D itsMap;
};

#endif // __Surveyor_HH
