// ------------------------------------------------------------------------
// $RCSfile: Surveyor.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Surveyor
//   The visitor class for printing the survey of a beamline.
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2000/05/03 08:16:05 $
// $Author: mad $
//
// ------------------------------------------------------------------------

#include "Algorithms/Surveyor.h"
#include "AbsBeamline/AttributeSet.h"
#include "AbsBeamline/Patch.h"


// Class Surveyor
// ------------------------------------------------------------------------

Surveyor::Surveyor(Beamline &beamline, bool revTrack):
  DefaultVisitor(beamline, false, revTrack),
  itsMap()
{}


Surveyor::Surveyor(Beamline &beamline, double x0, double y0, double z0,
		   double theta0, double phi0, double psi0,
		   bool revTrack):
  DefaultVisitor(beamline, false, revTrack),
  itsMap(x0, y0, z0, theta0, phi0, psi0)
{}


Surveyor::Surveyor(Beamline &beamline, const Euclid3D &init, bool revTrack):
  DefaultVisitor(beamline, false, revTrack), itsMap(init)
{}


Surveyor::~Surveyor()
{}


void Surveyor::getMap(Euclid3D &map) const
{
  map = itsMap;
}


void Surveyor::setMap(const Euclid3D &map)
{
  itsMap = map;
}


void Surveyor::visitPatch(const Patch &patch)
{
  Euclid3D elementMap = patch.getPatch();
  if (back_path) elementMap = Inverse(elementMap);
  itsMap.dotBy(elementMap);
}


void Surveyor::applyDefault(const ElementBase &element)
{
  Euclid3D elementMap = element.getGeometry().getTotalTransform();
  if (back_path) elementMap = Inverse(elementMap);
  itsMap.dotBy(elementMap);
}
