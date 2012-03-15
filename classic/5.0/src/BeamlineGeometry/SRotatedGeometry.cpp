// ------------------------------------------------------------------------
// $RCSfile: SRotatedGeometry.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: SRotatedGeometry
//    A Geometry which wraps another arbitrary geometry in two s-rotations.
//
// ------------------------------------------------------------------------
// Class category: BeamlineGeometry
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineGeometry/SRotatedGeometry.h"
#include "BeamlineGeometry/Euclid3D.h"


// Class SRotatedGeometry.
// ------------------------------------------------------------------------

SRotatedGeometry::SRotatedGeometry(const Geometry &g,
				   double sin, double sout):
  srotIn(sin), srotOut(sout), geom(g)
{}


SRotatedGeometry::SRotatedGeometry(const Geometry &g,
				   double srot, BalanceMode mode):
  srotIn(srot), srotOut(0), geom(g)
{
  balanceSrots(mode);
}

SRotatedGeometry::SRotatedGeometry(const SRotatedGeometry &rhs):
  srotIn(rhs.srotIn),srotOut(rhs.srotOut),geom(rhs.geom)
{}


SRotatedGeometry::~SRotatedGeometry()
{}


double SRotatedGeometry::getArcLength() const
{
  return geom.getArcLength();
}


double SRotatedGeometry::getElementLength() const
{
  return geom.getElementLength();
}


double SRotatedGeometry::getSrotIn()  const
{
  return srotIn;
}


double SRotatedGeometry::getSrotOut() const
{
  return srotOut;
}


void SRotatedGeometry::setSrotIn(double phi)
{
  srotIn=phi;
}


void SRotatedGeometry::setSrotOut(double phi)
{
  srotOut=phi;
}

void SRotatedGeometry::balanceSrots(BalanceMode mode)
{
//   switch(mode) {
//   tilt:
//     srotOut = -srotIn;
//     break;
//   balanceX:
//     // to be implemented
//     break;
//   }
  if(mode==tilt)
    srotOut = -srotIn;
  else {
    // balanceX to be implemented
  }
}

      
double SRotatedGeometry::getOrigin() const
{
  return geom.getOrigin();
}


double SRotatedGeometry::getEntrance() const
{
  return geom.getEntrance();
}


double SRotatedGeometry::getExit() const
{
  return geom.getExit();
}


Euclid3D SRotatedGeometry::getTransform(double fromS, double toS) const
{
  Euclid3D t = geom.getTransform(fromS,toS);

  if(fromS == geom.getEntrance())
    t = t*Euclid3D::ZRotation(srotIn);
  else if(fromS == geom.getExit())
    t = t*Euclid3D::ZRotation(-srotOut);  

  if(toS == geom.getEntrance())
    t *= Euclid3D::ZRotation(-srotIn);
  else if(toS == geom.getExit())
    t *= Euclid3D::ZRotation(srotOut);

  return t;
}


Euclid3D SRotatedGeometry::getTotalTransform() const
{
  Euclid3D tin = Euclid3D::ZRotation(srotIn);
  Euclid3D tout = Euclid3D::ZRotation(srotOut);
  Euclid3D t = geom.getTotalTransform();

  return tout*t*tin;
}


Euclid3D SRotatedGeometry::getTransform(double s) const
{
  return getTransform(0,s);
}


Euclid3D SRotatedGeometry::getExitFrame() const
{
  Euclid3D t = geom.getExitFrame();
  return Euclid3D::ZRotation(srotOut)*t;
}


Euclid3D SRotatedGeometry::getEntranceFrame() const
{
  Euclid3D t = geom.getEntranceFrame();
  return t*Euclid3D::ZRotation(srotIn);
}


Euclid3D SRotatedGeometry::getEntrancePatch() const
{
  return Euclid3D::ZRotation(srotIn);
}


Euclid3D SRotatedGeometry::getExitPatch() const
{
  return Euclid3D::ZRotation(srotOut);
}
