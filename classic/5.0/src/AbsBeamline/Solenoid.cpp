// ------------------------------------------------------------------------
// $RCSfile: Solenoid.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Solenoid
//   Defines the abstract interface for a solenoid magnet.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:32 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Solenoid.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include <iostream>
#include <fstream>


// Class Solenoid
// ------------------------------------------------------------------------

Solenoid::Solenoid():
  Component(),
  dx_m(0.0),
  dy_m(0.0),
  dz_m(0.0),
  lengthUnit_m(1.0),
  fast_m(false)
{ }


Solenoid::Solenoid(const Solenoid &right):
  Component(right),
  filename_m(right.filename_m),
  scale_m(right.scale_m),
  fast_m(right.fast_m),
  dx_m(right.dx_m),
  dy_m(right.dy_m),
  dz_m(right.dz_m),
  lengthUnit_m(right.lengthUnit_m),
  myFieldmap(right.myFieldmap)
{}


Solenoid::Solenoid(const string &name):
  Component(name)
{}


Solenoid::~Solenoid()
{
  Fieldmap::deleteFieldmap(filename_m);
}


void Solenoid::accept(BeamlineVisitor &visitor) const
{
  visitor.visitSolenoid(*this);
}

void Solenoid::setFieldMapFN(string fn)
{
  filename_m = fn;
}

void Solenoid::setFast(bool fast)
{
  fast_m = fast;
}


bool Solenoid::getFast() const
{
  return fast_m;
}

void Solenoid::setMisalignment(double x, double y, double z)
{
  dx_m = x;
  dy_m = y;
  dz_m = z;
}

void Solenoid::getMisalignment(double &x, double &y, double &z) const
{
  x = dx_m;
  y = dy_m;
  z = dz_m;
}

bool Solenoid::readFieldMap(double &startField, double &endField, double scaleFactor)
{

  Inform msg("Solenoid ");
  msg << getName() << " using file ";
  myFieldmap = Fieldmap::getFieldmap(filename_m,fast_m);
  myFieldmap->getInfo(&msg);
  if (dx_m != 0.0 || dy_m != 0.0 || dz_m != 0.0)
    msg << "misaligned by dx = " << dx_m << ", dy = " << dy_m << ", dz = " << dz_m << endl;

  double zBegin = 0.0, zEnd = 0.0, rBegin = 0.0, rEnd = 0.0;
  myFieldmap->getFieldDimensions(zBegin, zEnd, rBegin, rEnd);

  endField = startField + zEnd;
  startField += zBegin;
  startField_m = startField;

  return true;
}

void Solenoid::rescaleFieldMap(double scaleFactor) 
{
  startField_m *= scaleFactor/lengthUnit_m;
  dx_m *= scaleFactor/lengthUnit_m;
  dy_m *= scaleFactor/lengthUnit_m;
  dz_m *= scaleFactor/lengthUnit_m;
  myFieldmap->rescale(scaleFactor);
  lengthUnit_m = scaleFactor;
}

bool Solenoid::getFieldstrength(double R[], double t, double E[], double B[]) const
{
  Vector_t Rv(R[0],R[1],R[2]), Ev(0,0,0), Bv(0,0,0);
  if (getFieldstrength(Rv,t,Ev,Bv)) return true;
  
  E[0] = Ev(0); E[1] = Ev(1); E[2] = Ev(2);
  B[0] = Bv(0); B[1] = Bv(1); B[2] = Bv(2);

  return false;
}

bool Solenoid::getFieldstrength(Vector_t R, double t, Vector_t &E, Vector_t &B) const
{
  Vector_t tmpR(R(0) - dx_m, R(1) - dy_m, R(2) - startField_m - dz_m), 
    tmpE(0.0,0.0,0.0), 
    tmpB(0.0,0.0,0.0);
  bool out_of_bounds = myFieldmap->getFieldstrength(tmpR,tmpE,tmpB);
  B += scale_m * tmpB;
  return out_of_bounds;
}

void Solenoid::setKS(double ks)
{
  scale_m = ks;
}
