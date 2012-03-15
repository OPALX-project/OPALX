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
  lengthUnit_m(1.0),
  fast_m(false)
{ 
  setElType(isSolenoid);
}


Solenoid::Solenoid(const Solenoid &right):
  Component(right),
  filename_m(right.filename_m),
  scale_m(right.scale_m),
  fast_m(right.fast_m),
  lengthUnit_m(right.lengthUnit_m),
  myFieldmap(right.myFieldmap)
{
  setElType(isSolenoid);
}


Solenoid::Solenoid(const string &name):
  Component(name)
{
  setElType(isSolenoid);
}


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

bool Solenoid::apply(const int &i, const double &t, double E[], double B[])
{
  Vector_t Ev(0,0,0), Bv(0,0,0);
  if (apply(RefPartBunch_m->R[i],t,Ev,Bv)) return true;
      
  E[0] = Ev(0); E[1] = Ev(1); E[2] = Ev(2);
  B[0] = Bv(0); B[1] = Bv(1); B[2] = Bv(2);
      
  return false;
}

bool Solenoid::apply(const int &i, const double &t, Vector_t &E, Vector_t &B)
{
  const Vector_t tmpR(RefPartBunch_m->R[i](0) - dx_m, RefPartBunch_m->R[i](1) - dy_m, RefPartBunch_m->R[i](2) - startField_m - ds_m);
  Vector_t tmpE(0.0,0.0,0.0), tmpB(0.0,0.0,0.0);

  const bool out_of_bounds = myFieldmap->getFieldstrength(tmpR,tmpE,tmpB);
  B += scale_m * tmpB;

  return out_of_bounds;
}
  
bool Solenoid::apply( const Vector_t &R, const  double &t, Vector_t &E, Vector_t &B)
{
  const Vector_t tmpR(R(0) - ds_m, R(1) - ds_m, R(2) - startField_m - ds_m);
  Vector_t tmpE(0.0,0.0,0.0), tmpB(0.0,0.0,0.0);

  const bool out_of_bounds = myFieldmap->getFieldstrength(tmpR,tmpE,tmpB);
  B += scale_m * tmpB;

  return out_of_bounds;
}

void Solenoid::initialise(const PartBunch *bunch, double &startField, double &endField, const double &scaleFactor)
{

  Inform msg("Solenoid ");

  RefPartBunch_m = bunch;

  msg << getName() << " using file ";
  myFieldmap = Fieldmap::getFieldmap(filename_m,fast_m);
  if (myFieldmap != NULL)
    {
      myFieldmap->getInfo(&msg);
      if (dx_m != 0.0 || dy_m != 0.0 || ds_m != 0.0)
        msg << "misaligned by dx = " << dx_m << ", dy = " << dy_m << ", dz = " << ds_m << endl;

      double zBegin = 0.0, zEnd = 0.0, rBegin = 0.0, rEnd = 0.0;
      myFieldmap->getFieldDimensions(zBegin, zEnd, rBegin, rEnd);

      endField = startField + zEnd;
      startField += zBegin;
      startField_m = startField;
      endField_m = endField;

    }
  else
    {
      endField = startField;
    }
}

void Solenoid::finalise()
{}

void Solenoid::rescaleFieldMap(const double &scaleFactor) 
{
  startField_m *= scaleFactor/lengthUnit_m;
  dx_m *= scaleFactor/lengthUnit_m;
  dy_m *= scaleFactor/lengthUnit_m;
  ds_m *= scaleFactor/lengthUnit_m;
  myFieldmap->rescale(scaleFactor);
  lengthUnit_m = scaleFactor;
}

bool Solenoid::bends() const
{
  return false;
}


void Solenoid::goOnline()
{
  Fieldmap::readMap(filename_m);
  online_m = true;
}

void Solenoid::goOffline()
{
  Fieldmap::freeMap(filename_m);
  online_m = false;
}

void Solenoid::setKS(double ks)
{
  scale_m = ks;
}

void Solenoid::getDimensions(double &zBegin, double &zEnd) const
{
  zBegin = startField_m;
  zEnd = endField_m;
}
