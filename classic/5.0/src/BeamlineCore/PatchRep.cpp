// ------------------------------------------------------------------------
// $RCSfile: PatchRep.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: PatchRep
//   Defines a concrete representation for a geometry patch.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/PatchRep.h"
#include "AbsBeamline/ElementImage.h"
#include "Channels/IndirectChannel.h"


// Attribute access table.
// ------------------------------------------------------------------------

namespace {
  struct Entry {
    const char *name;
    double (PatchRep::*get)() const;
    void (PatchRep::*set)(double);
  };
  
  const Entry entries[] = {
    {
      "X",
      &PatchRep::getX,
      &PatchRep::setX
    },
    {
      "Y",
      &PatchRep::getY,
      &PatchRep::setY
    },
    {
      "Z",
      &PatchRep::getZ,
      &PatchRep::setZ
    },
    {
      "VX",
      &PatchRep::getVX,
      &PatchRep::setVX
    },
    {
      "VY",
      &PatchRep::getVY,
      &PatchRep::setVY
    },
    {
      "VZ",
      &PatchRep::getVZ,
      &PatchRep::setVZ
    },
    { 0, 0, 0 }
  };
}


// Class PatchRep
// ------------------------------------------------------------------------

PatchRep::PatchRep():
  Patch(), geometry(), patch()
{}


PatchRep::PatchRep(const PatchRep &rhs):
  Patch(rhs), geometry(), patch(rhs.patch)
{}


PatchRep::PatchRep(const string &name):
  Patch(name), geometry(), patch()
{}


PatchRep::~PatchRep()
{}


ElementBase *PatchRep::clone() const
{
  return new PatchRep(*this);
}


Channel *PatchRep::getChannel(const string &aKey)
{
  for (const Entry *entry = entries; entry->name != 0; ++entry) {
    if (aKey == entry->name) {
      return new IndirectChannel<PatchRep>(*this, entry->get, entry->set);
    }
  }
  
  return ElementBase::getChannel(aKey);
}


NullField &PatchRep::getField()
{
  return field;
}

const NullField &PatchRep::getField() const
{
  return field;
}


NullGeometry &PatchRep::getGeometry()
{
  return geometry;
}


const NullGeometry &PatchRep::getGeometry() const
{
  return geometry;
}


const Euclid3D &PatchRep::getPatch() const
{
  return patch;
}


void PatchRep::setPatch(const Euclid3D &p)
{
  patch = p;
}


void PatchRep::setPatch(double x, double y, double z,
			double vx, double vy, double vz)
{
  patch = Euclid3D(x, y, z, vx, vy, vz);
}


ElementImage *PatchRep::getImage() const
{
  ElementImage *image = ElementBase::getImage();

  for (const Entry *entry = entries; entry->name != 0; ++entry) {
    image->setAttribute(entry->name, (this->*(entry->get))());
  }

  return image;
}


double PatchRep::getX() const
{
  return patch.getX();
}


double PatchRep::getY() const
{
  return patch.getY();
}


double PatchRep::getZ() const
{
  return patch.getZ();
}


double PatchRep::getVX() const
{
  double vx, vy, vz;
  patch.getRotation().getAxis(vx, vy, vz);
  return vx;
}


double PatchRep::getVY() const
{
  double vx, vy, vz;
  patch.getRotation().getAxis(vx, vy, vz);
  return vy;
}


double PatchRep::getVZ() const
{
  double vx, vy, vz;
  patch.getRotation().getAxis(vx, vy, vz);
  return vz;
}


void PatchRep::setX(double x)
{
  patch.setX(x);
}


void PatchRep::setY(double y)
{
  patch.setY(y);
}


void PatchRep::setZ(double z)
{
  patch.setZ(z);
}


void PatchRep::setVX(double v)
{
  double vx, vy, vz;
  patch.getRotation().getAxis(vx, vy, vz);
  vx = v;
  patch.setRotation(Rotation3D(vx, vy, vz));
}


void PatchRep::setVY(double v)
{
  double vx, vy, vz;
  patch.getRotation().getAxis(vx, vy, vz);
  vy = v;
  patch.setRotation(Rotation3D(vx, vy, vz));
}


void PatchRep::setVZ(double v)
{
  double vx, vy, vz;
  patch.getRotation().getAxis(vx, vy, vz);
  vz = v;
  patch.setRotation(Rotation3D(vx, vy, vz));
}
