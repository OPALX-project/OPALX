// ------------------------------------------------------------------------
// $RCSfile: CorrectorRep.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: CorrectorRep
//   Defines a concrete representation for an orbit corrector acting on
//   both planes.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/CorrectorRep.h"
#include "AbsBeamline/ElementImage.h"
#include "Channels/IndirectChannel.h"
#include "ComponentWrappers/CorrectorWrapper.h"


// Attribute access table.
// ------------------------------------------------------------------------

namespace {
  struct Entry {
    const char *name;
    double (CorrectorRep::*get)() const;
    void (CorrectorRep::*set)(double);
  };
  
  const Entry entries[] = {
    {
      "L",
      &CorrectorRep::getElementLength,
      &CorrectorRep::setElementLength
    },
    {
      "BY",
      &CorrectorRep::getBy,
      &CorrectorRep::setBy
    },
    {
      "BX",
      &CorrectorRep::getBx,
      &CorrectorRep::setBx
    },
    { 0, 0, 0 }
  };
}


// Class CorrectorRep
// ------------------------------------------------------------------------

CorrectorRep::CorrectorRep():
  Corrector(), geometry(), field(), active(true)
{}
 

CorrectorRep::CorrectorRep(const CorrectorRep &right):
  Corrector(right), geometry(right.geometry), field(right.field), active(true)
{}


CorrectorRep::CorrectorRep(const string &name):
  Corrector(name), geometry(), field(), active(true)
{}


CorrectorRep::~CorrectorRep()
{}


ElementBase *CorrectorRep::clone() const
{
  return new CorrectorRep(*this);
}


Channel *CorrectorRep::getChannel(const string &aKey)
{
  for (const Entry *table = entries; table->name != 0; ++table) {
    if (aKey == table->name) {
      return new IndirectChannel<CorrectorRep>(*this, table->get, table->set);
    }
  }

  return ElementBase::getChannel(aKey);
}


StraightGeometry &CorrectorRep::getGeometry()
{
  return geometry;
}


const StraightGeometry &CorrectorRep::getGeometry() const
{
  return geometry;
}


ElementImage *CorrectorRep::getImage() const
{
  ElementImage *image = ElementBase::getImage();

  for (const Entry *table = entries; table->name != 0; ++table) {
    image->setAttribute(table->name, (this->*(table->get))());
  }

  return image;
}


Corrector::Plane CorrectorRep::getPlane() const
{
  return active ? XY : OFF;
}


const string &CorrectorRep::getType() const
{
  static const string type("CorrectorRep");
  return type;
}


double CorrectorRep::getBx() const
{
  return field.getBx();
}


double CorrectorRep::getBy() const
{
  return field.getBy();
}


BDipoleField &CorrectorRep::getField()
{
  return field;
}


const BDipoleField &CorrectorRep::getField() const
{
  return field;
}


void CorrectorRep::setBx(double Bx)
{
  field.setBx(Bx);
}


void CorrectorRep::setBy(double By)
{
  field.setBy(By);
}


ElementBase *CorrectorRep::makeFieldWrapper()
{
  ElementBase *wrap = new CorrectorWrapper(this);
  wrap->setName(getName());
  return wrap;
}


void CorrectorRep::setActive(bool flag)
{
  active = flag;
}
