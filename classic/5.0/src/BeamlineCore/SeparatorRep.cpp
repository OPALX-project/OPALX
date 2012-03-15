// ------------------------------------------------------------------------
// $RCSfile: SeparatorRep.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: SeparatorRep
//   Defines the abstract interface for an electrostatic separator.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/SeparatorRep.h"
#include "AbsBeamline/ElementImage.h"
#include "Channels/IndirectChannel.h"


// Attribute access table.
// ------------------------------------------------------------------------

namespace {
  struct Entry {
    const char *name;
    double (SeparatorRep::*get)() const;
    void (SeparatorRep::*set)(double);
  };
  
  const Entry entries[] = {
    {
      "L",
      &SeparatorRep::getElementLength,
      &SeparatorRep::setElementLength
    },
    {
      "EX",
      &SeparatorRep::getEx,
      &SeparatorRep::setEx
    },
    {
      "EY",
      &SeparatorRep::getEy,
      &SeparatorRep::setEy
    },
    { NULL, NULL, NULL }
  };
}


// Class SeparatorRep
// ------------------------------------------------------------------------

SeparatorRep::SeparatorRep():
  Separator(), 
  geometry(),
  field()
{}


SeparatorRep::SeparatorRep(const SeparatorRep &right):
  Separator(right),
  geometry(right.geometry),
  field(right.field)
{}


SeparatorRep::SeparatorRep(const string &name):
  Separator(name),
  geometry(),
  field()
{}


SeparatorRep::~SeparatorRep()
{}


ElementBase *SeparatorRep::clone() const
{
  return new SeparatorRep(*this);
}


Channel *SeparatorRep::getChannel(const string &aKey)
{
  for (const Entry *table = entries; table->name != 0; ++table) {
    if (aKey == table->name) {
      return new IndirectChannel<SeparatorRep>(*this, table->get, table->set);
    }
  }

  return ElementBase::getChannel(aKey);
}


EDipoleField &SeparatorRep::getField()
{
  return field;
}

const EDipoleField &SeparatorRep::getField() const
{
  return field;
}


StraightGeometry &SeparatorRep::getGeometry()
{
  return geometry;
}

const StraightGeometry &SeparatorRep::getGeometry() const
{
  return geometry;
}


ElementImage *SeparatorRep::getImage() const
{
  ElementImage *image = ElementBase::getImage();

  for (const Entry *table = entries; table->name != 0; ++table) {
    image->setAttribute(table->name, (this->*(table->get))());
  }

  return image;
}


double SeparatorRep::getEx() const
{
  return field.getEx();
}


double SeparatorRep::getEy() const
{
  return field.getEy();
}


void SeparatorRep::setEx(double value)
{
  field.setEx(value);
}


void SeparatorRep::setEy(double value)
{
  field.setEy(value);
}
