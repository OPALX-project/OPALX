// ------------------------------------------------------------------------
// $RCSfile: ProbeRep.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ProbeRep
//   Defines a representation for a septa.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2009/10/07 10:26:06 $
// $Author: bi $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/ProbeRep.h"
#include "AbsBeamline/ElementImage.h"
#include "Channels/IndirectChannel.h"


// Attribute access table.
// ------------------------------------------------------------------------

namespace {
  struct Entry {
    const char *name;
    double (ProbeRep::*get)() const;
    void (ProbeRep::*set)(double);
  };
  
  static const Entry entries[] = {
    {
      "L",
      &ProbeRep::getElementLength,
      &ProbeRep::setElementLength
    },
    { 0, 0, 0 }
  };
}


// Class ProbeRep
// ------------------------------------------------------------------------

ProbeRep::ProbeRep():
  Probe(), field(), geometry(), active(true)
{}


ProbeRep::ProbeRep(const ProbeRep &right):
  Probe(right), field(), geometry(right.geometry), active(true)
{}


ProbeRep::ProbeRep(const string &name):
  Probe(name), field(), geometry(), active(true)
{}


ProbeRep::~ProbeRep()
{}


ElementBase *ProbeRep::clone() const
{
  return new ProbeRep(*this);
}


Channel *ProbeRep::getChannel(const string &aKey)
{
  for (const Entry *entry = entries; entry->name != 0; ++entry) {
    if (aKey == entry->name) {
      return new IndirectChannel<ProbeRep>(*this, entry->get, entry->set);
    }
  }

  return ElementBase::getChannel(aKey);
}

NullField &ProbeRep::getField()
{
  return field;
}

const NullField &ProbeRep::getField() const
{
  return field;
}

StraightGeometry &ProbeRep::getGeometry()
{
  return geometry;
}

const StraightGeometry &ProbeRep::getGeometry() const
{
  return geometry;
}


ElementImage *ProbeRep::getImage() const
{
  ElementImage *image = ElementBase::getImage();

  for (const Entry *entry = entries; entry->name != 0; ++entry) {
    image->setAttribute(entry->name, (this->*(entry->get))());
  }

  return image;
}

void ProbeRep::setActive(bool flag)
{
  active = flag;
}

