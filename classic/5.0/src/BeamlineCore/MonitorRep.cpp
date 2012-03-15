// ------------------------------------------------------------------------
// $RCSfile: MonitorRep.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MonitorRep
//   Defines a concrete representation for an orbit monitor acting on
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

#include "BeamlineCore/MonitorRep.h"
#include "AbsBeamline/ElementImage.h"
#include "Channels/IndirectChannel.h"


// Attribute access table.
// ------------------------------------------------------------------------

namespace {
  struct Entry {
    const char *name;
    double (MonitorRep::*get)() const;
    void (MonitorRep::*set)(double);
  };
  
  const Entry entries[] = {
    {
      "L",
      &MonitorRep::getElementLength,
      &MonitorRep::setElementLength
    },
    { 0, 0, 0 }
  };
}


// Class MonitorRep
// ------------------------------------------------------------------------

MonitorRep::MonitorRep():
  Monitor(), field(), geometry(), active(true)
{}


MonitorRep::MonitorRep(const MonitorRep &right):
  Monitor(right), field(), geometry(right.geometry), active(true)
{}


MonitorRep::MonitorRep(const string &name):
  Monitor(name), field(), geometry(), active(true)
{}


MonitorRep::~MonitorRep()
{}


ElementBase *MonitorRep::clone() const
{
  return new MonitorRep(*this);
}


Channel *MonitorRep::getChannel(const string &aKey)
{
  for (const Entry *entry = entries; entry->name != 0; ++entry) {
    if (aKey == entry->name) {
      return new IndirectChannel<MonitorRep>(*this, entry->get, entry->set);
    }
  }
  
  return ElementBase::getChannel(aKey);
}


NullField &MonitorRep::getField()
{
  return field;
}

const NullField &MonitorRep::getField() const
{
  return field;
}


StraightGeometry &MonitorRep::getGeometry()
{
  return geometry;
}

const StraightGeometry &MonitorRep::getGeometry() const
{
  return geometry;
}


ElementImage *MonitorRep::getImage() const
{
  ElementImage *image = ElementBase::getImage();

  for (const Entry *entry = entries; entry->name != 0; ++entry) {
    image->setAttribute(entry->name, (this->*(entry->get))());
  }

  return image;
}


Monitor::Plane MonitorRep::getPlane() const
{
  return active ? XY : OFF;
}


const string &MonitorRep::getType() const
{
  static const string type("MonitorRep");
  return type;
}


void MonitorRep::setActive(bool flag)
{
  active = flag;
}

