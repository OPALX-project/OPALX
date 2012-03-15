// ------------------------------------------------------------------------
// $RCSfile: CyclotronRep.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: CyclotronRep
//   Defines a concrete representation for a sector (curved) bend.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 18:57:53 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/CyclotronRep.h"
#include "AbsBeamline/ElementImage.h"
#include "Channels/IndexedChannel.h"
#include "Channels/IndirectChannel.h"
#include "ComponentWrappers/CyclotronWrapper.h"
#include <cctype>

// Attribute access table.
// ------------------------------------------------------------------------

namespace {
  struct Entry {
    const char *name;
    double (CyclotronRep::*get)() const;
    void (CyclotronRep::*set)(double);
  };
  
  static const Entry entries[] = {
    {
      "RINIT",
      //      &CyclotronRep::getRadius,
      //&CyclotronRep::setRadius
    },
    { 0, 0, 0 }
  };
}


// Class CyclotronRep
// ------------------------------------------------------------------------

CyclotronRep::CyclotronRep():
  Cyclotron(),
  geometry(0.0, 0.0),
  field()
{
  rInit = 0.0;
}


CyclotronRep::CyclotronRep(const CyclotronRep &rhs):
  Cyclotron(rhs),
  geometry(rhs.geometry),
  field(rhs.field)
{
  rInit = rhs.rInit;
}


CyclotronRep::CyclotronRep(const string &name):
  Cyclotron(name),
  geometry(0.0, 0.0),
  field()
{
  rInit = 0.0;
}


CyclotronRep::~CyclotronRep()
{}


ElementBase *CyclotronRep::clone() const
{
  return new CyclotronRep(*this);
}


Channel *CyclotronRep::getChannel(const string &aKey)
{
  if (aKey[0] == 'a'  ||  aKey[0] == 'b') {
    int n = 0;
    
    for (string::size_type k = 1; k < aKey.length(); k++) {
      if (isdigit(aKey[k])) {
	n = 10 * n + aKey[k] - '0';
      } else {
	return 0;
      }
    }
  } else {
    for (const Entry *table = entries; table->name != 0; ++table) {
      if (aKey == table->name) {
	return new IndirectChannel<CyclotronRep>(*this, table->get, table->set);
      }
    }
    
    return ElementBase::getChannel(aKey);
  }
}


const string &CyclotronRep::getType() const
{
  static const string type("CyclotronRep");
  return type;
}



double CyclotronRep::getSlices() const
{
  return slices;
}

double CyclotronRep::getStepsize() const
{
  return stepsize;
}

void CyclotronRep::setSlices(double sl)
{
  slices = sl;
}

void CyclotronRep::setStepsize(double ds)
{
  stepsize = ds;
}

/*
void CyclotronRep::setRadius(double r)
{
  rInit = r;
}

double CyclotronRep::getRadius() const
{
  return rInit ;
}
*/

PlanarArcGeometry &CyclotronRep::getGeometry()
{
  return geometry;
}

const PlanarArcGeometry &CyclotronRep::getGeometry() const
{
  return geometry;
}

BMultipoleField &CyclotronRep::getField()
{
  return field;
}

const BMultipoleField &CyclotronRep::getField() const
{
  return field;
}

void CyclotronRep::setField(const BMultipoleField &f)
{
  field = f;
}



ElementBase *CyclotronRep::makeFieldWrapper()
{
  ElementBase *wrap = new CyclotronWrapper(this);
  wrap->setName(getName());
  return wrap;
}
