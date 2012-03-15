// ------------------------------------------------------------------------
// $RCSfile: XCorrectorRep.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: XCorrectorRep
//   Defines a concrete representation for a horizontal orbit corrector.
//
// ------------------------------------------------------------------------
// 02-Nov-98, Chris Iselin, CERN
//   Original release.
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// Revision History:
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/XCorrectorRep.h"
#include "AbsBeamline/ElementImage.h"
#include "Channels/IndirectChannel.h"


// Attribute access table.
// ------------------------------------------------------------------------

namespace {
  struct Entry {
    const char *name;
    double (XCorrectorRep::*get)() const;
    void (XCorrectorRep::*set)(double);
  };
  
  const Entry entries[] = {
    {
      "L",
      &XCorrectorRep::getElementLength,
      &XCorrectorRep::setElementLength
    },
    {
      "BY",
      &XCorrectorRep::getBy,
      &XCorrectorRep::setBy
    },
    { 0, 0, 0 }
  };
}


// Class XCorrectorRep
// ------------------------------------------------------------------------

XCorrectorRep::XCorrectorRep():
  CorrectorRep()
{}


XCorrectorRep::XCorrectorRep(const XCorrectorRep &rhs):
  CorrectorRep(rhs)
{}


XCorrectorRep::XCorrectorRep(const string &name):
  CorrectorRep(name)
{}


XCorrectorRep::~XCorrectorRep()
{}


ElementBase *XCorrectorRep::clone() const
{
  return new XCorrectorRep(*this);
}


Channel *XCorrectorRep::getChannel(const string &aKey)
{
  for (const Entry *entry = entries; entry->name != 0; ++entry) {
    if (aKey == entry->name) {
      return new IndirectChannel<XCorrectorRep>(*this, entry->get, entry->set);
    }
  }

  return ElementBase::getChannel(aKey);
}


ElementImage *XCorrectorRep::getImage() const
{
  ElementImage *image = ElementBase::getImage();

  for (const Entry *entry = entries; entry->name != 0; ++entry) {
    image->setAttribute(entry->name, (this->*(entry->get))());
  }

  return image;
}


Corrector::Plane XCorrectorRep::getPlane() const
{
  return active ? X : OFF;
}


const string &XCorrectorRep::getType() const
{
  static const string type("XCorrectorRep");
  return type;
}


double XCorrectorRep::getBx() const
{
  return 0.0;
}


void XCorrectorRep::setBx(double)
{
  // Do nothing.
}
