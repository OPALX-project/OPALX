// ------------------------------------------------------------------------
// $RCSfile: YCorrectorRep.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: YCorrectorRep
//   Defines a concrete representation for a vertical orbit corrector.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/YCorrectorRep.h"
#include "AbsBeamline/ElementImage.h"
#include "Channels/IndirectChannel.h"


// Attribute access table.
// ------------------------------------------------------------------------

namespace {
  struct Entry {
    const char *name;
    double (YCorrectorRep::*get)() const;
    void (YCorrectorRep::*set)(double);
  };
  
  const Entry entries[] = {
    {
      "L",
      &YCorrectorRep::getElementLength,
      &YCorrectorRep::setElementLength
    },
    {
      "BX",
      &YCorrectorRep::getBx,
      &YCorrectorRep::setBx
    },
    { 0, 0, 0 }
  };
}


// Class YCorrectorRep
// ------------------------------------------------------------------------

YCorrectorRep::YCorrectorRep():
  CorrectorRep()
{}


YCorrectorRep::YCorrectorRep(const YCorrectorRep &right):
  CorrectorRep(right)
{}


YCorrectorRep::YCorrectorRep(const string &name):
  CorrectorRep(name)
{}


YCorrectorRep::~YCorrectorRep()
{}


ElementBase *YCorrectorRep::clone() const
{
  return new YCorrectorRep(*this);
}


Channel *YCorrectorRep::getChannel(const string &aKey)
{
  for (const Entry *entry = entries; entry->name != 0; ++entry) {
    if (aKey == entry->name) {
      return new IndirectChannel<YCorrectorRep>(*this, entry->get, entry->set);
    }
  }

  return ElementBase::getChannel(aKey);
}


ElementImage *YCorrectorRep::getImage() const
{
  ElementImage *image = ElementBase::getImage();

  for (const Entry *entry = entries; entry->name != 0; ++entry) {
    image->setAttribute(entry->name, (this->*(entry->get))());
  }

  return image;
}


Corrector::Plane YCorrectorRep::getPlane() const
{
  return Y;
}


double YCorrectorRep::getBy() const
{
  return 0.0;
}


void YCorrectorRep::setBy(double)
{
  // Do nothing.
}

