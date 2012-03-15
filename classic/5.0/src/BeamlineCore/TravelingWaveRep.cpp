// ------------------------------------------------------------------------
// $RCSfile: RFCavityRep.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RFCavityRep
//   Defines a representation for a RF cavity.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/TravelingWaveRep.h"
#include "AbsBeamline/ElementImage.h"
#include "Channels/IndirectChannel.h"


// Attribute access table.
// ------------------------------------------------------------------------

bool TravelingWaveRep::ignoreCavities = false;

namespace {
  struct Entry {
    const char *name;
    double (TravelingWaveRep::*get)() const;
    void (TravelingWaveRep::*set)(double);
  };
  
  static const Entry entries[] = {
    {
      "L",
      &TravelingWaveRep::getElementLength,
      &TravelingWaveRep::setElementLength
    },
    {
      "AMPLITUDE",
      &TravelingWaveRep::getAmplitude,
      &TravelingWaveRep::setAmplitude
    },
    {
      "FREQUENCY",
      &TravelingWaveRep::getFrequency,
      &TravelingWaveRep::setFrequency
    },
    {
      "PHASE",
      &TravelingWaveRep::getPhase,
      &TravelingWaveRep::setPhase
    }
  };
}


// Class TravelingWaveRep
// ------------------------------------------------------------------------

TravelingWaveRep::TravelingWaveRep():
  TravelingWave()
{}


TravelingWaveRep::TravelingWaveRep(const TravelingWaveRep &right):
  TravelingWave(right)
{}


TravelingWaveRep::TravelingWaveRep(const string &name):
  TravelingWave(name)
{}


TravelingWaveRep::~TravelingWaveRep()
{}


ElementBase *TravelingWaveRep::clone() const
{
  return new TravelingWaveRep(*this);
}


Channel *TravelingWaveRep::getChannel(const string &aKey)
{
  for (const Entry *entry = entries; entry->name != 0; ++entry) {
    if (aKey == entry->name) {
      return new IndirectChannel<TravelingWaveRep>(*this, entry->get, entry->set);
    }
  }

  return ElementBase::getChannel(aKey);
}


AcceleratingField &TravelingWaveRep::getField()
{
  return field;
}

const AcceleratingField &TravelingWaveRep::getField() const
{
  return field;
}


StraightGeometry &TravelingWaveRep::getGeometry()
{
  return geometry;
}

const StraightGeometry &TravelingWaveRep::getGeometry() const
{
  return geometry;
}


ElementImage *TravelingWaveRep::getImage() const
{
  ElementImage *image = ElementBase::getImage();

  for (const Entry *entry = entries; entry->name != 0; ++entry) {
    image->setAttribute(entry->name, (this->*(entry->get))());
  }

  return image;
}


const string &TravelingWaveRep::getType() const
{
  static const string type("TravelingWaveRep");
  return type;
}


double TravelingWaveRep::getAmplitude() const
{
  return ignoreCavities ? 0.0 : field.getEz();
}


double TravelingWaveRep::getFrequency() const
{
  return field.getFrequency();
}


double TravelingWaveRep::getPhase() const
{
  return field.getPhase();
}


void TravelingWaveRep::setAmplitude(double amplitude)
{
  field.setEz(amplitude);
}


void TravelingWaveRep::setFrequency(double frequency)
{
  field.setFrequency(frequency);
}


void TravelingWaveRep::setPhase(double phase)
{
  field.setPhase(phase);
}

void TravelingWaveRep::setIgnore(bool ignore)
{
  ignoreCavities = ignore;
}
