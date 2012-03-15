// ------------------------------------------------------------------------
// $RCSfile: TrackCmd.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TrackCmd
//   The class for the OPAL TRACK command.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:47 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Track/TrackCmd.h"
#include "AbstractObjects/BeamSequence.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Structure/Beam.h"
#include "Track/Track.h"
#include "Track/TrackParser.h"


// Class Track
// ------------------------------------------------------------------------

namespace {

  // The attributes of class TrackRun.
  enum {
    LINE,         // The name of lattice to be tracked.
    BEAM,         // The name of beam to be used.
    DT,           // The integration timestep in second.
    MAXSTEPS,     // The maximum timesteps we integrate
    ZSTOP,        // Defines a z-location [m], after which the simulation stops when the last particles passes 
    STEPSPERTURN, // Return the timsteps per revolution period. ONLY available for OPAL-cycl.
    TIMEINTEGRATOR, // the name of time integrator
	NNB, // Number of neighbouring bunches in OPAL-cycl
    SIZE
  };
}

TrackCmd::TrackCmd():
  Action(SIZE, "TRACK",
	 "The \"TRACK\" command initiates tracking.")
{
  itsAttr[LINE] = Attributes::makeString
     ("LINE", "Name of lattice to be tracked");
  itsAttr[BEAM] = Attributes::makeString
     ("BEAM", "Name of beam to be used", "UNNAMED_BEAM");
  itsAttr[DT] = Attributes::makeReal
     ("DT","THE INTEGRATION TIMESTEP IN SECONDS",1e-12);
  itsAttr[MAXSTEPS] = Attributes::makeReal
     ("MAXSTEPS","THE MAXIMUM NUMBER OF INTEGRATION STEPS DT, should be larger ZSTOP/(beta*c average)",10);
  itsAttr[STEPSPERTURN] = Attributes::makeReal
     ("STEPSPERTURN","THE TIME STEPS PER REVOLUTION PERIOD, ONLY FOR OPAL-CYCL",720);
  itsAttr[ZSTOP] = Attributes::makeReal
     ("ZSTOP","Defines a z-location [m], after which the simulation stops when the last particles passes",1000000.0); 
  itsAttr[TIMEINTEGRATOR] = Attributes::makeString
    ("TIMEINTEGRATOR", "Name of time integrator to be used", "RK-4");
  itsAttr[NNB] = Attributes::makeReal
	("NNB","Number of neighbouring bunches in OPAL-cycl",0.0);
}

TrackCmd::TrackCmd(const string &name, TrackCmd *parent):
  Action(name, parent)
{}


TrackCmd::~TrackCmd()
{}


TrackCmd *TrackCmd::clone(const string &name)
{
  return new TrackCmd(name, this);
}

double TrackCmd::getDT() const
{
  return Attributes::getReal(itsAttr[DT]);
}

double TrackCmd::getZSTOP() const
{
  return Attributes::getReal(itsAttr[ZSTOP]);
}

int TrackCmd::getMAXSTEPS() const
{
  return (int) Attributes::getReal(itsAttr[MAXSTEPS]);
}

int TrackCmd::getSTEPSPERTURN() const
{
  return (int) Attributes::getReal(itsAttr[STEPSPERTURN]);
}

int TrackCmd::getNNB() const
{
  return (int) Attributes::getReal(itsAttr[NNB]);
}

// return int type rathor than string to improve the speed 
int TrackCmd::getTIMEINTEGRATOR() const
{
  string name = Attributes::getString(itsAttr[TIMEINTEGRATOR]);
  int  nameID;
  if (name == string("RK-4")) 
    nameID =  0;
  else if(name == string("LF-2"))
    nameID =  1;     
  else
    nameID = -1;

  return nameID;
}

void TrackCmd::execute()
{
  // Find BeamSequence and Beam definitions.
  BeamSequence *use = BeamSequence::find(Attributes::getString(itsAttr[LINE]));
  Beam *beam = Beam::find(Attributes::getString(itsAttr[BEAM]));

  double dt = getDT();
  int    maxsteps = getMAXSTEPS();
  int    stepsperturn = getSTEPSPERTURN();
  double zstop = getZSTOP();
  int timeintegrator = getTIMEINTEGRATOR();
  
  // Execute track block.
  Track::block = new Track(use, beam->getReference(), dt, maxsteps, stepsperturn, zstop, timeintegrator);
  Track::block->parser.run();

  // Clean up.
  delete Track::block;
  Track::block = 0;
}
