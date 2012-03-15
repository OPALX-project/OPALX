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
    LINE,     // The name of lattice to be tracked.
    BEAM,     // The name of beam to be used.
    // THE INTEGRATION TIMESTEP IN SEC
    DT,
    // THE MAXIMUM TIMESTEPS WE INTEGRATE
    MAXSTEPS,
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
    ("MAXSTEPS","THE MAXIMUM NUMBER OF INTEGRATION DT",10);
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

int TrackCmd::getMAXSTEPS() const
{
  return (int) Attributes::getReal(itsAttr[MAXSTEPS]);
}



void TrackCmd::execute()
{
  // Find BeamSequence and Beam definitions.
  BeamSequence *use = BeamSequence::find(Attributes::getString(itsAttr[LINE]));
  Beam *beam = Beam::find(Attributes::getString(itsAttr[BEAM]));

  double dt = getDT();
  int    maxsteps = getMAXSTEPS();


  // Execute track block.
  Track::block = new Track(use, beam->getReference(), dt, maxsteps);
  Track::block->parser.run();

  // Clean up.
  delete Track::block;
  Track::block = 0;
}
