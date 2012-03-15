// ------------------------------------------------------------------------
// $RCSfile: TrackStart.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TrackStart
//   The class for the OPAL START command.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:47 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Track/TrackStart.h"
#include "Attributes/Attributes.h"
#include "Track/Track.h"


// Class TrackStart
// ------------------------------------------------------------------------

// The attributes of class TrackStart.
namespace {
  enum {
    X,  // The initial horizontal   position in m.
    Y,  // The initial vertical     position in m.
    T,  // The initial longitudinal position in m.
    PX, // The initial horizontal   momentum in rad.
    PY, // The initial vertical     momentum in rad.
    PT, // The initial longitudinal momentum in rad.
    SIZE
  };
}


TrackStart::TrackStart():
  Action(SIZE, "START",
	 "The \"START\" sub-command defines one particle for tracking "
	 "through the given lattice.")
{
  itsAttr[X]  = Attributes::makeReal
    ("X",  "Initial horizontal position in m");
  itsAttr[PX] = Attributes::makeReal
    ("PX", "Initial horizontal momentum in 1");
  itsAttr[Y]  = Attributes::makeReal
    ("Y",  "Initial horizontal position in m");
  itsAttr[PY] = Attributes::makeReal
    ("PY", "Initial horizontal momentum in 1");
  itsAttr[T]  = Attributes::makeReal
    ("T",  "Initial horizontal position in m");
  itsAttr[PT] = Attributes::makeReal
    ("PT", "Initial horizontal momentum in 1");
}


TrackStart::TrackStart(const string &name, TrackStart *parent):
  Action(name, parent)
{}


TrackStart::~TrackStart()
{}


TrackStart *TrackStart::clone(const string &name)
{
  return new TrackStart(name, this);
}


void TrackStart::execute()
{
  double x  = Attributes::getReal(itsAttr[X]);
  double y  = Attributes::getReal(itsAttr[Y]);
  double t  = Attributes::getReal(itsAttr[T]);
  double px = Attributes::getReal(itsAttr[PX]);
  double py = Attributes::getReal(itsAttr[PY]);
  double pt = Attributes::getReal(itsAttr[PT]);
  Track::block->bunch->push_back(Particle(x, px, y, py, t, pt));
}
