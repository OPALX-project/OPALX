// ------------------------------------------------------------------------
// $RCSfile: BeamBeam.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: BeamBeam
//   Defines the abstract interface for a beam-beam interaction.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/BeamBeam.h"
#include "AbsBeamline/BeamlineVisitor.h"


// Class BeamBeam
// ------------------------------------------------------------------------

BeamBeam::BeamBeam():
  Component()
{}


BeamBeam::BeamBeam(const BeamBeam &right):
  Component(right)
{}


BeamBeam::BeamBeam(const string &name):
  Component(name)
{}


BeamBeam::~BeamBeam()
{}


void BeamBeam::accept(BeamlineVisitor &visitor) const
{
  visitor.visitBeamBeam(*this);
}
