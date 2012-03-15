// ------------------------------------------------------------------------
// $RCSfile: Monitor.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Monitor
//   Defines the abstract interface for a beam position monitor.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Monitor.h"
#include "AbsBeamline/BeamlineVisitor.h"


// Class Monitor
// ------------------------------------------------------------------------

Monitor::Monitor():
  Component()
{}


Monitor::Monitor(const Monitor &right):
  Component(right)
{}


Monitor::Monitor(const string &name):
  Component(name)
{}


Monitor::~Monitor()
{}


void Monitor::accept(BeamlineVisitor &visitor) const
{
  visitor.visitMonitor(*this);
}

