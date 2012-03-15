// ------------------------------------------------------------------------
// $RCSfile: Component.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Component
//   An abstract base class which defines the common interface for all
//   CLASSIC components, i.e. beamline members which are not themselves
//   beamlines.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Component.h"
#include "Utilities/LogicalError.h"


// Class Component
// ------------------------------------------------------------------------
//   Represents an arbitrary component in an accelerator.  A component is
//   the basic element in the accelerator model, and can be thought of as
//   acting as a leaf in the Composite pattern.  A Component is associated
//   with an electromagnetic field.

Component::Component():
  ElementBase()
{}


Component::Component(const Component &right):
  ElementBase(right)
{}


Component::Component(const string &name):
  ElementBase(name)
{

}


Component::~Component()
{}


const ElementBase &Component::getDesign() const
{
  return *this;
}

bool Component::getFieldstrength(double R[], double t, double E[], double B[]) const
{

}

bool Component::getFieldstrength(Vector_t R, double t, Vector_t &E, Vector_t &B) const
{

}


bool Component::readFieldMap(double &startField, double &endField, double scaleFactor)

{

}

void Component::rescaleFieldMap(double scaleFactor)
{

}

void Component::trackBunch(PartBunch &, const PartData &, bool, bool) const
{
  throw LogicalError("Component::trackBunch()",
		     "Called for component \"" + getName() + "\".");
}


void Component::
trackMap(FVps<double,6> &, const PartData &, bool, bool) const
{
  throw LogicalError("Component::trackMap()",
		     "Called for component \"" + getName() + "\".");
}

void Component::readFieldMap(double scaleFactor)
{
  
}
