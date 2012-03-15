// ------------------------------------------------------------------------
// $RCSfile: Beamline.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Beamline
//   Represents an abstract sequence of elements.
//
// ------------------------------------------------------------------------
// Class category: Beamlines
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Beamlines/Beamline.h"


// Class Beamline
// ------------------------------------------------------------------------


Beamline::Beamline():
  ElementBase()
{
  shareFlag = false;
}


Beamline::Beamline(const Beamline &):
  ElementBase()
{
  shareFlag = false;
}


Beamline::Beamline(const string &name):
  ElementBase(name)
{
  shareFlag = false;
}


Beamline::~Beamline()
{}

