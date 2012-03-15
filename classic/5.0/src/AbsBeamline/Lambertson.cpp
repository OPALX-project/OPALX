// ------------------------------------------------------------------------
// $RCSfile: Lambertson.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Lambertson
//   Defines the abstract interface for a Lambertson septum magnet.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Lambertson.h"
#include "AbsBeamline/BeamlineVisitor.h"


// Class Lambertson
// ------------------------------------------------------------------------

Lambertson::Lambertson():
  Component()
{}


Lambertson::Lambertson(const Lambertson &rhs):
  Component(rhs)
{}


Lambertson::Lambertson(const string &name):
  Component(name)
{}


Lambertson::~Lambertson()
{}


void Lambertson::accept(BeamlineVisitor &visitor) const
{
  visitor.visitLambertson(*this);
}
