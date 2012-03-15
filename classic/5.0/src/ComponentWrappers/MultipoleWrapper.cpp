// ------------------------------------------------------------------------
// $RCSfile: MultipoleWrapper.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MultipoleWrapper
//   Defines a representation for a modified multipole.
//
// ------------------------------------------------------------------------
// Class category: ComponentWrappers
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:35 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "ComponentWrappers/MultipoleWrapper.h"
#include "AbsBeamline/BeamlineVisitor.h"


// Class MultipoleWrapper
//   A MultipoleWrapper represents a unique instance of a multipole magnet
//   in the accelerator model.
// ------------------------------------------------------------------------

MultipoleWrapper::MultipoleWrapper(Multipole *ideal):
  Multipole(*ideal),
  itsDesign(ideal),
  itsError(),
  tempField()
{
  setName(ideal->getName());
  shareFlag = false;
}


MultipoleWrapper::MultipoleWrapper(const MultipoleWrapper &rhs):
  Multipole(rhs),
  itsDesign(rhs.itsDesign),
  itsError(rhs.itsError),
  tempField()
{
  shareFlag = false;
}


MultipoleWrapper::~MultipoleWrapper()
{}


void MultipoleWrapper::accept(BeamlineVisitor &visitor) const
{
  visitor.visitMultipoleWrapper(*this);
}


ElementBase *MultipoleWrapper::clone() const
{
  MultipoleWrapper *mw = new MultipoleWrapper(*this);
  mw->itsDesign = dynamic_cast<Multipole *>(itsDesign->clone());
  return mw;
}


ElementBase *MultipoleWrapper::copyStructure()
{
  if (isSharable()) {
    return this;
  } else {
    MultipoleWrapper *mw = new MultipoleWrapper(*this);
    mw->itsDesign = dynamic_cast<Multipole *>(itsDesign->copyStructure());
    return mw;
  }
}


BMultipoleField &MultipoleWrapper::errorField() const
{
  return itsError;

}


void MultipoleWrapper::makeSharable()
{
  shareFlag = true;
  itsDesign->makeSharable();
}


BMultipoleField &MultipoleWrapper::getField()
{
  tempField = itsDesign->getField();
  tempField.addField(itsError);
  return tempField;
}


const BMultipoleField &MultipoleWrapper::getField() const
{
  tempField = itsDesign->getField();
  tempField.addField(itsError);
  return tempField;
}


StraightGeometry &MultipoleWrapper::getGeometry()
{
  return itsDesign->getGeometry();
}


const StraightGeometry &MultipoleWrapper::getGeometry() const
{
  return itsDesign->getGeometry();
}


const string &MultipoleWrapper::getType() const
{
  static const string type("MultipoleWrapper");
  return type;
}


const Multipole &MultipoleWrapper::getDesign() const
{
  return *itsDesign;
}


Multipole &MultipoleWrapper::getDesign()
{
  return *itsDesign;
}


ElementBase *MultipoleWrapper::makeFieldWrapper()
{
  return this;
}


ElementBase *MultipoleWrapper::removeFieldWrapper()
{
  return &*itsDesign;
}


const ElementBase *MultipoleWrapper::removeFieldWrapper() const
{
  return &*itsDesign;
}


ElementBase *MultipoleWrapper::removeWrappers()
{
  return &*itsDesign;
}


const ElementBase *MultipoleWrapper::removeWrappers() const
{
  return &*itsDesign;
}
