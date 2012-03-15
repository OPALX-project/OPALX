// ------------------------------------------------------------------------
// $RCSfile: CyclotronWrapper.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: CyclotronWrapper
//   Defines a representation for a modified Cyclotron.
//
// ------------------------------------------------------------------------
// Class category: ComponentWrappers
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 18:57:54 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "ComponentWrappers/CyclotronWrapper.h"
#include "AbsBeamline/BeamlineVisitor.h"


// Class CyclotronWrapper
// ------------------------------------------------------------------------

CyclotronWrapper::CyclotronWrapper(Cyclotron *ideal):
  Cyclotron(*ideal),
  itsDesign(ideal),
  itsError(),
  tempField()
{
  shareFlag = false;
}


CyclotronWrapper::CyclotronWrapper(const CyclotronWrapper &rhs):
  Cyclotron(rhs),
  itsDesign(rhs.itsDesign),
  itsError(rhs.itsError),
  tempField()
{
  shareFlag = false;
}


CyclotronWrapper::~CyclotronWrapper()
{}


void CyclotronWrapper::accept(BeamlineVisitor &visitor) const
{
  visitor.visitCyclotronWrapper(*this);
}


ElementBase *CyclotronWrapper::clone() const
{
  CyclotronWrapper *sbw = new CyclotronWrapper(*this);
  sbw->itsDesign = dynamic_cast<Cyclotron *>(itsDesign->clone());
  return sbw;
}


ElementBase *CyclotronWrapper::copyStructure()
{
  if (isSharable()) {
    return this;
  } else {
    CyclotronWrapper *sbw = new CyclotronWrapper(*this);
    sbw->itsDesign = dynamic_cast<Cyclotron *>(itsDesign->copyStructure());
    return sbw;
  }
}


BMultipoleField &CyclotronWrapper::errorField() const
{
  return itsError;
}


void CyclotronWrapper::makeSharable()
{
  shareFlag = true;
  itsDesign->makeSharable();
}


BMultipoleField &CyclotronWrapper::getField()
{
  /*
  tempField = itsDesign->getField();
  tempField.addField(itsError);
  return tempField;
  */
  BMultipoleField tf;
  return tf;
}


const BMultipoleField &CyclotronWrapper::getField() const
{
  /*
  tempField = itsDesign->getField();
  tempField.addField(itsError);
  return tempField;
  */
  BMultipoleField tf;
  return tf;
}


PlanarArcGeometry &CyclotronWrapper::getGeometry()
{
  /*
  return itsDesign->getGeometry();
  */
  PlanarArcGeometry gm(.0);
  return gm;
}


const PlanarArcGeometry &CyclotronWrapper::getGeometry() const
{
  /*
  return itsDesign->getGeometry();
  */
  PlanarArcGeometry gm(.0);
  return gm;
}


double CyclotronWrapper::getSlices()  const
{
  return itsDesign->getSlices();
}


double CyclotronWrapper::getStepsize()  const
{
  return itsDesign->getStepsize();
}


const string &CyclotronWrapper::getType() const
{
  static const string type("CyclotronWrapper");
  return type;
}


const Cyclotron &CyclotronWrapper::getDesign() const
{
  return *itsDesign;
}


Cyclotron &CyclotronWrapper::getDesign()
{
  return *itsDesign;
}


ElementBase *CyclotronWrapper::makeFieldWrapper()
{
  return this;
}


ElementBase *CyclotronWrapper::removeFieldWrapper()
{
  return &*itsDesign;
}


const ElementBase *CyclotronWrapper::removeFieldWrapper() const
{
  return &*itsDesign;
}


ElementBase *CyclotronWrapper::removeWrappers()
{
  return &*itsDesign;
}


const ElementBase *CyclotronWrapper::removeWrappers() const
{
  return &*itsDesign;
}
