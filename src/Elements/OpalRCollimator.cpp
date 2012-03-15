// ------------------------------------------------------------------------
// $RCSfile: OpalRCollimator.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalRCollimator
//   The class of OPAL rectangular collimators.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:40 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalRCollimator.h"
#include "AbstractObjects/Attribute.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/CollimatorRep.h"


// Class OpalRCollimator
// ------------------------------------------------------------------------

OpalRCollimator::OpalRCollimator():
  OpalElement(SIZE, "RCOLLIMATOR",
	     "The \"RCOLLIMATOR\" element defines a rectangular collimator.")
{
  itsAttr[XSIZE] = Attributes::makeReal
    ("XSIZE", "Horizontal half-aperture in m");
  itsAttr[YSIZE] = Attributes::makeReal
    ("YSIZE", "Vertical half-aperture in m");

  registerRealAttribute("XSIZE");
  registerRealAttribute("YSIZE");

  setElement((new CollimatorRep("RCOLLIMATOR"))->makeAlignWrapper());
}


OpalRCollimator::OpalRCollimator(const string &name, OpalRCollimator *parent):
  OpalElement(name, parent)
{
  setElement((new CollimatorRep(name))->makeAlignWrapper());
}


OpalRCollimator::~OpalRCollimator()
{}


OpalRCollimator *OpalRCollimator::clone(const string &name)
{
  return new OpalRCollimator(name, this);
}


void OpalRCollimator::
fillRegisteredAttributes(const ElementBase &base, ValueFlag flag)
{
  OpalElement::fillRegisteredAttributes(base, flag);

  const CollimatorRep *coll =
    dynamic_cast<const CollimatorRep *>(base.removeWrappers());
  attributeRegistry["XSIZE"]->setReal(coll->getXsize());
  attributeRegistry["YSIZE"]->setReal(coll->getYsize());
}


void OpalRCollimator::update()
{
  CollimatorRep *coll =
    dynamic_cast<CollimatorRep *>(getElement()->removeWrappers());
  coll->setElementLength(Attributes::getReal(itsAttr[LENGTH]));
  coll->setXsize(Attributes::getReal(itsAttr[XSIZE]));
  coll->setYsize(Attributes::getReal(itsAttr[YSIZE]));

  // Transmit "unknown" attributes.
  OpalElement::updateUnknown(coll);
}
