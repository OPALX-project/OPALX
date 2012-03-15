// ------------------------------------------------------------------------
// $RCSfile: OpalCollimator.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalECollimator
//   The class of OPAL elliptic collimators.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalECollimator.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/CollimatorRep.h"


// Class OpalECollimator
// ------------------------------------------------------------------------

OpalECollimator::OpalECollimator():
  OpalElement(SIZE, "ECOLLIMATOR",
	     "The \"ECOLLIMATOR\" element defines an elliptic collimator.")
{
  itsAttr[XSIZE] = Attributes::makeReal
    ("XSIZE", "Horizontal half-aperture in m");
  itsAttr[YSIZE] = Attributes::makeReal
    ("YSIZE", "Vertical half-aperture in m");
  itsAttr[OUTFN] = Attributes::makeString
    ("OUTFN", "Monitor output filename");

  registerStringAttribute("OUTFN");
  registerRealAttribute("XSIZE");
  registerRealAttribute("YSIZE");

  setElement((new CollimatorRep("ECOLLIMATOR"))->makeAlignWrapper());
}


OpalECollimator::OpalECollimator(const string &name, OpalECollimator *parent):
  OpalElement(name, parent)
{
  setElement((new CollimatorRep(name))->makeAlignWrapper());
}


OpalECollimator::~OpalECollimator()
{}


OpalECollimator *OpalECollimator::clone(const string &name)
{
  return new OpalECollimator(name, this);
}


void OpalECollimator::fillRegisteredAttributes(const ElementBase &base, ValueFlag flag)
{
  OpalElement::fillRegisteredAttributes(base, flag);

  const CollimatorRep *coll =
    dynamic_cast<const CollimatorRep *>(base.removeWrappers());
  attributeRegistry["XSIZE"]->setReal(coll->getXsize());
  attributeRegistry["YSIZE"]->setReal(coll->getYsize());
}


void OpalECollimator::update()
{
  CollimatorRep *coll =
    dynamic_cast<CollimatorRep *>(getElement()->removeWrappers());
  double length = Attributes::getReal(itsAttr[LENGTH]);
  coll->setElementLength(length);
  coll->setXsize(Attributes::getReal(itsAttr[XSIZE]));
  coll->setYsize(Attributes::getReal(itsAttr[YSIZE]));
  coll->setOutputFN(Attributes::getString(itsAttr[OUTFN]));

  // Transmit "unknown" attributes.
  OpalElement::updateUnknown(coll);
}
