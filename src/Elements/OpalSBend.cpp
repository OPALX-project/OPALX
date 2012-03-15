// ------------------------------------------------------------------------
// $RCSfile: OpalSBend.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalSBend
//   The class of OPAL rectangular bend magnets.
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:11 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalSBend.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/SBendRep.h"
#include "Fields/BMultipoleField.h"
#include "ComponentWrappers/SBendWrapper.h"
#include "Physics/Physics.h"
#include <cmath>


// Class OpalSBend
// ------------------------------------------------------------------------

OpalSBend::OpalSBend():
  OpalBend("SBEND",
	  "The \"SBEND\" element defines a sector bending magnet.")
{
  setElement((new SBendRep("SBEND"))->makeWrappers());
}


OpalSBend::OpalSBend(const string &name, OpalSBend *parent):
  OpalBend(name, parent)
{
  setElement((new SBendRep(name))->makeWrappers());
}


OpalSBend::~OpalSBend()
{}


OpalSBend *OpalSBend::clone(const string &name)
{
  return new OpalSBend(name, this);
}


void OpalSBend::
fillRegisteredAttributes(const ElementBase &base, ValueFlag flag)
{
  OpalElement::fillRegisteredAttributes(base, flag);

  // Get the desired field.
  const SBendWrapper *bend =
    dynamic_cast<const SBendWrapper *>(base.removeAlignWrapper());
  BMultipoleField field;

  // Get the desired field.
  if (flag == ERROR_FLAG) {
    field = bend->errorField();
  } else if (flag == ACTUAL_FLAG) {
    field = bend->getField();
  } else if (flag == IDEAL_FLAG) {
    field = bend->getDesign().getField();
  }

  double length = getLength();
  double scale = Physics::c / OPAL.getP0();
  if (length != 0.0) scale *= length;

  for (int i = 1; i <= field.order(); ++i) {
    string normName("K0L");
    normName[1] += (i - 1);
    attributeRegistry[normName]->setReal(scale * field.normal(i));

    string skewName("K0SL");
    skewName[1] += (i - 1);
    attributeRegistry[skewName]->setReal(scale * field.skew(i));
    scale *= double(i);
  }

  // Store pole face information.
  attributeRegistry["E1"]->setReal(bend->getEntryFaceRotation());
  attributeRegistry["E2"]->setReal(bend->getExitFaceRotation());
  attributeRegistry["H1"]->setReal(bend->getEntryFaceCurvature());
  attributeRegistry["H2"]->setReal(bend->getExitFaceCurvature());

  // Store integration parameters.
  attributeRegistry["SLICES"]->setReal(bend->getSlices());
  attributeRegistry["STEPSIZE"]->setReal(bend->getStepsize());
}


void OpalSBend::update()
{
  // Define geometry.
  SBendRep *bend =
    dynamic_cast<SBendRep *>(getElement()->removeWrappers());
  double length = Attributes::getReal(itsAttr[LENGTH]);
  double angle  = Attributes::getReal(itsAttr[ANGLE]);
  PlanarArcGeometry &geometry = bend->getGeometry();

  if (length) {
    geometry = PlanarArcGeometry(length, angle / length);
  } else {
    geometry = PlanarArcGeometry(angle);
  }

  // Define pole face angles.
  bend->setEntryFaceRotation(Attributes::getReal(itsAttr[E1]));
  bend->setExitFaceRotation(Attributes::getReal(itsAttr[E2]));
  bend->setEntryFaceCurvature(Attributes::getReal(itsAttr[H1]));
  bend->setExitFaceCurvature(Attributes::getReal(itsAttr[H2]));

  // Define integration parameters.
  bend->setSlices(Attributes::getReal(itsAttr[SLICES]));
  bend->setStepsize(Attributes::getReal(itsAttr[STEPSIZE]));

  // Define field.
  double factor = OPAL.getP0() / Physics::c;
  BMultipoleField field;
  double k0 = itsAttr[K0] ? Attributes::getReal(itsAttr[K0]) :
    length ? angle / length : angle;
  field.setNormalComponent(1, factor*k0);
  field.setSkewComponent  (1, factor*Attributes::getReal(itsAttr[K0S]));
  field.setNormalComponent(2, factor*Attributes::getReal(itsAttr[K1]));
  field.setSkewComponent  (2, factor*Attributes::getReal(itsAttr[K1S]));
  field.setNormalComponent(3, factor*Attributes::getReal(itsAttr[K2]) /2.0);
  field.setSkewComponent  (3, factor*Attributes::getReal(itsAttr[K2S])/2.0);
  field.setNormalComponent(4, factor*Attributes::getReal(itsAttr[K3]) /6.0);
  field.setSkewComponent  (4, factor*Attributes::getReal(itsAttr[K3S])/6.0);
  bend->setField(field);

  // Transmit "unknown" attributes.
  OpalElement::updateUnknown(bend);
}
