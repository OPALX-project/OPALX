// ------------------------------------------------------------------------
// $RCSfile: OpalRBend.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalRBend
//   The class of OPAL rectangular bend magnets.
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:11 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalRBend.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/RBendRep.h"
#include "Fields/BMultipoleField.h"
#include "ComponentWrappers/RBendWrapper.h"
#include "Physics/Physics.h"
#include "Structure/Wake.h"
#include <cmath>


// Class OpalRBend
// ------------------------------------------------------------------------

OpalRBend::OpalRBend():
  OpalBend("RBEND",
	  "The \"RBEND\" element defines a rectangular bending magnet.")
{
  setElement((new RBendRep("RBEND"))->makeWrappers());
}


OpalRBend::OpalRBend(const string &name, OpalRBend *parent):
  OpalBend(name, parent)
{
  setElement((new RBendRep(name))->makeWrappers());
}


OpalRBend::~OpalRBend()
{}


OpalRBend *OpalRBend::clone(const string &name)
{
  return new OpalRBend(name, this);
}


void OpalRBend::
fillRegisteredAttributes(const ElementBase &base, ValueFlag flag)
{
  OpalElement::fillRegisteredAttributes(base, flag);

  // Get the desired field.
  const RBendWrapper *bend =
    dynamic_cast<const RBendWrapper *>(base.removeAlignWrapper());
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


void OpalRBend::update()
{
  // Define geometry.
  RBendRep *bend =
    dynamic_cast<RBendRep *>(getElement()->removeWrappers());
  double length = Attributes::getReal(itsAttr[LENGTH]);
  double angle  = Attributes::getReal(itsAttr[ANGLE]);
  RBendGeometry &geometry = bend->getGeometry();
  geometry.setElementLength(length);
  geometry.setBendAngle(angle);

  bend->setLength(length);

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
  double k0 =
    itsAttr[K0] ? Attributes::getReal(itsAttr[K0]) :
    length ? 2*sin(angle/2) / length : angle;
                //JMJ 4/10/2000: above line replaced
                //    length ? angle / length : angle;
                // to avoid closed orbit created by RBEND with defalt K0.
  field.setNormalComponent(1, factor*k0);
  field.setSkewComponent  (1, factor*Attributes::getReal(itsAttr[K0S]));
  field.setNormalComponent(2, factor*Attributes::getReal(itsAttr[K1]));
  field.setSkewComponent  (2, factor*Attributes::getReal(itsAttr[K1S]));
  field.setNormalComponent(3, factor*Attributes::getReal(itsAttr[K2]) /2.0);
  field.setSkewComponent  (3, factor*Attributes::getReal(itsAttr[K2S])/2.0);
  field.setNormalComponent(4, factor*Attributes::getReal(itsAttr[K3]) /6.0);
  field.setSkewComponent  (4, factor*Attributes::getReal(itsAttr[K3S])/6.0);
  bend->setField(field);

//   const std::vector<double> enge = Attributes::getRealArray(itsAttr[ENGECOEF]);
  bend->setAmplitudem(k0);
  bend->setGapWidth(Attributes::getReal(itsAttr[HGAP]));
  
  if (itsAttr[FMAPFN])
    bend->setFieldMapFN(Attributes::getString(itsAttr[FMAPFN]));
  else if (bend->getName() != "RBEND")
    cerr << bend->getName() << ": No filename for a field map given" << endl;

  if (itsAttr[ALPHA])
    bend->setAlpha(Attributes::getReal(itsAttr[ALPHA]));
  else
    bend->setAlpha(0.0);

  if (itsAttr[BETA])
    bend->setBeta(Attributes::getReal(itsAttr[BETA]));
  else
    bend->setBeta(0.0);

  bend->setDesignEnergy(Attributes::getReal(itsAttr[DESIGNENERGY])*1e6); // convert from MeV to eV

  if (itsAttr[WAKEF])
    bend->setWake(Wake::find(Attributes::getString(itsAttr[WAKEF])));
                    
  // Transmit "unknown" attributes.
  OpalElement::updateUnknown(bend);
}
