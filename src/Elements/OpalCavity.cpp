// ------------------------------------------------------------------------
// $RCSfile: OpalCavity.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalCavity
//   The class of OPAL RF cavities.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalCavity.h"
#include "AbstractObjects/Attribute.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/RFCavityRep.h"
#include "Structure/Wake.h"
#include "Physics/Physics.h"


// Class OpalCavity
// ------------------------------------------------------------------------

OpalCavity::OpalCavity():
  OpalElement(SIZE, "RFCAVITY",
              "The \"RFCAVITY\" element defines an RF cavity.")
{
  itsAttr[VOLT] = Attributes::makeReal
    ("VOLT", "RF voltage in MV");
  itsAttr[FREQ] = Attributes::makeReal
    ("FREQ", "RF frequency in MHz");
  itsAttr[LAG] = Attributes::makeReal
    ("LAG", "Phase lag, in multiples of (2*pi)");
  itsAttr[HARMON] = Attributes::makeReal
    ("HARMON", "Harmonic number");
  itsAttr[BETARF] = Attributes::makeReal
    ("BETRF", "beta_RF");
  itsAttr[PG] = Attributes::makeReal
    ("PG", "RF power in MW");
  itsAttr[ZSHUNT] = Attributes::makeReal
    ("SHUNT", "Shunt impedance in MOhm");
  itsAttr[TFILL] = Attributes::makeReal
    ("TFILL", "Fill time in microseconds");
  itsAttr[FMAPFN] = Attributes::makeString
    ("FMAPFN", "Filename for the fieldmap");
  itsAttr[FAST] = Attributes::makeBool
    ("FAST", "Faster but less accurate");
  itsAttr[CAVITYTYPE] = Attributes::makeString
    ("CAVITYTYPE", "STANDING or SINGLEGAP cavity in photoinjector and LINAC; SINGLEGAP or DOUBLEGAP cavity in cyclotron");

  itsAttr[RMIN] = Attributes::makeReal
    ("RMIN", " Minimal Radius of a cyclotron cavity");
  itsAttr[RMAX] = Attributes::makeReal
    ("RMAX", " Maximal Radius of a cyclotron cavity");
  itsAttr[ANGLE] = Attributes::makeReal
    ("ANGLE", "Azimuth position of a cyclotron cavity");
  itsAttr[PDIS] = Attributes::makeReal
    ("PDIS", "Shift distance of cavity gap from center of cyclotron");
  itsAttr[GAPWIDTH] = Attributes::makeReal
    ("GAPWIDTH", "Gap width of a cyclotron cavity");
  itsAttr[PHI0] = Attributes::makeReal
    ("PHI0","initial phase of cavity");

  itsAttr[DX] = Attributes::makeReal
    ("DX", "Misalignment in x direction");
  itsAttr[DY] = Attributes::makeReal
    ("DY", "Misalignment in y direction");
  itsAttr[DZ] = Attributes::makeReal
    ("DZ", "Misalignment in z direction");

  registerRealAttribute("VOLT");
  registerRealAttribute("FREQ");
  registerRealAttribute("LAG");
  registerStringAttribute("FMAPFN");
  registerStringAttribute("CAVITYTYPE");
  registerRealAttribute("RMIN");
  registerRealAttribute("RMAX");
  registerRealAttribute("ANGLE");
  registerRealAttribute("PDIS");
  registerRealAttribute("GAPWIDTH");
  registerRealAttribute("PHI0");
  registerRealAttribute("DX");
  registerRealAttribute("DY");
  registerRealAttribute("DZ");

  setElement((new RFCavityRep("RFCAVITY"))->makeAlignWrapper());
}


OpalCavity::OpalCavity(const string &name, OpalCavity *parent):
  OpalElement(name, parent)
{
  setElement((new RFCavityRep(name))->makeAlignWrapper());
}


OpalCavity::~OpalCavity()
{}


OpalCavity *OpalCavity::clone(const string &name)
{
  return new OpalCavity(name, this);
}


void OpalCavity::fillRegisteredAttributes(const ElementBase &base, ValueFlag flag)
{
  OpalElement::fillRegisteredAttributes(base, flag);

  if (flag != ERROR_FLAG) {
    const RFCavityRep *rfc =
      dynamic_cast<const RFCavityRep *>(base.removeWrappers());
    attributeRegistry["VOLT"]->setReal(rfc->getAmplitude());
    attributeRegistry["FREQ"]->setReal(rfc->getFrequency());
    attributeRegistry["LAG"]->setReal(rfc->getPhase());
    attributeRegistry["FMAPFN"]->setString(rfc->getFieldMapFN());
    attributeRegistry["CAVITYTYPE"]->setString(rfc->getCavityType());
    double dx, dy, dz;
    rfc->getMisalignment(dx,dy,dz);
    cout << "Cav dx: " << dx << " Cav dy: " << dy << " Cav dz: " << dz << endl;
    attributeRegistry["DX"]->setReal(dx);
    attributeRegistry["DY"]->setReal(dy);
    attributeRegistry["DZ"]->setReal(dz);
    
  }
}


void OpalCavity::update()
{
  using Physics::two_pi;
  RFCavityRep *rfc =
    dynamic_cast<RFCavityRep *>(getElement()->removeWrappers());

  double length = Attributes::getReal(itsAttr[LENGTH]);
  double vPeak  = 1.0e6 * Attributes::getReal(itsAttr[VOLT]);
  double phase  = two_pi * Attributes::getReal(itsAttr[LAG]);
  double freq   = (1.0e6 * two_pi) * Attributes::getReal(itsAttr[FREQ]);
  string fmapfm = Attributes::getString(itsAttr[FMAPFN]);
  string type = Attributes::getString(itsAttr[TYPE]);
  bool fast = Attributes::getBool(itsAttr[FAST]);

  double rmin = Attributes::getReal(itsAttr[RMIN]);
  double rmax = Attributes::getReal(itsAttr[RMAX]);
  double angle = Attributes::getReal(itsAttr[ANGLE]);
  double pdis = Attributes::getReal(itsAttr[PDIS]);
  double gapwidth = Attributes::getReal(itsAttr[GAPWIDTH]);
  double phi0 = Attributes::getReal(itsAttr[PHI0]);
  double dx = Attributes::getReal(itsAttr[DX]);
  double dy = Attributes::getReal(itsAttr[DY]);
  double dz = Attributes::getReal(itsAttr[DZ]);

  
  if (itsAttr[WAKEF])
    rfc->setWake(Wake::find(Attributes::getString(itsAttr[WAKEF])));

  rfc->setMisalignment(dx, dy, dz);

  rfc->setElementLength(length);
  rfc->setAmplitude(vPeak);
  rfc->setFrequency(freq);
  rfc->setPhase(phase);

  rfc->setFieldMapFN(fmapfm);
  rfc->setFast(fast);
  rfc->setAmplitudem(vPeak);
  rfc->setFrequencym(freq);
  rfc->setPhasem(phase);
  rfc->setCavityType(type);

  rfc->setComponentType(type);
  rfc->setRmin(rmin);
  rfc->setRmax(rmax);
  rfc->setAzimuth(angle);
  rfc->setPerpenDistance(pdis);
  rfc->setGapWidth(gapwidth);
  rfc->setPhi0(phi0);
  
  // Transmit "unknown" attributes.
  OpalElement::updateUnknown(rfc);
}
