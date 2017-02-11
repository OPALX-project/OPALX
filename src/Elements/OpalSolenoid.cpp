// ------------------------------------------------------------------------
// $RCSfile: OpalSolenoid.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalSolenoid
//   The class of OPAL solenoids.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:40 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalSolenoid.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/SolenoidRep.h"
#include "Physics/Physics.h"


// Class OpalSolenoid
// ------------------------------------------------------------------------

OpalSolenoid::OpalSolenoid():
    OpalElement(SIZE, "SOLENOID",
                "The \"SOLENOID\" element defines a Solenoid.") {
    itsAttr[KS] = Attributes::makeReal
                  ("KS", "Normalised solenoid strength in m**(-1)");
    itsAttr[DKS] = Attributes::makeReal
                  ("DKS", "Normalised solenoid strength error in m**(-1)");

    registerRealAttribute("KS");
    registerRealAttribute("DKS");

    itsAttr[FMAPFN] = Attributes::makeString
                      ("FMAPFN", "Solenoid field map filename ");
    itsAttr[FAST] = Attributes::makeBool
                    ("FAST", "Faster but less accurate", true);

    registerStringAttribute("FMAPFN");

    setElement((new SolenoidRep("SOLENOID"))->makeAlignWrapper());
}


OpalSolenoid::OpalSolenoid(const std::string &name, OpalSolenoid *parent):
    OpalElement(name, parent) {
    setElement((new SolenoidRep(name))->makeAlignWrapper());
}


OpalSolenoid::~OpalSolenoid()
{}


OpalSolenoid *OpalSolenoid::clone(const std::string &name) {
    return new OpalSolenoid(name, this);
}


void OpalSolenoid::
fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);

    if(flag != ERROR_FLAG) {
        const SolenoidRep *sol =
            dynamic_cast<const SolenoidRep *>(base.removeWrappers());
        double length = sol->getElementLength();
        double ks = length * sol->getBz() * Physics::c / OpalData::getInstance()->getP0();
        attributeRegistry["KS"]->setReal(ks);
    }
}


void OpalSolenoid::update() {
    OpalElement::update();

    SolenoidRep *sol =
        dynamic_cast<SolenoidRep *>(getElement()->removeWrappers());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    double Bz = Attributes::getReal(itsAttr[KS]) * OpalData::getInstance()->getP0() / Physics::c;
    bool fast = Attributes::getBool(itsAttr[FAST]);

    sol->setElementLength(length);
    sol->setFieldMapFN(Attributes::getString(itsAttr[FMAPFN]));
    sol->setFast(fast);
    sol->setBz(Bz);
    sol->setKS(Attributes::getReal(itsAttr[KS]));
    sol->setDKS(Attributes::getReal(itsAttr[DKS]));

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(sol);
}
