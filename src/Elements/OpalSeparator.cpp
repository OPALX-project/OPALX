// ------------------------------------------------------------------------
// $RCSfile: OpalSeparator.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalSeparator
//   The class of OPAL electrostatic separators.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:40 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalSeparator.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/SeparatorRep.h"


// Class OpalSeparator
// ------------------------------------------------------------------------

OpalSeparator::OpalSeparator():
    OpalElement(SIZE, "SEPARATOR",
                "The \"SEPARATOR\" element defines an electrostatic separator.") {
    itsAttr[EX] = Attributes::makeReal
                  ("EX", "The horizontal electrostatic field in MV");
    itsAttr[EY] = Attributes::makeReal
                  ("EY", "The vertical electrostatic field in MV");

    registerRealAttribute("EXL");
    registerRealAttribute("EYL");

    registerOwnership();

    setElement((new SeparatorRep("SEPARATOR"))->makeAlignWrapper());
}


OpalSeparator::OpalSeparator(const std::string &name, OpalSeparator *parent):
    OpalElement(name, parent) {
    setElement((new SeparatorRep(name))->makeAlignWrapper());
}


OpalSeparator::~OpalSeparator()
{}


OpalSeparator *OpalSeparator::clone(const std::string &name) {
    return new OpalSeparator(name, this);
}


void OpalSeparator::
fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);

    if(flag != ERROR_FLAG) {
        const SeparatorRep *sep =
            dynamic_cast<const SeparatorRep *>(base.removeWrappers());
        double length = sep->getElementLength();
        attributeRegistry["EXL"]->setReal(length * sep->getEx());
        attributeRegistry["EYL"]->setReal(length * sep->getEy());
    }
}


void OpalSeparator::update() {
    OpalElement::update();

    SeparatorRep *sep =
        dynamic_cast<SeparatorRep *>(getElement()->removeWrappers());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    double Ex     = Attributes::getReal(itsAttr[EX]) * 1.0e6;
    double Ey     = Attributes::getReal(itsAttr[EY]) * 1.0e6;
    sep->setElementLength(length);
    sep->setEx(Ex);
    sep->setEy(Ey);

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(sep);
}