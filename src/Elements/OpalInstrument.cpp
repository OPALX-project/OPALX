// ------------------------------------------------------------------------
// $RCSfile: OpalInstrument.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalInstrument
//   The class of OPAL generic beam instrumentation.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalInstrument.h"
#include "AbstractObjects/Attribute.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/DriftRep.h"


// Class OpalInstrument
// ------------------------------------------------------------------------

OpalInstrument::OpalInstrument():
    OpalElement(COMMON, "INSTRUMENT",
                "The \"INSTRUMENT\" element defines a generic "
                "beam observation device.") {
    setElement(new DriftRep("INSTRUMENT"));
}


OpalInstrument::OpalInstrument(const std::string &name, OpalInstrument *parent):
    OpalElement(name, parent) {
    setElement(new DriftRep(name));
}


OpalInstrument::~OpalInstrument()
{}


OpalInstrument *OpalInstrument::clone(const std::string &name) {
    return new OpalInstrument(name, this);
}


void OpalInstrument::update() {
    DriftRep *ins = dynamic_cast<DriftRep *>(getElement());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    ins->setElementLength(length);

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(ins);
}
