// ------------------------------------------------------------------------
// $RCSfile: OpalDrift.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalDrift
//   The class of OPAL drift spaces.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalDrift.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/DriftRep.h"
#include "Structure/Wake.h"

// Class OpalDrift
// ------------------------------------------------------------------------

OpalDrift::OpalDrift():
    OpalElement(COMMON, "DRIFT",
            "The \"DRIFT\" element defines a drift space.")
{
    itsAttr[LENGTH] = Attributes::makeReal
        ("LENGTH", "Drift length");

    registerRealAttribute("LENGTH");

    setElement(new DriftRep("DRIFT"));
}


OpalDrift::OpalDrift(const string &name, OpalDrift *parent):
    OpalElement(name, parent)
{
    setElement(new DriftRep(name));
}


OpalDrift::~OpalDrift()
{}


OpalDrift *OpalDrift::clone(const string &name)
{
    return new OpalDrift(name, this);
}


bool OpalDrift::isDrift() const
{
    return true;
}


void OpalDrift::update()
{
    DriftRep *drf = static_cast<DriftRep *>(getElement());

    drf->setElementLength(Attributes::getReal(itsAttr[LENGTH]));

    if (itsAttr[WAKEF]) 
        drf->setWake(Wake::find(Attributes::getString(itsAttr[WAKEF])));

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(drf);
}
