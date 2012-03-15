// ------------------------------------------------------------------------
// $RCSfile: OpalCyclotronValley.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalCyclotronValley
//   The class of OPAL CyclotronValley for Multipacting Simulation.
//
// ------------------------------------------------------------------------
//
// $Date: 2010/12/8 14:47:39 $
// $Author: Chuan Wang $
//
// ------------------------------------------------------------------------

#include "Elements/OpalCyclotronValley.h"
#include "AbstractObjects/Attribute.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/CyclotronValleyRep.h"
#include "Physics/Physics.h"

extern Inform *gmsg;

// Class OpalCyclotronValley
// ------------------------------------------------------------------------

OpalCyclotronValley::OpalCyclotronValley():
    OpalElement(SIZE, "CYCLOTRONVALLEY",
                "The \"CYCLOTRONVALLEY\" element defines a CyclotronValley for Multipacting Simulation.") {
    
    itsAttr[FMAPFN] = Attributes::makeString
                      ("FMAPFN", "Filename for the fieldmap");
    itsAttr[DX] = Attributes::makeReal
                  ("DX", "Misalignment in x direction");
    itsAttr[DY] = Attributes::makeReal
                  ("DY", "Misalignment in y direction");
    itsAttr[DZ] = Attributes::makeReal
                  ("DZ", "Misalignment in z direction");
    itsAttr[BFLG] = Attributes::makeReal
                  ("BFLG", "B flag");
    registerStringAttribute("FMAPFN");
    registerRealAttribute("DX");
    registerRealAttribute("DY");
    registerRealAttribute("DZ");
    registerRealAttribute("BFLG");
    setElement((new CyclotronValleyRep("CyclotronValley"))->makeAlignWrapper());
}


OpalCyclotronValley::OpalCyclotronValley(const string &name, OpalCyclotronValley *parent):
    OpalElement(name, parent) {
    setElement((new CyclotronValleyRep(name))->makeAlignWrapper());
}


OpalCyclotronValley::~OpalCyclotronValley() {
   
}


OpalCyclotronValley *OpalCyclotronValley::clone(const string &name) {
    return new OpalCyclotronValley(name, this);
}


void OpalCyclotronValley::fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);

    if(flag != ERROR_FLAG) {
        const CyclotronValleyRep *cv =
            dynamic_cast<const CyclotronValleyRep *>(base.removeWrappers());
        attributeRegistry["FMAPFN"]->setString(cv->getFieldMapFN());
        double dx, dy, dz, bflg=1.0;
        cv->getMisalignment(dx, dy, dz);
        attributeRegistry["DX"]->setReal(dx);
        attributeRegistry["DY"]->setReal(dy);
        attributeRegistry["DZ"]->setReal(dz);

        attributeRegistry["BFLG"]->setReal(bflg);
    }
}


void OpalCyclotronValley::update() {
    using Physics::two_pi;
    CyclotronValleyRep *cv =
        dynamic_cast<CyclotronValleyRep *>(getElement()->removeWrappers());

   
    string fmapfm = Attributes::getString(itsAttr[FMAPFN]);
    double dx = Attributes::getReal(itsAttr[DX]);
    double dy = Attributes::getReal(itsAttr[DY]);
    double dz = Attributes::getReal(itsAttr[DZ]);
  
    cv->setMisalignment(dx, dy, dz);
  
    cv->setFieldMapFN(fmapfm);

    cv->setFast(false);//fast flag for cyclotronvalley has not been implemented yet.
    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(cv);
}
