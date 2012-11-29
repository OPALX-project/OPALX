// ------------------------------------------------------------------------
// $RCSfile: OpalParallelPlate.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalParallelPlate
//   The class of OPAL  cavities.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalParallelPlate.h"
#include "Structure/BoundaryGeometry.h"
#include "AbstractObjects/Attribute.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/ParallelPlateRep.h"
#include "Physics/Physics.h"

extern Inform *gmsg;

// Class OpalParallelPlate
// ------------------------------------------------------------------------

OpalParallelPlate::OpalParallelPlate():
    OpalElement(SIZE, "PARALLELPLATE",
                "The \"PARALLELPLATE\" element defines an  cavity.") {
    itsAttr[VOLT] = Attributes::makeReal
                    ("VOLT", " voltage in MV");
    itsAttr[FREQ] = Attributes::makeReal
                    ("FREQ", " frequency in MHz");
    itsAttr[LAG] = Attributes::makeReal
                   ("LAG", "Phase lag (rad), !!!! was before in multiples of (2*pi) !!!!");

    itsAttr[GEOMETRY] = Attributes::makeString
                        ("GEOMETRY", "BoundaryGeometry for ParallelPlate");

    itsAttr[PLENGTH] = Attributes::makeReal
                       ("PLENGTH", " Gap length in Meter");
    itsAttr[DX] = Attributes::makeReal
      ("DX", "Misalignment in x direction",0.0);
    itsAttr[DY] = Attributes::makeReal
      ("DY", "Misalignment in y direction",0.0);
    itsAttr[DZ] = Attributes::makeReal
      ("DZ", "Misalignment in z direction",0.0);


    registerRealAttribute("VOLT");
    registerRealAttribute("FREQ");
    registerRealAttribute("LAG");
    registerStringAttribute("GEOMETRY");
    registerRealAttribute("PLENGTH");

    registerRealAttribute("DX");
    registerRealAttribute("DY");
    registerRealAttribute("DZ");
    setElement((new ParallelPlateRep("ParallelPlate"))->makeAlignWrapper());
}


OpalParallelPlate::OpalParallelPlate(const string &name, OpalParallelPlate *parent):
    OpalElement(name, parent) {
    setElement((new ParallelPlateRep(name))->makeAlignWrapper());
}


OpalParallelPlate::~OpalParallelPlate() {

}


OpalParallelPlate *OpalParallelPlate::clone(const string &name) {
    return new OpalParallelPlate(name, this);
}


void OpalParallelPlate::fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);

    if(flag != ERROR_FLAG) {
        const ParallelPlateRep *pplate =
            dynamic_cast<const ParallelPlateRep *>(base.removeWrappers());
        attributeRegistry["VOLT"]->setReal(pplate->getAmplitude());
        attributeRegistry["FREQ"]->setReal(pplate->getFrequency());
        attributeRegistry["LAG"]->setReal(pplate->getPhase());
        attributeRegistry["PLENGTH"]->setReal(pplate->getElementLength());
        double dx, dy, dz;
        pplate->getMisalignment(dx, dy, dz);
        attributeRegistry["DX"]->setReal(dx);
        attributeRegistry["DY"]->setReal(dy);
        attributeRegistry["DZ"]->setReal(dz);
    }
}


void OpalParallelPlate::update() {
    using Physics::two_pi;
    ParallelPlateRep *pplate =
        dynamic_cast<ParallelPlateRep *>(getElement()->removeWrappers());

    double vPeak  = Attributes::getReal(itsAttr[VOLT]);
    //    double phase  = two_pi * Attributes::getReal(itsAttr[LAG]);
    double phase  = Attributes::getReal(itsAttr[LAG]);
    double freq   = (1.0e6 * two_pi) * Attributes::getReal(itsAttr[FREQ]);
    double length = Attributes::getReal(itsAttr[PLENGTH]);
    double dx = Attributes::getReal(itsAttr[DX]);
    double dy = Attributes::getReal(itsAttr[DY]);
    double dz = Attributes::getReal(itsAttr[DZ]);

    if(itsAttr[GEOMETRY] && obgeo_m == NULL) {
        obgeo_m = (BoundaryGeometry::find(Attributes::getString(itsAttr[GEOMETRY])))->clone(getOpalName() + string("_geometry"));
        if(obgeo_m) {
            //obgeo_m->initialize();

            pplate->setBoundaryGeometry(obgeo_m);
        }
    }

    pplate->setMisalignment(dx, dy, dz);
    pplate->setAmplitude(1.0e6 * vPeak);
    pplate->setFrequency(freq);
    pplate->setPhase(phase);

    pplate->setElementLength(length);

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(pplate);
}
