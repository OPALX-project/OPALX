// ------------------------------------------------------------------------
// $RCSfile: OpalBeamBeam.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalBeamBeam
//   The class for OPAL BEAMBEAM elements.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/08/06 20:32:28 $
// $Author: dbruhwil $
//
// ------------------------------------------------------------------------

#include "Elements/OpalBeamBeam.h"
#include "BeamlineCore/BeamBeamRep.h"
#include "Attributes/Attributes.h"


// Class OpalBeamBeam
// ------------------------------------------------------------------------


OpalBeamBeam::OpalBeamBeam():
    OpalElement(SIZE, "BEAMBEAM",
                "The \"BEAMBEAM\" element defines a beam-beam interaction.") {
    itsAttr[HORD] = Attributes::makeReal
                  ("HORD", "Horizontal displacement of opposite beam");
    itsAttr[VERTD] = Attributes::makeReal
                  ("VERTD", "Vertical displacement of opposite beam");
    itsAttr[SIGX] = Attributes::makeReal
                    ("SIGX", "Horizontal half-width of opposite beam");
    itsAttr[SIGY] = Attributes::makeReal
                    ("SIGY", "Vertical half-width of opposite beam");
    itsAttr[CHARGE] = Attributes::makeReal
                      ("CHARGE",
                       "Charge per particle in opposite beam, expressed in proton charges",
                       1.0);
    itsAttr[NPART] = Attributes::makeReal
                     ("NPART", "Number of particles in opposite beam");

    registerRealAttribute("HORD");
    registerRealAttribute("VERTD");
    registerRealAttribute("SIGX");
    registerRealAttribute("SIGY");
    registerRealAttribute("CHARGE");
    registerRealAttribute("NPART");

    registerOwnership();

    setElement(new BeamBeamRep("BEAMBEAM"));
}


OpalBeamBeam::OpalBeamBeam(const std::string &name, OpalBeamBeam *parent):
    OpalElement(name, parent) {
    setElement(new BeamBeamRep(name));
}


OpalBeamBeam::~OpalBeamBeam()
{}


OpalBeamBeam *OpalBeamBeam::clone(const std::string &name) {
    return new OpalBeamBeam(name, this);
}


void OpalBeamBeam::
fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);
    // *** MISSING: BeamBeam
}


void OpalBeamBeam::update() {
    OpalElement::update();

    BeamBeamRep *bb = dynamic_cast<BeamBeamRep *>(getElement());

    bb->setElementLength(0.0);
    Vector3D displacement(Attributes::getReal(itsAttr[HORD]),
                          Attributes::getReal(itsAttr[VERTD]),
                          0.0);
    bb->setBunchDisplacement(displacement);
    double sigx = Attributes::getReal(itsAttr[SIGX]);
    double sigy = Attributes::getReal(itsAttr[SIGY]);
    Matrix3D moments(sigx * sigx, 0.0, 0.0, 0.0, sigy * sigy, 0.0, 0.0, 0.0, 0.0);
    bb->setBunchMoment(moments);
    bb->setBunchCharge(Attributes::getReal(itsAttr[CHARGE]) *
                       Attributes::getReal(itsAttr[NPART]));

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(bb);
}