// ------------------------------------------------------------------------
// $RCSfile: OpalTravelingWave.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalTravelingWave
//   The class of OPAL RF cavities.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalTravelingWave.h"
#include "AbstractObjects/Attribute.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/TravelingWaveRep.h"
#include "Structure/OpalWake.h"
#include "Physics/Physics.h"


// Class OpalTravelingWave
// ------------------------------------------------------------------------

OpalTravelingWave::OpalTravelingWave():
    OpalElement(SIZE, "TRAVELINGWAVE",
                "The \"TRAVELINGWAVE\" element defines a traveling wave structure."),
    owk_m(NULL) {
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
                    ("FAST", "Faster but less accurate", true);
    itsAttr[APVETO] = Attributes::makeBool
                    ("APVETO", "Do not use this cavity in the Autophase procedure", false);
    itsAttr[CAVITYTYPE] = Attributes::makeString
                          ("CAVITYTYPE", "STANDING or TRAVELING wave cavity in photoinjector and LINAC; SINGLEGAP or DOUBLEGAP cavity in cyclotron");
    itsAttr[NUMCELLS] = Attributes::makeReal
                        ("NUMCELLS", "Number of cells in a TW structure");
    itsAttr[DX] = Attributes::makeReal
      ("DX", "Misalignment in x direction",0.0);
    itsAttr[DY] = Attributes::makeReal
      ("DY", "Misalignment in y direction",0.0);

    registerRealAttribute("VOLT");
    registerRealAttribute("FREQ");
    registerRealAttribute("LAG");
    registerStringAttribute("FMAPFN");
    registerStringAttribute("CAVITYTYPE");
    registerRealAttribute("NUMCELLS");
    registerRealAttribute("DX");
    registerRealAttribute("DY");

    setElement((new TravelingWaveRep("TRAVELINGWAVE"))->makeAlignWrapper());
}


OpalTravelingWave::OpalTravelingWave(const string &name, OpalTravelingWave *parent):
    OpalElement(name, parent),
    owk_m(NULL) {
    setElement((new TravelingWaveRep(name))->makeAlignWrapper());
}


OpalTravelingWave::~OpalTravelingWave() {
    if(owk_m)
        delete owk_m;
}


OpalTravelingWave *OpalTravelingWave::clone(const string &name) {
    return new OpalTravelingWave(name, this);
}


void OpalTravelingWave::
fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);

    if(flag != ERROR_FLAG) {
        const TravelingWaveRep *rfc =
            dynamic_cast<const TravelingWaveRep *>(base.removeWrappers());
        attributeRegistry["VOLT"]->setReal(rfc->getAmplitude());
        attributeRegistry["FREQ"]->setReal(rfc->getFrequency());
        attributeRegistry["LAG"]->setReal(rfc->getPhase());
        attributeRegistry["FMAPFN"]->setString(rfc->getFieldMapFN());
        double dx, dy, dz;
        rfc->getMisalignment(dx, dy, dz);
        attributeRegistry["DX"]->setReal(dx);
        attributeRegistry["DY"]->setReal(dy);
    }
}


void OpalTravelingWave::update() {
    using Physics::two_pi;
    TravelingWaveRep *rfc =
        dynamic_cast<TravelingWaveRep *>(getElement()->removeWrappers());

    double length = Attributes::getReal(itsAttr[LENGTH]);
    double vPeak  = Attributes::getReal(itsAttr[VOLT]);
    double phase  = Attributes::getReal(itsAttr[LAG]);
    double freq   = (1.0e6 * two_pi) * Attributes::getReal(itsAttr[FREQ]);
    string fmapfm = Attributes::getString(itsAttr[FMAPFN]);
    bool fast = Attributes::getBool(itsAttr[FAST]);
    bool apVeto = Attributes::getBool(itsAttr[APVETO]);

    string type = Attributes::getString(itsAttr[TYPE]);
    double dx = Attributes::getReal(itsAttr[DX]);
    double dy = Attributes::getReal(itsAttr[DY]);

    rfc->setMisalignment(dx, dy, 0.0);

    rfc->setElementLength(length);
    rfc->setAmplitude(1.0e6 * vPeak);
    rfc->setFrequency(freq);
    rfc->setPhase(phase);

    rfc->setFieldMapFN(fmapfm);
    rfc->setFast(fast);
    rfc->setAutophaseVeto(apVeto);
    rfc->setAmplitudem(vPeak);
    rfc->setFrequencym(freq);
    rfc->setPhasem(phase);
    rfc->setNumCells((int)Attributes::getReal(itsAttr[NUMCELLS]));

    if(itsAttr[WAKEF] && owk_m == NULL) {
        owk_m = (OpalWake::find(Attributes::getString(itsAttr[WAKEF])))->clone(getOpalName() + string("_wake"));
        owk_m->initWakefunction(*rfc);
        rfc->setWake(owk_m->wf_m);
    }
    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(rfc);
}
