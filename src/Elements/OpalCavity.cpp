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
#include "Structure/OpalWake.h"
#include "Structure/BoundaryGeometry.h"
#include "Physics/Physics.h"

extern Inform *gmsg;

// Class OpalCavity
// ------------------------------------------------------------------------

OpalCavity::OpalCavity():
    OpalElement(SIZE, "RFCAVITY",
                "The \"RFCAVITY\" element defines an RF cavity."),
    owk_m(NULL),
    obgeo_m(NULL) {
    itsAttr[VOLT] = Attributes::makeRealArray
                    ("VOLT", "RF voltage in MV");
    itsAttr[FREQ] = Attributes::makeRealArray
	            ("FREQ", "RF frequency in MHz");
    itsAttr[LAG] = Attributes::makeRealArray
                   ("LAG", "Phase lag (rad), !!!! was before in multiples of (2*pi) !!!!");
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
    itsAttr[FMAPFN] = Attributes::makeStringArray
                      ("FMAPFN", "Filenames of the fieldmap(s)");
    itsAttr[GEOMETRY] = Attributes::makeString
                        ("GEOMETRY", "BoundaryGeometry for Cavities");
    itsAttr[FAST] = Attributes::makeBool
                    ("FAST", "Faster but less accurate", true);
    itsAttr[APVETO] = Attributes::makeBool
                    ("APVETO", "Do not use this cavity in the Autophase procedure", false);
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
                    ("PHI0", "initial phase of cavity");

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
    registerStringAttribute("GEOMETRY");
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
    OpalElement(name, parent),
    owk_m(NULL),
    obgeo_m(NULL) {
    setElement((new RFCavityRep(name))->makeAlignWrapper());
}


OpalCavity::~OpalCavity() {
    if(owk_m)
        delete owk_m;
}


OpalCavity *OpalCavity::clone(const string &name) {
    return new OpalCavity(name, this);
}


void OpalCavity::fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);

    if(flag != ERROR_FLAG) {
        const RFCavityRep *rfc =
            dynamic_cast<const RFCavityRep *>(base.removeWrappers());
        attributeRegistry["VOLT"]->setReal(rfc->getAmplitude());
        attributeRegistry["FREQ"]->setReal(rfc->getFrequency());
        attributeRegistry["LAG"]->setReal(rfc->getPhase());
        attributeRegistry["FMAPFN"]->setString(rfc->getFieldMapFN());
        attributeRegistry["CAVITYTYPE"]->setString(rfc->getCavityType());
        double dx, dy, dz;
        rfc->getMisalignment(dx, dy, dz);
        attributeRegistry["DX"]->setReal(dx);
        attributeRegistry["DY"]->setReal(dy);
        attributeRegistry["DZ"]->setReal(dz);
    }
}


void OpalCavity::update() {
    using Physics::two_pi;
    RFCavityRep *rfc =
        dynamic_cast<RFCavityRep *>(getElement()->removeWrappers());

    double length = Attributes::getReal(itsAttr[LENGTH]);
    std::vector<double> vPeak  = Attributes::getRealArray(itsAttr[VOLT]);
    std::vector<double> phase  = Attributes::getRealArray(itsAttr[LAG]);
    std::vector<double> freq   = Attributes::getRealArray(itsAttr[FREQ]);
    std::vector<string> fmapfns = Attributes::getStringArray(itsAttr[FMAPFN]);
    string type = Attributes::getString(itsAttr[TYPE]);
    bool fast = Attributes::getBool(itsAttr[FAST]);
    double max_freq = 0.0;
    for (size_t i = 0; i < freq.size(); ++ i) max_freq = std::max(std::abs(freq[i]), max_freq);
    bool apVeto = (Attributes::getBool(itsAttr[APVETO]) || max_freq < 1e-9);

    double rmin = Attributes::getReal(itsAttr[RMIN]);
    double rmax = Attributes::getReal(itsAttr[RMAX]);
    double angle = Attributes::getReal(itsAttr[ANGLE]);
    double pdis = Attributes::getReal(itsAttr[PDIS]);
    double gapwidth = Attributes::getReal(itsAttr[GAPWIDTH]);
    double phi0 = Attributes::getReal(itsAttr[PHI0]);
    double dx = Attributes::getReal(itsAttr[DX]);
    double dy = Attributes::getReal(itsAttr[DY]);
    double dz = Attributes::getReal(itsAttr[DZ]);

    if(itsAttr[WAKEF] && owk_m == NULL) {
        owk_m = (OpalWake::find(Attributes::getString(itsAttr[WAKEF])))->clone(getOpalName() + string("_wake"));
        owk_m->initWakefunction(*rfc);
        rfc->setWake(owk_m->wf_m);
    }
    
    if(itsAttr[GEOMETRY] && obgeo_m == NULL) {
        obgeo_m = (BoundaryGeometry::find(Attributes::getString(itsAttr[GEOMETRY])))->clone(getOpalName() + string("_geometry"));
        if(obgeo_m) {
	    rfc->setBoundaryGeometry(obgeo_m);
        }
    }

    rfc->setMisalignment(dx, dy, dz);

    rfc->setElementLength(length);

    double peak = 0.0, frequency = 0.0, phi = 0.0;
    if (vPeak.size() > 0) peak = vPeak[0];
    rfc->setAmplitude(1.0e6 * peak);
    if (freq.size() > 0) frequency = 1.0e6 * two_pi * freq[0];
    rfc->setFrequency(frequency);
    if (phase.size() > 0) phi = phase[0];
    rfc->setPhase(phi);

    rfc->dropFieldmaps();

    std::vector<double>::iterator peak_it;
    for (peak_it = vPeak.begin(); peak_it != vPeak.end(); ++ peak_it) {
        rfc->setAmplitudem(*peak_it);
    }
    std::vector<double>::iterator freq_it;
    for (freq_it = freq.begin(); freq_it != freq.end(); ++ freq_it) {
        rfc->setFrequencym(1.0e6 * two_pi * (*freq_it));
    }
    std::vector<double>::iterator phase_it;
    for (phase_it = phase.begin(); phase_it != phase.end(); ++ phase_it) {
        rfc->setPhasem(*phase_it);
    }
    std::vector<string>::iterator fmap_it;
    for (fmap_it = fmapfns.begin(); fmap_it != fmapfns.end(); ++ fmap_it) {
        rfc->setFieldMapFN(*fmap_it);
    }

    rfc->setFast(fast);
    rfc->setAutophaseVeto(apVeto);
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
