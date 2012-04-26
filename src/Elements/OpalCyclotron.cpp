// ------------------------------------------------------------------------
// $RCSfile: OpalCyclotron.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalCyclotron
//   The class of OPAL cyclotron.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalCyclotron.h"
#include "AbstractObjects/Attribute.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/CyclotronRep.h"
#include "Physics/Physics.h"


// Class OpalCyclotron
// ------------------------------------------------------------------------

OpalCyclotron::OpalCyclotron():
    OpalElement(SIZE, "CYCLOTRON",
                "The \"CYCLOTRON\" defines an cyclotron") {
    itsAttr[CYHARMON] = Attributes::makeReal
                        ("CYHARMON", "the harmonic number of the cyclotron");
    itsAttr[SYMMETRY] = Attributes::makeReal
                        ("SYMMETRY", "defines how the field is stored");
    itsAttr[RINIT] = Attributes::makeReal
                     ("RINIT", "Initial radius of the reference particle [m]");
    itsAttr[PRINIT] = Attributes::makeReal
                      ("PRINIT", "Initial radial momentum of the reference particle, pr=beta_r*gamma");
    itsAttr[PHIINIT] = Attributes::makeReal
                       ("PHIINIT", "Initial azimuth of the reference particle [deg]");
    itsAttr[FMAPFN] = Attributes::makeString
                      ("FMAPFN", "Filename for the B fieldmap");
    itsAttr[BSCALE] = Attributes::makeReal
                      ("BSCALE", "Scale factor for the B-field", 1.0);

    itsAttr[RFFREQ] = Attributes::makeRealArray
                      ("RFFREQ", "RF Frequency(ies) [MHz]");

    itsAttr[ESCALE] = Attributes::makeRealArray
                      ("ESCALE", "Scale factor for the RF field(s)");
    
    itsAttr[SUPERPOSE] = Attributes::makeBool
 	                 ("SUPERPOSE", "Option Whether all of the electric field maps are superposed, only used when TYPE = BANDRF", true);

    itsAttr[RFMAPFN] = Attributes::makeStringArray
                       ("RFMAPFN", "Filename for the RF fieldmap(s)");

    itsAttr[RFPHI] = Attributes::makeRealArray
                     ("RFPHI", "Initial phase(s) of the electric field map(s) [deg]");

    itsAttr[TYPE] = Attributes::makeString
                    ("TYPE", "Used to identify special cyclotron types");    
    itsAttr[TCR1] = Attributes::makeReal
                    ("TCR1", "trim coil r1 [mm]");
    itsAttr[TCR2] = Attributes::makeReal
                    ("TCR2", "trim coil r2 [mm]");
    itsAttr[MBTC] = Attributes::makeReal
                    ("MBTC", "max bfield of trim coil [kG]");
    itsAttr[SLPTC] = Attributes::makeReal
                    ("SLPTC", "slope of the rising edge");

    itsAttr[MINZ] = Attributes::makeReal
                    ("MINZ","Minimal vertical extend of the machine [mm]",-10000.0);
    itsAttr[MAXZ] = Attributes::makeReal
                    ("MAXZ","Maximal vertical extend of the machine [mm]",10000.0);
    itsAttr[MINR] = Attributes::makeReal
                    ("MINR","Minimal radial extend of the machine [mm]", 10.0);
    itsAttr[MAXR] = Attributes::makeReal
                   ("MAXR","Maximal radial extend of the machine [mm]", 10000.0);
    
    registerStringAttribute("FMAPFN");
    registerStringAttribute("RFMAPFN");
    registerStringAttribute("TYPE");
    registerRealAttribute("CYHARMON");
    registerRealAttribute("RINIT");
    registerRealAttribute("PRINIT");
    registerRealAttribute("PHIINIT");
    registerRealAttribute("SYMMETRY");
    registerRealAttribute("RFFREQ");
    registerRealAttribute("BSCALE");
    registerRealAttribute("ESCALE");
    registerRealAttribute("TCR1");
    registerRealAttribute("TCR2");
    registerRealAttribute("MBTC");
    registerRealAttribute("SLPTC");
    registerRealAttribute("RFPHI");
    registerRealAttribute("MINZ");
    registerRealAttribute("MAXZ");
    registerRealAttribute("MINR");
    registerRealAttribute("MAXR");

    setElement((new CyclotronRep("CYCLOTRON"))->makeAlignWrapper());
}


OpalCyclotron::OpalCyclotron(const string &name, OpalCyclotron *parent):
    OpalElement(name, parent) {
    setElement((new CyclotronRep(name))->makeAlignWrapper());
}


OpalCyclotron::~OpalCyclotron()
{}


OpalCyclotron *OpalCyclotron::clone(const string &name) {
    return new OpalCyclotron(name, this);
}


void OpalCyclotron::
fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);
}


void OpalCyclotron::update() {
    using Physics::two_pi;
    CyclotronRep *cycl =
        dynamic_cast<CyclotronRep *>(getElement()->removeWrappers());

    string fmapfm = Attributes::getString(itsAttr[FMAPFN]);

    string type = Attributes::getString(itsAttr[TYPE]);

    double harmnum = Attributes::getReal(itsAttr[CYHARMON]);
    double symmetry = Attributes::getReal(itsAttr[SYMMETRY]);
    double rinit = Attributes::getReal(itsAttr[RINIT]);
    double prinit = Attributes::getReal(itsAttr[PRINIT]);
    double phiinit = Attributes::getReal(itsAttr[PHIINIT]);
    double bscale = Attributes::getReal(itsAttr[BSCALE]);

    double tcr1 = Attributes::getReal(itsAttr[TCR1]);
    double tcr2 = Attributes::getReal(itsAttr[TCR2]);
    double mbtc = Attributes::getReal(itsAttr[MBTC]);
    double slptc = Attributes::getReal(itsAttr[SLPTC]);
    bool   superpose = Attributes::getBool(itsAttr[SUPERPOSE]);
    
    double minz = Attributes::getReal(itsAttr[MINZ]);
    double maxz = Attributes::getReal(itsAttr[MAXZ]);
    double minr = Attributes::getReal(itsAttr[MINR]);
    double maxr = Attributes::getReal(itsAttr[MAXR]);

    cycl->setFieldMapFN(fmapfm);
    cycl->setSymmetry(symmetry);

    cycl->setRinit(rinit);
    cycl->setPRinit(prinit);
    cycl->setPHIinit(phiinit);
    cycl->setBScale(bscale);

    cycl->setType(type);
    cycl->setCyclHarm(harmnum);

    cycl->setTCr1(tcr1);
    cycl->setTCr2(tcr2);
    cycl->setMBtc(mbtc);
    cycl->setSLPtc(slptc);
    cycl->setSuperpose(superpose);


    cycl->setMinR(minr);
    cycl->setMaxR(maxr);
    cycl->setMinZ(minz);
    cycl->setMaxZ(maxz);

    std::vector<string> fm_str = Attributes::getStringArray(itsAttr[RFMAPFN]);
    std::vector<double> scale_str = Attributes::getRealArray(itsAttr[ESCALE]);
    std::vector<double> phi_str = Attributes::getRealArray(itsAttr[RFPHI]);
    std::vector<double> rff_str = Attributes::getRealArray(itsAttr[RFFREQ]);

    // if ((fm_str.size() == scale_str.size()) && 
    //     (fm_str.size() == phi_str.size()) && 
    //     (fm_str.size() == rff_str.size())) {

    //     std::vector<string>::const_iterator fm    = fm_str.begin();
    //     std::vector<double>::const_iterator scale = scale_str.begin();        
    //     std::vector<double>::const_iterator phi   = phi_str.begin();
    //     std::vector<double>::const_iterator rff   = rff_str.begin();
        
    //     for( ; fm != fm_str.end(); ++fm,++scale,++phi,++rff)
    //         INFOMSG(" FM= " << *fm << "\t SCALE= " << *scale << "\t PHI= " << *phi << "\t FREQU= " << *rff << endl;);
    // }

    cycl->setRfPhi(phi_str);
    cycl->setEScale(scale_str);
    cycl->setRfFieldMapFN(fm_str);
    cycl->setRfFrequ(rff_str);
    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(cycl);
}

//  LocalWords:  rfphi
