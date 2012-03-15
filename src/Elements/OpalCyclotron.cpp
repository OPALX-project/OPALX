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
                "The \"CYCLOTRON\" defines an cyclotron")
{
    itsAttr[CYHARMON] = Attributes::makeReal
        ("CYHARMON", "the harmonic number of the cyclotron");
    itsAttr[SYMMETRY] = Attributes::makeReal
        ("SYMMETRY", "defines how the field is stored");
    itsAttr[RINIT] = Attributes::makeReal
        ("RINIT", "Initial radius [m]");
    itsAttr[PRINIT] = Attributes::makeReal
        ("PRINIT", "Initial radial momenta [pr/p0] []");
    itsAttr[PHIINIT] = Attributes::makeReal
        ("PHIINIT", "the initial phase [deg]");
    itsAttr[RFFREQ] = Attributes::makeReal
        ("RFFREQ", "First hamonic of the RF system");
    itsAttr[BSCALE] = Attributes::makeReal
        ("BSCALE", "Scale factor for the B-field",1.0);
    itsAttr[FMAPFN] = Attributes::makeString
        ("FMAPFN", "Filename for the fieldmap");

    itsAttr[TYPE] = Attributes::makeString
        ("TYPE", "Used to identify special cyclotron types");
  
  
    itsAttr[TCR1] = Attributes::makeReal
        ("TCR1", "trim coil r1 (mm)"); 
   
    itsAttr[TCR2] = Attributes::makeReal
        ("TCR2", "trim coil r2 (mm)");
    itsAttr[MBTC] = Attributes::makeReal
        ("MBTC", "max bfield of trim coil (kG)");
    itsAttr[SLPTC] = Attributes::makeReal
    ("SLPTC", "slope of the rising edge"); 
  
 

    registerStringAttribute("FMAPFN");
    registerStringAttribute("TYPE");
    registerRealAttribute("CYHARMON");
    registerRealAttribute("RINIT");
    registerRealAttribute("PRINIT");
    registerRealAttribute("PHIINIT");
    registerRealAttribute("SYMMETRY");
    registerRealAttribute("RFFREQ");
    registerRealAttribute("BSCALE");
    registerRealAttribute("TCR1");
    registerRealAttribute("TCR2");
    registerRealAttribute("MBTC");
    registerRealAttribute("SLPTC");
    setElement((new CyclotronRep("CYCLOTRON"))->makeAlignWrapper());
}


OpalCyclotron::OpalCyclotron(const string &name, OpalCyclotron *parent):
    OpalElement(name, parent)
{
    setElement((new CyclotronRep(name))->makeAlignWrapper());
}


OpalCyclotron::~OpalCyclotron()
{}


OpalCyclotron *OpalCyclotron::clone(const string &name)
{
    return new OpalCyclotron(name, this);
}


void OpalCyclotron::
fillRegisteredAttributes(const ElementBase &base, ValueFlag flag)
{
    OpalElement::fillRegisteredAttributes(base, flag);
  
    if (flag != ERROR_FLAG) {
        const CyclotronRep *cycl =
            dynamic_cast<const CyclotronRep *>(base.removeWrappers());

        //    attributeRegistry["RINIT"]->setReal(Attributes::getReal(itsAttr[RINIT]));
        //attributeRegistry["PHIINIT"]->setReal(cycl->getPhase());
        //attributeRegistry["FMAPFN"]->setString(cycl->getFieldMapFN());
        //attributeRegistry["RFFREQ"]->setReal(cycl->getRfFrequ());
        //*gmsg << " attr rinir = " << Attributes::getReal(itsAttr[RINIT]) << endl;

    }
}


void OpalCyclotron::update()
{
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
    double rffrequ = Attributes::getReal(itsAttr[RFFREQ]);
    double bscale = Attributes::getReal(itsAttr[BSCALE]);


    double tcr1 = Attributes::getReal(itsAttr[TCR1]);
    double tcr2 = Attributes::getReal(itsAttr[TCR2]);
    double mbtc = Attributes::getReal(itsAttr[MBTC]);
    double slptc = Attributes::getReal(itsAttr[SLPTC]);
    cycl->setFieldMapFN(fmapfm);
    cycl->setSymmetry(symmetry);
    cycl->setRfFrequ(rffrequ);

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

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(cycl);
}

//  LocalWords:  OpalCyclotron
//  attributeRegistry["RINIT"]->setReal(Attributes::getReal(itsAttr[RINIT]));
