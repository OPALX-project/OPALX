// ------------------------------------------------------------------------
// $RCSfile: OpalBend.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalBend
//   The parent class of all OPAL bending magnets.
//   This class factors out all special behaviour for the DOOM interface
//   and the printing in OPAL format, as well as the bend attributes.
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:11 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalBend.h"
#include "AbstractObjects/DoomDB.h"
#include "AbstractObjects/DoomWriter.h"
#include "Attributes/Attributes.h"
#include "Utilities/Options.h"
#include <iostream>


// Class OpalBend
// ------------------------------------------------------------------------

OpalBend::OpalBend(const char *name, const char *help):
    OpalElement(SIZE, name, help) {
    itsAttr[ANGLE] = Attributes::makeReal
                     ("ANGLE", "Upright dipole coefficient in m^(-1)");
    itsAttr[K0] = Attributes::makeReal
                  ("K0", "Normal dipole coefficient in m^(-1)");
    itsAttr[K0S] = Attributes::makeReal
                   ("K0S", "Skew dipole coefficient in m^(-1)");
    itsAttr[K1] = Attributes::makeReal
                  ("K1", "Upright quadrupole coefficient in m^(-2)");
    itsAttr[K1S] = Attributes::makeReal
                   ("K1S", "Skew quadrupole coefficient in m^(-2)");
    itsAttr[K2] = Attributes::makeReal
                  ("K2", "Upright sextupole coefficient in m^(-3)");
    itsAttr[K2S] = Attributes::makeReal
                   ("K2S", "Skew sextupole coefficient in m^(-3)");
    itsAttr[K3] = Attributes::makeReal
                  ("K3", "Upright octupole coefficient in m^(-4)");
    itsAttr[K3S] = Attributes::makeReal
                   ("K3S", "Skew octupole coefficient in m^(-4)");
    itsAttr[E1] = Attributes::makeReal
                  ("E1", "Entry pole face angle in rad");
    itsAttr[E2] = Attributes::makeReal
                  ("E2", "Exit pole face angle in rad");
    itsAttr[H1] = Attributes::makeReal
                  ("H1", "Entry pole face curvature in m^(-1)");
    itsAttr[H2] = Attributes::makeReal
                  ("H2", "Exit pole face curvature in m^(-1)");
    itsAttr[HGAP] = Attributes::makeReal
                    ("HGAP", "Half gap width m");
    itsAttr[FINT] = Attributes::makeReal
                    ("FINT", "Field integral (no dimension)", 0.5);
    itsAttr[SLICES] = Attributes::makeReal
                      ("SLICES", "Number of slices to use", 1.0);
    itsAttr[STEPSIZE] = Attributes::makeReal
                        ("STEPSIZE", "Step-size to use for integration");
    itsAttr[FMAPFN] = Attributes::makeString
                      ("FMAPFN", "Filename for the fieldmap");
    itsAttr[ALPHA] = Attributes::makeReal
                     ("ALPHA", "Pole face angle in degree");
    itsAttr[BETA] = Attributes::makeReal
                    ("BETA", "Pole face angle in degree");
    itsAttr[EXITANGLE] = Attributes::makeReal
                         ("EXITANGLE", "Angle between the entry and exit face (SBEND only)");
    itsAttr[DESIGNENERGY] = Attributes::makeReal
                            ("DESIGNENERGY", "the mean energy of the particles");

    registerRealAttribute("ANGLE");
    registerRealAttribute("K0L");
    registerRealAttribute("K0SL");
    registerRealAttribute("K1L");
    registerRealAttribute("K1SL");
    registerRealAttribute("K2L");
    registerRealAttribute("K2SL");
    registerRealAttribute("K3L");
    registerRealAttribute("K3SL");
    registerRealAttribute("E1");
    registerRealAttribute("E2");
    registerRealAttribute("H1");
    registerRealAttribute("H2");
    registerRealAttribute("HGAP");
    registerRealAttribute("FINT");
    registerRealAttribute("SLICES");
    registerRealAttribute("STEPSIZE");
    registerStringAttribute("FMAPFN");
    registerRealAttribute("ALPHA");
    registerRealAttribute("BETA");
    registerRealAttribute("EXITANGLE");
    registerRealAttribute("DESIGNENERGY");
}


OpalBend::OpalBend(const string &name, OpalBend *parent):
    OpalElement(name, parent)
{}


OpalBend::~OpalBend()
{}


void OpalBend::doomPut(DoomWriter &writer) const {
    // Store the OPAL-9 data.
    OpalElement::doomPut(writer);

    // Store curvature.
    double length = Attributes::getReal(itsAttr[LENGTH]);
    double angle  = Attributes::getReal(itsAttr[ANGLE]);
    int index = DoomDB::getAttributeIndex("RHOINV");
    if(itsAttr[K0]) {
        writer.putReal(index, Attributes::getReal(itsAttr[K0]));
    } else {
        writer.putReal(index, length ? (angle / length) : 0.0);
    }
    writer.putInt(index, 1);
}


void OpalBend::print(std::ostream &os) const {

  OpalElement::print(os);

}
