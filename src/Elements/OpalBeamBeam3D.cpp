// ------------------------------------------------------------------------
// $RCSfile: OpalBeamBeam3D.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalBeamBeam3D
//   The class for OPAL BEAMINT elements.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:38 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalBeamBeam3D.h"
#include "Algorithms/BeamBeam3D.h"
#include "Attributes/Attributes.h"
#include "Physics/Physics.h"
#include "Utilities/Round.h"
#include <cmath>


// Class OpalBeamBeam3D
// ------------------------------------------------------------------------

namespace {
    // Auxiliary inline function to square a double variable.
    inline double sqr(double x) { return x*x; }
}


OpalBeamBeam3D::OpalBeamBeam3D():
    OpalElement(SIZE, "BEAMINT",
                "The \"BEAMINT\" element defines a 3-dimensional beam-beam "
                "interaction point.")
{
    // Data for 4-D element.
    itsAttr[SIGX] =
        Attributes::makeReal
        ("SIGX", "Horizontal standard deviation of strong beam.");
    itsAttr[SIGY] =
        Attributes::makeReal
        ("SIGY", "Vertical standard deviation of strong beam.");
    itsAttr[XMA] =
        Attributes::makeReal
        ("XMA", "Horizontal displacement of strong beam.");
    itsAttr[YMA] =
        Attributes::makeReal
        ("YMA", "Vertical displacement of strong beam.");
    itsAttr[ZMA] =
        Attributes::makeReal
        ("ZMA", "Longitudinal displacement of strong beam.");

    registerRealAttribute("SIGX");
    registerRealAttribute("SIGY");
    registerRealAttribute("XMA");
    registerRealAttribute("YMA");
    registerRealAttribute("ZMA");

    // Data for 6-D element.
    itsAttr[SLICES] =
        Attributes::makeReal
        ("SLICES", "Number of slices in strong beam.", 1.0);
    itsAttr[ANGLE] =
        Attributes::makeReal
        ("ANGLE", "Horizontal crossing angle.");
    itsAttr[XIYN] =
        Attributes::makeReal
        ("XIYN", "Tune shift factor.");
    itsAttr[FAST] = Attributes::makeBool
        ("FAST", "If true, use tables for error function.");
    itsAttr[ALFXS] = Attributes::makeReal
        ("ALFXS", "Horizontal alpha* for the strong beam.");
    itsAttr[ALFYS] = Attributes::makeReal
        ("ALFYS", "Vertical alpha* for the strong beam.");
    itsAttr[DXS] = Attributes::makeReal
        ("DXS", "Horizontal dispersion for strong beam.");
    itsAttr[DPXS] = Attributes::makeReal
        ("DPXS", "Derivative of horizontal dispersion for strong beam.");
    itsAttr[DYS] = Attributes::makeReal
        ("DYS", "Vertical dispersion for strong beam.");
    itsAttr[DPYS] = Attributes::makeReal
        ("DPYS", "Derivative of vertical dispersion for strong beam.");
    itsAttr[EXS] = Attributes::makeReal
        ("EXS", "Horizontal emittance for strong beam in m.");
    itsAttr[EYS] = Attributes::makeReal
        ("EYS", "Vertical emittance for strong beam in m.");
    itsAttr[SIGTS] = Attributes::makeReal
        ("SIGTS", "Bunch length of strong beam.");
    itsAttr[SIGES] = Attributes::makeReal
        ("SIGES", "Energy spread of strong beam.");

    setElement(new BeamBeam3D("BEAMINT"));
}


OpalBeamBeam3D::OpalBeamBeam3D(const string &name, OpalBeamBeam3D *parent):
    OpalElement(name, parent)
{
    setElement(new BeamBeam3D(name));
}


OpalBeamBeam3D::~OpalBeamBeam3D()
{}


OpalBeamBeam3D *OpalBeamBeam3D::clone(const string &name)
{
    return new OpalBeamBeam3D(name, this);
}


void OpalBeamBeam3D::
fillRegisteredAttributes(const ElementBase &base, ValueFlag flag)
{
    OpalElement::fillRegisteredAttributes(base, flag);
    // *** MISSING: BeamBeam3D
}


void OpalBeamBeam3D::update()
{BeamBeam3D *bb = dynamic_cast<BeamBeam3D *>(getElement());

    bb->setElementLength(0.0);

    // Crossing angle.
    bb->setCrossingAngle(Attributes::getReal(itsAttr[ANGLE]));

    // Beam displacement.
    Vector3D disp(Attributes::getReal(itsAttr[XMA]),
                  Attributes::getReal(itsAttr[YMA]),
                  Attributes::getReal(itsAttr[ZMA]));

    // Lattice functions for strong beam.
    BeamBeam3D::Beta beta;
    beta.emitx = Attributes::getReal(itsAttr[EXS]);
    if (beta.emitx == 0.0) beta.emitx = 1.0;
    beta.emity = Attributes::getReal(itsAttr[EYS]);
    if (beta.emity == 0.0) beta.emity = 1.0;
    beta.betax = pow(Attributes::getReal(itsAttr[SIGX]), 2) / beta.emitx;
    if (beta.betax == 0.0) beta.betax = 1.0;
    beta.betay = pow(Attributes::getReal(itsAttr[SIGY]), 2) / beta.emity;
    if (beta.betay == 0.0) beta.betay = 1.0;
    beta.alphax = Attributes::getReal(itsAttr[ALFXS]);
    beta.alphay = Attributes::getReal(itsAttr[ALFYS]);
    beta.etax  = Attributes::getReal(itsAttr[DXS]);
    beta.etapx = Attributes::getReal(itsAttr[DPXS]);
    beta.etay  = Attributes::getReal(itsAttr[DYS]);
    beta.etapy = Attributes::getReal(itsAttr[DPYS]);
    beta.sige  = Attributes::getReal(itsAttr[SIGES]);
    beta.sigt  = Attributes::getReal(itsAttr[SIGTS]);
    bb->setBeamDescription(disp, beta);

    // Beam-beam parameter and number of slices.
    bb->setBeamBeamParameter(Attributes::getReal(itsAttr[XIYN]));
    bb->setSlices(int(Round(Attributes::getReal(itsAttr[SLICES]))));

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(bb);
}
