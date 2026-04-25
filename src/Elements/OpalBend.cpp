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

#include "Attributes/Attributes.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"

#include <cmath>
#include <iostream>
#include <limits>

// Class OpalBend
// ------------------------------------------------------------------------

OpalBend::OpalBend(const char* name, const char* help) : OpalElement(SIZE, name, help) {
    itsAttr[ANGLE] = Attributes::makeReal("ANGLE", "Total analytic bend angle in rad");
    itsAttr[K0]    = Attributes::makeReal(
            "K0", "Deprecated normal dipole coefficient in m^(-1), checked against ANGLE / L");
    itsAttr[K0S]      = Attributes::makeReal("K0S", "Skew dipole coefficient in m^(-1)");
    itsAttr[K1]       = Attributes::makeReal("K1", "Upright quadrupole coefficient in m^(-2)");
    itsAttr[K1S]      = Attributes::makeReal("K1S", "Skew quadrupole coefficient in m^(-2)");
    itsAttr[K2]       = Attributes::makeReal("K2", "Upright sextupole coefficient in m^(-3)");
    itsAttr[K2S]      = Attributes::makeReal("K2S", "Skew sextupole coefficient in m^(-3)");
    itsAttr[K3]       = Attributes::makeReal("K3", "Upright octupole coefficient in m^(-4)");
    itsAttr[K3S]      = Attributes::makeReal("K3S", "Skew octupole coefficient in m^(-4)");
    itsAttr[E1]       = Attributes::makeReal("E1", "Entry pole face angle in rad", 0.0);
    itsAttr[E2]       = Attributes::makeReal("E2", "Exit pole face angle in rad", 0.0);
    itsAttr[H1]       = Attributes::makeReal("H1", "Entry pole face curvature in m^(-1)");
    itsAttr[H2]       = Attributes::makeReal("H2", "Exit pole face curvature in m^(-1)");
    itsAttr[HGAP]     = Attributes::makeReal("HGAP", "Half gap width m");
    itsAttr[FINT]     = Attributes::makeReal("FINT", "Field integral (no dimension)", 0.5);
    itsAttr[SLICES]   = Attributes::makeReal("SLICES", "Number of slices to use", 1.0);
    itsAttr[STEPSIZE] = Attributes::makeReal("STEPSIZE", "Step-size to use for integration");
    itsAttr[FMAPFN]   = Attributes::makeString("FMAPFN", "Filename for the fieldmap");
    itsAttr[GAP]      = Attributes::makeReal("GAP", "Full gap height of the magnet (m)", 0.0);
    itsAttr[HAPERT]   = Attributes::makeReal("HAPERT", "Bend plane magnet aperture (m)", 0.0);
    itsAttr[ROTATION] =
            Attributes::makeReal("ROTATION", "-- not supported any more; use PSI instead --");
    itsAttr[DESIGNENERGY] =
            Attributes::makeReal("DESIGNENERGY", "the mean energy of the particles in MeV");
    itsAttr[GREATERTHANPI] = Attributes::makeBool("GREATERTHANPI", "-- not supported any more --");
    itsAttr[NSLICES]       = Attributes::makeReal(
            "NSLICES", "The number of slices/ steps for this element in Map Tracking", 1);
}

OpalBend::OpalBend(const std::string& name, OpalBend* parent) : OpalElement(name, parent) {}

OpalBend::~OpalBend() {}

void OpalBend::print(std::ostream& os) const { OpalElement::print(os); }

double OpalBend::deriveAnalyticDipoleCoefficient(const double arcLength, const double angle) {
    if (std::abs(arcLength) <= std::numeric_limits<double>::epsilon()) {
        return angle;
    }

    return angle / arcLength;
}

void OpalBend::validateAnalyticBendDefinition(
        const std::string& elementName, const bool hasAngle, const bool hasK0,
        const double arcLength, const double angle, const double k0Input) {
    if (!hasAngle) {
        throw OpalException(
                elementName,
                "ANGLE is mandatory for analytic bend geometry. K0 is only supported as a "
                "consistency-check quantity.");
    }

    if (!hasK0) {
        return;
    }

    if (arcLength <= 0.0) {
        throw OpalException(
                elementName,
                "K0 consistency checking requires a positive LENGTH so that K0 can be compared "
                "against ANGLE / L.");
    }

    const double derivedK0 = deriveAnalyticDipoleCoefficient(arcLength, angle);
    const double tolerance = 1.0e-12 + 1.0e-9 * std::max(std::abs(derivedK0), std::abs(k0Input));
    if (std::abs(k0Input - derivedK0) > tolerance) {
        throw OpalException(
                elementName,
                "Inconsistent analytic bend definition: K0 must satisfy K0 = ANGLE / L within "
                "tolerance for SBEND/RBEND without fringe fields.");
    }
}
