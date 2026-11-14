//
// Class OpalCavity
//   The RFCAVITY element.
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#include "Elements/OpalCavity.h"

#include <memory>

#include "AbstractObjects/Attribute.h"
#include "Algorithms/AbstractTimeDependence.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/RFCavityRep.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"
#include "Utilities/OpalException.h"

OpalCavity::OpalCavity()
    : OpalElement(SIZE, "RFCAVITY", "The \"RFCAVITY\" element defines an RF cavity.") {
    itsAttr[TYPE] = Attributes::makePredefinedString(
            "TYPE",
            "RF cavity model: mapped STANDING (default), mapped SINGLEGAP, or uniform IDEAL.",
            {"STANDING", "SINGLEGAP", "IDEAL"});
    itsAttr[VOLT] = Attributes::makeReal(
            "VOLT", "Peak RF electric field in MV/m (STANDING/IDEAL); voltage in MV (SINGLEGAP)");
    itsAttr[DVOLT]    = Attributes::makeReal("DVOLT", "Peak RF electric field error in MV/m");
    itsAttr[FREQ]     = Attributes::makeReal("FREQ", "RF frequency in MHz");
    itsAttr[LAG]      = Attributes::makeReal("LAG", "Phase lag (rad)");
    itsAttr[DLAG]     = Attributes::makeReal("DLAG", "Phase lag error (rad)");
    itsAttr[FMAPFN]   = Attributes::makeString("FMAPFN", "Filename of the fieldmap");
    itsAttr[GEOMETRY] = Attributes::makeString("GEOMETRY", "BoundaryGeometry for Cavities");
    itsAttr[FAST]     = Attributes::makeBool("FAST", "Faster but less accurate", true);
    itsAttr[APVETO]   = Attributes::makeBool(
            "APVETO", "Do not use this cavity in the Autophase procedure", false);
    itsAttr[RMIN]  = Attributes::makeReal("RMIN", "SINGLEGAP lower local-x profile bound [m]");
    itsAttr[RMAX]  = Attributes::makeReal("RMAX", "SINGLEGAP upper local-x profile bound [m]");
    itsAttr[ANGLE] = Attributes::makeReal("ANGLE", "Azimuth position of a cyclotron cavity [deg]");
    itsAttr[PDIS]  = Attributes::makeReal(
            "PDIS", "Shift distance of cavity gap from center of cyclotron [mm]");
    itsAttr[GAPWIDTH] = Attributes::makeReal("GAPWIDTH", "SINGLEGAP transit-time width [m]");
    itsAttr[PHI0]     = Attributes::makeReal("PHI0", "Initial phase of cavity [deg]");
    itsAttr[DESIGNENERGY] =
            Attributes::makeReal("DESIGNENERGY", "the mean energy of the particles at exit", -1.0);
    // attibutes for timedependent values
    itsAttr[PHASE_MODEL] =
            Attributes::makeString("PHASE_MODEL", "The name of the phase time dependence model.");
    itsAttr[AMPLITUDE_MODEL] = Attributes::makeString(
            "AMPLITUDE_MODEL", "The name of the amplitude time dependence model.");
    itsAttr[FREQUENCY_MODEL] = Attributes::makeString(
            "FREQUENCY_MODEL", "The name of the frequency time dependence model.");

    registerOwnership();
}

OpalCavity::OpalCavity(const std::string& name, OpalCavity* parent) : OpalElement(name, parent) {}

OpalCavity::~OpalCavity() {}

OpalCavity* OpalCavity::clone(const std::string& name) { return new OpalCavity(name, this); }

std::shared_ptr<ElementBase> OpalCavity::build(const ElementValues& values) {
    auto element = std::make_shared<RFCavityRep>(values.getName());
    fill(*element, values);
    return element;
}

void OpalCavity::fill(ElementBase& element, const ElementValues& values) {
    setCommonFields(element, values);

    RFCavityRep* rfc = dynamic_cast<RFCavityRep*>(&element);

    double length      = values.getReal("L");
    double peak        = values.getReal("VOLT");
    double peakError   = values.getReal("DVOLT");
    double phase       = values.getReal("LAG");
    double phaseError  = values.getReal("DLAG");
    double freq        = Physics::two_pi * values.getReal("FREQ") * Units::MHz2Hz;
    std::string fmapfn = values.getString("FMAPFN");
    std::string type   = values.getString("TYPE");
    bool fast          = values.getBool("FAST");
    bool apVeto        = (values.getBool("APVETO"));

    double rmin          = values.getReal("RMIN");
    double rmax          = values.getReal("RMAX");
    double angle         = values.getReal("ANGLE");
    double pdis          = values.getReal("PDIS");
    double gapwidth      = values.getReal("GAPWIDTH");
    double phi0          = values.getReal("PHI0");
    double kineticEnergy = values.getReal("DESIGNENERGY");

    rfc->clearCyclotronGap();
    if (type == "SINGLEGAP" && !fmapfn.empty()) {
        if (angle != 0 || pdis != 0 || peakError != 0 || phase != 0 || phaseError != 0
            || !values.getString("PHASE_MODEL").empty()
            || !values.getString("AMPLITUDE_MODEL").empty()
            || !values.getString("FREQUENCY_MODEL").empty())
            throw OpalException(
                    "OpalCavity::build",
                    "SINGLEGAP uses SI geometry and X/Y/Z/THETA/PHI/PSI poses, PHI0, and static RF "
                    "parameters only.");
        rfc->configureCyclotronGap(
                fmapfn, rmin, rmax,
                CyclotronRFKick{peak * 1e6, freq, phi0 * Units::deg2rad, gapwidth, rmax - rmin});
    }

    rfc->getGeometry().setElementLength(length);

    rfc->dropFieldmaps();

    rfc->setAmplitudem(peak);
    rfc->setAmplitudeError(peakError);
    rfc->setFrequencym(freq);
    rfc->setPhasem(phase);
    rfc->setPhaseError(phaseError);
    rfc->setFieldMapFN(fmapfn);

    rfc->setFast(fast);
    rfc->setAutophaseVeto(apVeto || rfc->isCyclotronGap());
    rfc->setCavityType(type);
    rfc->setRmin(rmin);
    rfc->setRmax(rmax);
    rfc->setAzimuth(angle);
    rfc->setPerpenDistance(pdis);
    rfc->setGapWidth(gapwidth);
    rfc->setPhi0(phi0);
    rfc->setDesignEnergy(kineticEnergy);

    rfc->setPhaseModelName(values.getString("PHASE_MODEL"));
    rfc->setAmplitudeModelName(values.getString("AMPLITUDE_MODEL"));
    rfc->setFrequencyModelName(values.getString("FREQUENCY_MODEL"));

    rfc->validateIdeal();
}
