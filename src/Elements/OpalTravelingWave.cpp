//
// Class OpalTravelingWave
//   The TRAVELINGWAVE element.
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
#include "Elements/OpalTravelingWave.h"

#include <memory>
#include "AbstractObjects/Attribute.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/TravelingWaveRep.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"

OpalTravelingWave::OpalTravelingWave()
    : OpalElement(
              SIZE, "TRAVELINGWAVE",
              "The \"TRAVELINGWAVE\" element defines a traveling wave structure.") {
    itsAttr[VOLT]   = Attributes::makeReal("VOLT", "RF voltage in MV/m");
    itsAttr[DVOLT]  = Attributes::makeReal("DVOLT", "RF voltage error in MV/m");
    itsAttr[FREQ]   = Attributes::makeReal("FREQ", "RF frequency in MHz");
    itsAttr[LAG]    = Attributes::makeReal("LAG", "Phase lag in rad");
    itsAttr[DLAG]   = Attributes::makeReal("DLAG", "Phase lag error in rad");
    itsAttr[FMAPFN] = Attributes::makeString("FMAPFN", "Filename for the fieldmap");
    itsAttr[FAST]   = Attributes::makeBool("FAST", "Faster but less accurate", true);
    itsAttr[APVETO] = Attributes::makeBool(
            "APVETO", "Do not use this cavity in the Autophase procedure", false);
    itsAttr[NUMCELLS] = Attributes::makeReal("NUMCELLS", "Number of cells in a TW structure");
    itsAttr[DESIGNENERGY] =
            Attributes::makeReal("DESIGNENERGY", "the mean energy of the particles at exit", -1.0);
    itsAttr[MODE] = Attributes::makeReal(
            "MODE", "The phase shift between neighboring cells in 2*pi", 1.0 / 3.0);

    registerOwnership();

    setElement(new TravelingWaveRep("TRAVELINGWAVE"));
}

OpalTravelingWave::OpalTravelingWave(const std::string& name, OpalTravelingWave* parent)
    : OpalElement(name, parent) {
    setElement(new TravelingWaveRep(name));
}

OpalTravelingWave::~OpalTravelingWave() {}

OpalTravelingWave* OpalTravelingWave::clone(const std::string& name) {
    return new OpalTravelingWave(name, this);
}

std::shared_ptr<ElementBase> OpalTravelingWave::build(const ElementValues& values) {
    auto element = std::make_shared<TravelingWaveRep>(values.getName());
    fill(*element, values);
    return element;
}

void OpalTravelingWave::update() { fill(*getElement(), getValues()); }

void OpalTravelingWave::fill(ElementBase& element, const ElementValues& values) {
    setCommonFields(element, values);

    TravelingWaveRep* rfc = dynamic_cast<TravelingWaveRep*>(&element);

    double length      = values.getReal("L");
    double vPeak       = values.getReal("VOLT");
    double vPeakError  = values.getReal("DVOLT");
    double phase       = values.getReal("LAG");
    double phaseError  = values.getReal("DLAG");
    double freq        = Physics::two_pi * values.getReal("FREQ") * Units::MHz2Hz;
    std::string fmapfm = values.getString("FMAPFN");
    bool fast          = values.getBool("FAST");
    bool apVeto        = values.getBool("APVETO");

    //    std::string type = values.getString("TYPE");
    double kineticEnergy = values.getReal("DESIGNENERGY");

    rfc->getGeometry().setElementLength(length);

    rfc->setFieldMapFN(fmapfm);
    rfc->setFast(fast);
    rfc->setAutophaseVeto(apVeto);
    rfc->setAmplitudem(vPeak);
    rfc->setAmplitudeError(vPeakError);
    rfc->setFrequencym(freq);
    rfc->setPhasem(phase);
    rfc->setPhaseError(phaseError);
    rfc->setNumCells((int)values.getReal("NUMCELLS"));
    rfc->setMode(values.getReal("MODE"));
    rfc->setDesignEnergy(kineticEnergy);
}
