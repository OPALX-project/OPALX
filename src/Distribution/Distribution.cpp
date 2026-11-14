//
// Class Distribution
//   This class defines the initial beam that is injected or emitted into the simulation.
//
// Copyright (c) 2008 - 2022, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "Distribution/Distribution.h"

#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Distribution/DistributionSettings.h"
#include "Track/FreezeRun.h"
#include "Utilities/OpalException.h"

namespace DISTRIBUTION {
    enum {
        TYPE,
        FNAME,
        SIGMAX,
        SIGMAY,
        SIGMAZ,
        SIGMAPX,
        SIGMAPY,
        SIGMAPZ,
        CORR,
        CUTOFFPX,
        CUTOFFPY,
        CUTOFFPZ,
        CUTOFFX,
        CUTOFFY,
        CUTOFFLONG,
        CORRX,
        CORRY,
        CORRZ,
        CORRT,
        SIGMAT,
        TPULSEFWHM,
        TRISE,
        TFALL,
        FTOSCAMPLITUDE,
        FTOSCPERIODS,
        EMITTED,
        EMISSIONSTEPS,
        /// Optional per-distribution particle count (macroparticles).
        /// If <= 0, this distribution does not specify its own N and
        /// TrackRun/BEAM logic is used instead.
        NPARTDIST,
        SIZE
    };
}

Distribution::Distribution()
    : Definition(
              DISTRIBUTION::SIZE, "DISTRIBUTION",
              "The DISTRIBUTION statement defines data for the 6D particle distribution.") {
    itsAttr[DISTRIBUTION::TYPE] = Attributes::makePredefinedString(
            "TYPE", "Distribution type.",
            {"UNIFORM", "GAUSS", "MULTIVARIATEGAUSS", "FLATTOP", "OPALFLATTOP", "FROMFILE",
             "EMITTEDFROMFILE"});

    itsAttr[DISTRIBUTION::FNAME] =
            Attributes::makeString("FNAME", "File for reading in 6D particle coordinates.", "");

    // Parameters for defining an initial distribution.
    itsAttr[DISTRIBUTION::SIGMAX]  = Attributes::makeReal("SIGMAX", "SIGMAx (m)", 0.0);
    itsAttr[DISTRIBUTION::SIGMAY]  = Attributes::makeReal("SIGMAY", "SIGMAy (m)", 0.0);
    itsAttr[DISTRIBUTION::SIGMAZ]  = Attributes::makeReal("SIGMAZ", "SIGMAz (m)", 0.0);
    itsAttr[DISTRIBUTION::SIGMAPX] = Attributes::makeReal("SIGMAPX", "SIGMApx", 0.0);
    itsAttr[DISTRIBUTION::SIGMAPY] = Attributes::makeReal("SIGMAPY", "SIGMApy", 0.0);
    itsAttr[DISTRIBUTION::SIGMAPZ] = Attributes::makeReal("SIGMAPZ", "SIGMApz", 0.0);

    itsAttr[DISTRIBUTION::CORR] = Attributes::makeRealArray("CORR", "r correlation");

    itsAttr[DISTRIBUTION::CUTOFFPX] = Attributes::makeReal(
            "CUTOFFPX", "Distribution cutoff px dimension in units of sigma.", 3.0);
    itsAttr[DISTRIBUTION::CUTOFFPY] = Attributes::makeReal(
            "CUTOFFPY", "Distribution cutoff py dimension in units of sigma.", 3.0);
    itsAttr[DISTRIBUTION::CUTOFFPZ] = Attributes::makeReal(
            "CUTOFFPZ", "Distribution cutoff pz dimension in units of sigma.", 3.0);

    itsAttr[DISTRIBUTION::CUTOFFX] = Attributes::makeReal(
            "CUTOFFX", "Distribution cutoff x direction in units of sigma.", 3.0);
    itsAttr[DISTRIBUTION::CUTOFFY] = Attributes::makeReal(
            "CUTOFFY", "Distribution cutoff r direction in units of sigma.", 3.0);
    itsAttr[DISTRIBUTION::CUTOFFLONG] = Attributes::makeReal(
            "CUTOFFLONG", "Distribution cutoff z or t direction in units of sigma.", 3.0);

    itsAttr[DISTRIBUTION::CORRX] =
            Attributes::makeReal("CORRX", "x/px correlation, (R12 in transport notation).", 0.0);
    itsAttr[DISTRIBUTION::CORRY] =
            Attributes::makeReal("CORRY", "y/py correlation, (R34 in transport notation).", 0.0);
    itsAttr[DISTRIBUTION::CORRZ] =
            Attributes::makeReal("CORRZ", "z/pz correlation, (R56 in transport notation).", 0.0);
    itsAttr[DISTRIBUTION::CORRT] =
            Attributes::makeReal("CORRT", "t/pt correlation, (R56 in transport notation).", 0.0);

    itsAttr[DISTRIBUTION::SIGMAT] = Attributes::makeReal("SIGMAT", "SIGMAt (m)", 0.0);
    itsAttr[DISTRIBUTION::TPULSEFWHM] =
            Attributes::makeReal("TPULSEFWHM", "Pulse FWHM for emitted distribution.", 0.0);
    itsAttr[DISTRIBUTION::TRISE] =
            Attributes::makeReal("TRISE", "Rise time for emitted distribution.", 0.0);
    itsAttr[DISTRIBUTION::TFALL] =
            Attributes::makeReal("TFALL", "Fall time for emitted distribution.", 0.0);

    itsAttr[DISTRIBUTION::FTOSCAMPLITUDE] = Attributes::makeReal(
            "FTOSCAMPLITUDE",
            "Amplitude of oscillations superimposed "
            "on flat top portion of emitted GAUSS "
            "distribtuion (in percent of flat top "
            "amplitude)",
            0.0);

    itsAttr[DISTRIBUTION::FTOSCPERIODS] = Attributes::makeReal(
            "FTOSCPERIODS",
            "Number of oscillations superimposed on "
            "flat top portion of emitted GAUSS "
            "distribution",
            0.0);

    itsAttr[DISTRIBUTION::EMITTED] = Attributes::makeBool(
            "EMITTED",
            "Emitted beam, from cathode, as opposed to "
            "an injected beam.",
            false);

    itsAttr[DISTRIBUTION::EMISSIONSTEPS] = Attributes::makeReal(
            "EMISSIONSTEPS",
            "Number of OPAL-like time steps to use during OPALFLATTOP or EMITTEDFROMFILE "
            "emission.",
            100.0);

    itsAttr[DISTRIBUTION::NPARTDIST] = Attributes::makeReal(
            "NPARTDIST",
            "Number of macroparticles for this DISTRIBUTION. "
            "If 0 or negative, this distribution does not specify its own count.",
            0.0);

    registerOwnership(AttributeHandler::STATEMENT);
}

Distribution::~Distribution() {}

bool Distribution::canReplaceBy(Object* object) { return dynamic_cast<Distribution*>(object) != 0; }

void Distribution::execute() {
    // Check the values when the DISTRIBUTION statement is read; a RUN works them out again.
    static_cast<void>(DistributionSettings(ObjectValues(opalx::run::describeObject(*this)), 0.0));
}

Distribution* Distribution::find(const std::string& name) {
    Distribution* dist = dynamic_cast<Distribution*>(OpalData::getInstance()->find(name));

    if (dist == 0) {
        throw OpalException("Distribution::find()", "Distribution \"" + name + "\" not found.");
    }

    return dist;
}
