//
// Class TrackCmd
//   The class for the OPAL TRACK command.
//
// Copyright (c) 200x - 2022, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "Track/TrackCmd.h"

#include <algorithm>
#include <optional>

#include "AbstractObjects/BeamSequence.h"
#include "Attributes/Attributes.h"
#include "Structure/Beam.h"
#include "Track/RunDescription.h"
#include "Utilities/OpalException.h"
#include "Utility/Inform.h"

extern Inform* gmsg;

namespace {
    // The attributes of class TrackCmd
    enum {
        LINE,      // The name of lattice to be tracked.
        SOURCES,   // The name of the emission sources list (EMISSIONSOURCELIST).
        BEAM,      // The name of beam to be used.
        BEAMS,     // The names of beams to be used (replaces BEAM for input).
        DT,        // The integration timestep in second.
                   // In case of the adaptive integrator, time step guideline for
                   // external field integration.
        DTSCINIT,  // Only for adaptive integrator: Initial time step for space charge integration.
        DTAU,      // Only for adaptive integrator: Alternative way to set accuracy of space
                   // charge integration. Has no direct interpretation like DTSCINIT, but lower
                   // means smaller steps and more accurate. If given, DTSCINIT is not used. Useful
                   // for continuing with same step size in follow-up tracks.
        T0,        // The elapsed time (sec) of the bunch
        MAXSTEPS,  // The maximum timesteps we integrate
        ZSTART,    // Defines a z-location [m] where the reference particle starts
        ZSTOP,     // Defines a z-location [m], after which the simulation stops when the last
                   // particles passes
        TIMEINTEGRATOR,  // the name of time integrator
        SIZE
    };

    /// The settings of the TRACK that has no ENDTRACK yet.
    std::optional<opalx::run::TrackDescription> openTrack_m;
}  // namespace

TrackCmd::TrackCmd() : Action(SIZE, "TRACK", "The \"TRACK\" command initiates tracking.") {
    itsAttr[LINE] = Attributes::makeString("LINE", "Name of lattice to be tracked.");

    itsAttr[SOURCES] = Attributes::makeString(
            "SOURCES", "Name of the emission sources list (EMISSIONSOURCELIST).");

    itsAttr[BEAM] = Attributes::makeString("BEAM", "Name of beam to be used.", "UNNAMED_BEAM");

    itsAttr[BEAMS] =
            Attributes::makeStringArray("BEAMS", "Names of beams to be used (replaces BEAM).");

    itsAttr[DT] = Attributes::makeRealArray("DT", "The integration timestep in [s].");

    itsAttr[DTSCINIT] = Attributes::makeReal(
            "DTSCINIT",
            "Only for adaptive integrator: Initial time step for space charge integration.", 1e-12);

    itsAttr[DTAU] = Attributes::makeReal(
            "DTAU",
            "Only for adaptive integrator: Alternative way to set accuracy of space integration.",
            -1.0);

    itsAttr[T0] = Attributes::makeReal("T0", "The elapsed time of the bunch in seconds", 0.0);

    itsAttr[MAXSTEPS] = Attributes::makeRealArray(
            "MAXSTEPS",
            "The maximum number of integration steps dt, should be larger ZSTOP/(beta*c average).");

    itsAttr[ZSTART] = Attributes::makeReal(
            "ZSTART", "Defines a z-location [m] where the reference particle starts.", 0.0);

    itsAttr[ZSTOP] = Attributes::makeRealArray(
            "ZSTOP",
            "Defines a z-location [m], after which the simulation stops when the last particles "
            "passes.");

    itsAttr[TIMEINTEGRATOR] = Attributes::makePredefinedString(
            "TIMEINTEGRATOR", "Name of time integrator to be used.",
            {"RK-4", "RK4", "LF-2", "LF2", "MTS"}, "RK4");

    registerOwnership(AttributeHandler::COMMAND);
    AttributeHandler::addAttributeOwner("TRACK", AttributeHandler::COMMAND, "RUN");
    AttributeHandler::addAttributeOwner("TRACK", AttributeHandler::COMMAND, "ENDTRACK");
}

TrackCmd::~TrackCmd() {}

std::vector<double> TrackCmd::getDT() const {
    std::vector<double> dTs = Attributes::getRealArray(itsAttr[DT]);
    if (dTs.size() == 0) {
        dTs.push_back(1e-12);
    }
    for (double dt : dTs) {
        if (dt < 0.0) {
            throw OpalException(
                    "TrackCmd::getDT", "The time steps provided with DT have to be positive");
        }
    }
    return dTs;
}

double TrackCmd::getT0() const { return Attributes::getReal(itsAttr[T0]); }

double TrackCmd::getZStart() const { return Attributes::getReal(itsAttr[ZSTART]); }

std::vector<double> TrackCmd::getZStop() const {
    std::vector<double> zstop = Attributes::getRealArray(itsAttr[ZSTOP]);
    if (zstop.size() == 0) {
        zstop.push_back(1000000.0);
    }
    return zstop;
}

std::vector<unsigned long long> TrackCmd::getMaxSteps() const {
    std::vector<double> maxsteps_d = Attributes::getRealArray(itsAttr[MAXSTEPS]);
    std::vector<unsigned long long> maxsteps_i;
    if (maxsteps_d.size() == 0) {
        maxsteps_i.push_back(10ul);
    }
    for (double numSteps : maxsteps_d) {
        if (numSteps < 0) {
            throw OpalException(
                    "TrackCmd::getMAXSTEPS",
                    "The number of steps provided with MAXSTEPS has to be positive");
        } else {
            unsigned long long value = numSteps;
            maxsteps_i.push_back(value);
        }
    }

    return maxsteps_i;
}

const opalx::run::TrackDescription* TrackCmd::openTrack() {
    return openTrack_m ? &*openTrack_m : nullptr;
}

void TrackCmd::closeTrack() {
    if (!openTrack_m) {
        throw OpalException("TrackCmd::closeTrack", "ENDTRACK without TRACK.");
    }
    openTrack_m.reset();
}

void TrackCmd::execute() {
    if (openTrack_m) {
        throw OpalException("TrackCmd::execute", "TRACK inside an open TRACK; add ENDTRACK first.");
    }

    opalx::run::TrackDescription track;
    track.line = Attributes::getString(itsAttr[LINE]);
    BeamSequence::find(track.line);  // fail fast

    // TRACK,SOURCES is optional (and ignored); emission sources come from the selected BEAM.

    // BEAMS has precedence over BEAM.
    track.beamNames = Attributes::getStringArray(itsAttr[BEAMS]);
    if (track.beamNames.empty()) {
        track.beamNames.push_back(Attributes::getString(itsAttr[BEAM]));
    }
    for (const auto& name : track.beamNames) {
        if (name.empty()) {
            throw OpalException("TrackCmd::execute", "Empty beam name in \"BEAMS\".");
        }
        Beam::find(name);  // fail fast
    }

    std::vector<double> dt                   = getDT();
    std::vector<unsigned long long> maxsteps = getMaxSteps();
    std::vector<double> zstop                = getZStop();

    // Accepted so that old input files read unchanged; nothing uses them.
    for (const int attr : {TIMEINTEGRATOR, DTSCINIT, DTAU}) {
        if (!itsAttr[attr].defaultUsed()) {
            *gmsg << level1 << "* TRACK " << itsAttr[attr].getName() << " has no effect in OPALX."
                  << endl;
        }
    }

    size_t numTracks = dt.size();
    numTracks        = std::max(numTracks, maxsteps.size());
    numTracks        = std::max(numTracks, zstop.size());
    for (size_t i = dt.size(); i < numTracks; ++i) {
        dt.push_back(dt.back());
    }
    for (size_t i = maxsteps.size(); i < numTracks; ++i) {
        maxsteps.push_back(maxsteps.back());
    }
    for (size_t i = zstop.size(); i < numTracks; ++i) {
        zstop.push_back(zstop.back());
    }

    track.dT         = dt;
    track.maxSteps   = maxsteps;
    track.zStop      = zstop;
    track.zStopGiven = !Attributes::getRealArray(itsAttr[ZSTOP]).empty();
    track.zStart     = getZStart();
    track.t0         = getT0();
    openTrack_m      = std::move(track);
}
