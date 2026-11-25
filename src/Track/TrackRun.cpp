// Class TrackRun
//   The RUN command.
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
#include "Track/TrackRun.h"

#include "AbstractObjects/OpalData.h"
#include "Algorithms/ParallelTracker.h"
#include "Attributes/Attributes.h"
#include "Track/BuildElements.h"
#include "Track/BuildTracking.h"
#include "Track/FreezeRun.h"
#include "Track/RunDescriptionWriter.h"
#include "Track/TrackCmd.h"
#include "Track/ValidateRun.h"
#include "Utilities/OpalException.h"
#include "Utility/Inform.h"

#include <fstream>
#include <memory>
#include <vector>

extern Inform* gmsg;

namespace TRACKRUN {
    // The attributes of class TrackRun.
    enum {
        METHOD,            // Tracking method to use.
        FIELDSOLVER,       // The field solver attached
        BOUNDARYGEOMETRY,  // The boundary geometry
        TRACKBACK,         // In case we run the beam backwards
        SCFIELDUPDATE,     // Time within a drift-kick-drift step at which self-fields are solved.
        SIZE
    };
}  // namespace TRACKRUN

TrackRun::TrackRun()
    : Action(TRACKRUN::SIZE, "RUN",
             "The \"RUN\" sub-command tracks the defined particles through "
             "the given lattice.") {
    itsAttr[TRACKRUN::METHOD] = Attributes::makePredefinedString(
            "METHOD", "Name of tracking algorithm to use.", {"PARALLEL"});

    itsAttr[TRACKRUN::SCFIELDUPDATE] = Attributes::makePredefinedString(
            "SCFIELDUPDATE",
            "Space-charge field evaluation point: MIDPOINT uses positions after the first half "
            "drift; PRESTEP reproduces historical OPAL by solving before that drift.",
            {"MIDPOINT", "PRESTEP"}, "MIDPOINT");

    itsAttr[TRACKRUN::FIELDSOLVER] =
            Attributes::makeString("FIELDSOLVER", "Field solver to be used.");

    itsAttr[TRACKRUN::BOUNDARYGEOMETRY] = Attributes::makeString(
            "BOUNDARYGEOMETRY", "Boundary geometry to be used NONE (default).", "NONE");

    itsAttr[TRACKRUN::TRACKBACK] =
            Attributes::makeBool("TRACKBACK", "Track in reverse direction, default: false.", false);

    registerOwnership(AttributeHandler::SUB_COMMAND);
}

TrackRun::~TrackRun() {}

void TrackRun::execute() {
    const opalx::run::TrackDescription* track = TrackCmd::openTrack();
    if (track == nullptr) {
        throw OpalException("TrackRun::execute", "RUN is only valid inside TRACK.");
    }

    OpalData* opal = OpalData::getInstance();

    // build run description
    opalx::run::RunDescription description =
            opalx::run::freeze(*this, *track, ++opal->numberOfRuns);
    
    // validate parsed description
    std::vector<std::string> errors = opalx::run::validate(description);

    // build runtime elements
    std::vector<std::shared_ptr<ElementBase>> elements =
            opalx::run::buildElements(description.line, description.context.p0, errors);

    // --check run
    if (opal->checkOnly) {
        const std::string fileName = opalx::run::dumpFileName(description);
        if (ippl::Comm->rank() == 0) {
            std::ofstream out(fileName);
            opalx::run::write(description, out);
        }
        *gmsg << level1 << "* --check: wrote the run description to " << fileName
              << "; RUN does not track." << endl;
    }

    // most input errors caught here
    if (!errors.empty()) {
        throw OpalException("TrackRun::execute", opalx::run::joinErrors(errors));
    }

    // --check run ends here
    if (opal->checkOnly){
        return;
    } 

    /**
     * @brief Builds tracking objects
     * Tracking objects:
     *  - PartBunch
     *  - Containers
     *  - Solver
     *  - Beamline
     *  - Global Processes: decay, ...
     */
    opalx::run::RunObjects objects;
    TrackingSetup setup = opalx::run::buildTracking(description, std::move(elements), objects);

    // Build tracker
    auto tracker = std::make_unique<ParallelTracker>(std::move(setup));

    // Run tracker -> Start simulation
    tracker->execute();
}
