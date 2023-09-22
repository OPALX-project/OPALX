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

#include "AbstractObjects/BeamSequence.h"

#include "AbstractObjects/OpalData.h"

#include "Attributes/Attributes.h"

#include "Beamlines/TBeamline.h"

#include "BasicActions/Option.h"

#include "Distribution/Distribution.h"

#include "Physics/Physics.h"
#include "Physics/Units.h"

#include "Track/Track.h"

#include "Utilities/OpalException.h"

#include "Structure/Beam.h"
#include "Structure/BoundaryGeometry.h"
#include "Structure/DataSink.h"
#include "Structure/H5PartWrapper.h"
#include "Structure/H5PartWrapperForPT.h"

#include "OPALconfig.h"
#include "changes.h"

#include <boost/assign.hpp>

#include <cmath>
#include <fstream>
#include <iomanip>

extern Inform* gmsg;

namespace TRACKRUN {
    // The attributes of class TrackRun.
    enum {
        METHOD,            // Tracking method to use.
        TURNS,             // The number of turns to be tracked, we keep that for the moment
        BEAM,              // The beam to track
        FIELDSOLVER,       // The field solver attached
        BOUNDARYGEOMETRY,  // The boundary geometry
        DISTRIBUTION,      // The particle distribution
        TRACKBACK,         // In case we run the beam backwards
        SIZE
    };
}  // namespace TRACKRUN

const std::string TrackRun::defaultDistribution_m("DISTRIBUTION");

const boost::bimap<TrackRun::RunMethod, std::string> TrackRun::stringMethod_s =
    boost::assign::list_of<const boost::bimap<TrackRun::RunMethod, std::string>::relation>(
        RunMethod::PARALLEL, "PARALLEL");

TrackRun::TrackRun()
    : Action(
        TRACKRUN::SIZE, "RUN",
        "The \"RUN\" sub-command tracks the defined particles through "
        "the given lattice."),
      itsTracker_m(nullptr),
      dist_m(nullptr),
      fs_m(nullptr),
      ds_m(nullptr),
      phaseSpaceSink_m(nullptr),
      isFollowupTrack_m(false),
      method_m(RunMethod::NONE),
      macromass_m(0.0),
      macrocharge_m(0.0) {
    itsAttr[TRACKRUN::METHOD] = Attributes::makePredefinedString(
        "METHOD", "Name of tracking algorithm to use.", {"PARALLEL"});

    itsAttr[TRACKRUN::TURNS] = Attributes::makeReal(
        "TURNS",
        "Number of turns to be tracked; Number of neighboring bunches to be tracked in cyclotron.",
        1.0);

    itsAttr[TRACKRUN::BEAM] = Attributes::makeString("BEAM", "Name of beam.");

    itsAttr[TRACKRUN::FIELDSOLVER] =
        Attributes::makeString("FIELDSOLVER", "Field solver to be used.");

    itsAttr[TRACKRUN::BOUNDARYGEOMETRY] = Attributes::makeString(
        "BOUNDARYGEOMETRY", "Boundary geometry to be used NONE (default).", "NONE");

    itsAttr[TRACKRUN::DISTRIBUTION] =
        Attributes::makeStringArray("DISTRIBUTION", "List of particle distributions to be used.");

    itsAttr[TRACKRUN::TRACKBACK] =
        Attributes::makeBool("TRACKBACK", "Track in reverse direction, default: false.", false);

    registerOwnership(AttributeHandler::SUB_COMMAND);
    opal = OpalData::getInstance();
}

TrackRun::TrackRun(const std::string& name, TrackRun* parent)
    : Action(name, parent),
      itsTracker_m(nullptr),
      dist_m(nullptr),
      fs_m(nullptr),
      ds_m(nullptr),
      phaseSpaceSink_m(nullptr),
      isFollowupTrack_m(false),
      method_m(RunMethod::NONE),
      macromass_m(0.0),
      macrocharge_m(0.0) {
    opal = OpalData::getInstance();
}

TrackRun::~TrackRun() {
    delete phaseSpaceSink_m;
}

TrackRun* TrackRun::clone(const std::string& name) {
    return new TrackRun(name, this);
}

void TrackRun::execute() {
    const int currentVersion = ((OPAL_VERSION_MAJOR * 100) + OPAL_VERSION_MINOR) * 100;
    if (Options::version < currentVersion) {
        unsigned int fileVersion = Options::version / 100;
        bool newerChanges        = false;
        for (auto it = Versions::changes.begin(); it != Versions::changes.end(); ++it) {
            if (it->first > fileVersion) {
                newerChanges = true;
                break;
            }
        }
        if (newerChanges) {
            Inform errorMsg("Error");
            errorMsg << "\n******************** V E R S I O N   M I S M A T C H "
                        "***********************\n"
                     << endl;
            for (auto it = Versions::changes.begin(); it != Versions::changes.end(); ++it) {
                if (it->first > fileVersion) {
                    errorMsg << it->second << endl;
                }
            }
            errorMsg << "\n"
                     << "* Make sure you do understand these changes and adjust your input file \n"
                     << "* accordingly. Then add\n"
                     << "* OPTION, VERSION = " << currentVersion << ";\n"
                     << "* to your input file. " << endl;
            errorMsg << "\n************************************************************************"
                        "****\n"
                     << endl;
            throw OpalException("TrackRun::execute", "Version mismatch");
        }
    }

    isFollowupTrack_m = opal->hasBunchAllocated();
    if (!itsAttr[TRACKRUN::DISTRIBUTION] && !isFollowupTrack_m) {
        throw OpalException(
            "TrackRun::execute", "\"DISTRIBUTION\" must be set in \"RUN\" command.");
    }
    if (!itsAttr[TRACKRUN::FIELDSOLVER]) {
        throw OpalException("TrackRun::execute", "\"FIELDSOLVER\" must be set in \"RUN\" command.");
    }
    if (!itsAttr[TRACKRUN::BEAM]) {
        throw OpalException("TrackRun::execute", "\"BEAM\" must be set in \"RUN\" command.");
    }

    // Get algorithm to use.
    setRunMethod();
    switch (method_m) {
        case RunMethod::PARALLEL: {
            setupTracker();
            break;
        }
        default: {
            throw OpalException("TrackRun::execute", "Unknown \"METHOD\" for the \"RUN\" command");
        }
    }

    /// \todo itsTracker_m->execute();

    opal->setRestartRun(false);
    opal->bunchIsAllocated();

    /// \todo do we delete here itsTracker_m;
}

void TrackRun::setRunMethod() {
    if (!itsAttr[TRACKRUN::METHOD]) {
        throw OpalException(
            "TrackRun::setRunMethod", "The attribute \"METHOD\" isn't set for the \"RUN\" command");
    } else {
        auto it = stringMethod_s.right.find(Attributes::getString(itsAttr[TRACKRUN::METHOD]));
        if (it != stringMethod_s.right.end()) {
            method_m = it->second;
        }
    }
}

std::string TrackRun::getRunMethodName() const {
    return stringMethod_s.left.at(method_m);
}

void TrackRun::setupTracker() {
    OpalData::getInstance()->setInOPALTMode();

    if (isFollowupTrack_m) {
        Track::block->bunch->setLocalTrackStep(0);
    }

    Beam* beam = Beam::find(Attributes::getString(itsAttr[TRACKRUN::BEAM]));
    Track::block->bunch->setBeamFrequency(beam->getFrequency() * Units::MHz2Hz);
    Track::block->bunch->setPType(beam->getParticleName());

    setBoundaryGeometry();

    setupFieldsolver();

    if (opal->inRestartRun()) {
        phaseSpaceSink_m = new H5PartWrapperForPT(
            opal->getInputBasename() + std::string(".h5"), opal->getRestartStep(),
            OpalData::getInstance()->getRestartFileName(), H5_O_WRONLY);
    } else if (isFollowupTrack_m) {
        phaseSpaceSink_m = new H5PartWrapperForPT(
            opal->getInputBasename() + std::string(".h5"), -1,
            opal->getInputBasename() + std::string(".h5"), H5_O_WRONLY);
    } else {
        phaseSpaceSink_m =
            new H5PartWrapperForPT(opal->getInputBasename() + std::string(".h5"), H5_O_WRONLY);
    }

    macrocharge_m = setDistributionParallelT(beam);
    macromass_m   = beam->getMassPerParticle();

    *gmsg << *this << endl;

    Track::block->bunch->setdT(Track::block->dT.front());
    // ada Track::block->bunch->dtScInit_m = Track::block->dtScInit;
    // ada Track::block->bunch->deltaTau_m = Track::block->deltaTau;

    if (!isFollowupTrack_m && !opal->inRestartRun()) {
        Track::block->bunch->setT(Track::block->t0_m);
    }

    Track::block->bunch->setCharge(macrocharge_m);
    Track::block->bunch->setMass(macromass_m);

    // set coupling constant
    double coefE = 1.0 / (4 * Physics::pi * Physics::epsilon_0);
    Track::block->bunch->setCouplingConstant(coefE);

    // statistical data are calculated (rms, eps etc.)
    Track::block->bunch->calcBeamParameters();

    initDataSink();

    if (!isFollowupTrack_m) {
        *gmsg << std::scientific;
        *gmsg << *dist_m << endl;
    }

    if (Track::block->bunch->getTotalNum() > 0) {
        double spos = Track::block->zstart;
        auto& zstop = Track::block->zstop;
        auto it     = Track::block->dT.begin();

        unsigned int i = 0;
        while (i + 1 < zstop.size() && zstop[i + 1] < spos) {
            ++i;
            ++it;
        }

        Track::block->bunch->setdT(*it);
    } else {
        Track::block->zstart = 0.0;
    }

    //*gmsg << *beam << endl;
    //*gmsg << *fs_m << endl;

    /*
      This needs to come back as soon as we have RF

      findPhasesForMaxEnergy();

    */

    /*
    itsTracker_m = new ParallelTracker(
        *Track::block->use->fetchLine(), Track::block->bunch, *ds_m, Track::block->reference, false,
        Attributes::getBool(itsAttr[TRACKRUN::TRACKBACK]), Track::block->localTimeSteps,
        Track::block->zstart, Track::block->zstop, Track::block->dT);
    */
}

void TrackRun::setupFieldsolver() {
    /*
    fs_m = FieldSolver::find(Attributes::getString(itsAttr[TRACKRUN::FIELDSOLVER]));

    if (fs_m->getFieldSolverType() != FieldSolverType::NONE) {
        size_t numGridPoints =
            fs_m->getMX() * fs_m->getMY() * fs_m->getMZ();  // total number of gridpoints
        Beam* beam          = Beam::find(Attributes::getString(itsAttr[TRACKRUN::BEAM]));
        size_t numParticles = beam->getNumberOfParticles();

        if (!opal->inRestartRun() && numParticles < numGridPoints) {
            throw OpalException(
                "TrackRun::setupFieldsolver()",
                "The number of simulation particles (" + std::to_string(numParticles) + ") \n"
                    + "is smaller than the number of gridpoints (" + std::to_string(numGridPoints)
                    + ").\n"
                    + "Please increase the number of particles or reduce the size of the mesh.\n");
        }

        OpalData::getInstance()->addProblemCharacteristicValue("MX", fs_m->getMX());
        OpalData::getInstance()->addProblemCharacteristicValue("MY", fs_m->getMY());
        OpalData::getInstance()->addProblemCharacteristicValue("MT", fs_m->getMZ());
    }
    */
    // fs_m->initCartesianFields();
    // Track::block->bunch->setSolver(fs_m);
    /*
    if (fs_m->hasPeriodicZ()) {
        Track::block->bunch->setBCForDCBeam();
    } else {
        Track::block->bunch->setBCAllOpen();
    }
    */
}

void TrackRun::initDataSink() {
    if (!opal->inRestartRun()) {
        if (!opal->hasDataSinkAllocated()) {
            opal->setDataSink(new DataSink(phaseSpaceSink_m, false));
        } else {
            ds_m = opal->getDataSink();
            ds_m->changeH5Wrapper(phaseSpaceSink_m);
        }
    } else {
        opal->setDataSink(new DataSink(phaseSpaceSink_m, true));
    }
    ds_m = opal->getDataSink();
}

void TrackRun::setBoundaryGeometry() {
    if (Attributes::getString(itsAttr[TRACKRUN::BOUNDARYGEOMETRY]) != "NONE") {
        // Ask the dictionary if BoundaryGeometry is allocated.
        // If it is allocated use the allocated BoundaryGeometry
        if (!OpalData::getInstance()->hasGlobalGeometry()) {
            const std::string geomDescriptor =
                Attributes::getString(itsAttr[TRACKRUN::BOUNDARYGEOMETRY]);
            BoundaryGeometry* bg = BoundaryGeometry::find(geomDescriptor)->clone(geomDescriptor);
            OpalData::getInstance()->setGlobalGeometry(bg);
        }
    }
}

double TrackRun::setDistributionParallelT(Beam* beam) {
    /*
     * Distribution(s) can be set via a single distribution or a list
     * (array) of distributions. If an array is defined the first in the
     * list is treated as the primary distribution. All others are added to
     * it to create the full distribution.
     */
    std::vector<std::string> distributionArray =
        Attributes::getStringArray(itsAttr[TRACKRUN::DISTRIBUTION]);
    const size_t numberOfDistributions = distributionArray.size();

    dist_m = Distribution::find(defaultDistribution_m);

    /*
     * Initialize distributions.
     */
    size_t numberOfParticles = beam->getNumberOfParticles();

    // Track::block->bunch->setDistribution(dist_m, distrs_m, numberOfParticles);
    // Return charge per macroparticle.
    return beam->getChargePerParticle();
}

Inform& TrackRun::print(Inform& os) const {
    os << endl;
    os << "* ************* T R A C K  R U N *************************************************** "
       << endl;
    if (!isFollowupTrack_m) {
        os << "* Selected Tracking Method == " << getRunMethodName() << ", NEW TRACK" << '\n'
           << "* "
              "********************************************************************************** "
           << '\n';
    } else {
        os << "* Selected Tracking Method == " << getRunMethodName() << ", FOLLOWUP TRACK" << '\n'
           << "* "
              "********************************************************************************** "
           << '\n';
    }
    os << "* Phase space dump frequency    = " << Options::psDumpFreq << '\n'
       << "* Statistics dump frequency     = " << Options::statDumpFreq << " w.r.t. the time step."
       << '\n'
       << "* DT                            = " << Track::block->dT.front() << " [s]\n"
       << "* MAXSTEPS                      = " << Track::block->localTimeSteps.front() << '\n'
       << "* Mass of simulation particle   = " << macromass_m << " [GeV/c^2]" << '\n'
       << "* Charge of simulation particle = " << macrocharge_m << " [C]" << '\n';
    os << "* ********************************************************************************** ";
    return os;
}
