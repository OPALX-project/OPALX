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

#include "Algorithms/ParallelTracker.h"

#include "AbstractObjects/BeamSequence.h"

#include "AbstractObjects/OpalData.h"

#include "Attributes/Attributes.h"

#include "Beamlines/TBeamline.h"

#include "BasicActions/Option.h"

#include "Distribution/Distribution.h"

#include "Distribution/Gaussian.h"

#include "Distribution/MultiVariateGaussian.h"

#include "Distribution/FlatTop.h"

#include "Distribution/FromFile.h"

#include "Physics/Physics.h"
#include "Physics/Units.h"

#include "Track/Track.h"

#include "Utilities/OpalException.h"

#include "Lines/EmissionSourceList.h"
#include "Structure/Beam.h"
#include "Structure/BoundaryGeometry.h"
#include "Structure/DataSink.h"
#include "Structure/EmissionSource.h"
#include "Structure/H5PartWrapper.h"
#include "Structure/H5PartWrapperForPT.h"

#include "BuildInfo.h"
#include "Utility/Inform.h"
#include "changes.h"

#include "Utilities/BiMap.h"

#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <memory>
#include <vector>

#include <unistd.h>

#include <filesystem>

extern Inform* gmsg;

namespace {

/** Restart source path for container @p index when using per-container H5 names. */
std::string h5RestartSourceForContainer(
    const std::string& restartFile, const std::string& containerH5FileName, size_t numContainers) {
    if (numContainers <= 1) {
        return restartFile;
    }
    namespace fs = std::filesystem;
    fs::path rf(restartFile);
    fs::path leaf = fs::path(containerH5FileName).filename();
    if (rf.has_parent_path()) {
        return (rf.parent_path() / leaf).string();
    }
    return leaf.string();
}

}  // namespace

namespace TRACKRUN {
    // The attributes of class TrackRun.
    enum {
        METHOD,            // Tracking method to use.
        TURNS,             // The number of turns to be tracked, we keep that for the moment
        FIELDSOLVER,       // The field solver attached
        BOUNDARYGEOMETRY,  // The boundary geometry
        TRACKBACK,         // In case we run the beam backwards
        SIZE
    };
}  // namespace TRACKRUN

const BiMap<TrackRun::RunMethod, std::string> TrackRun::stringMethod_s = []() {
    BiMap<TrackRun::RunMethod, std::string> bimap;
    bimap.insert(TrackRun::RunMethod::PARALLEL, "PARALLEL");
    return bimap;
}();

TrackRun::TrackRun()
    : Action(
        TRACKRUN::SIZE, "RUN",
        "The \"RUN\" sub-command tracks the defined particles through "
        "the given lattice."),
      itsTracker_m(nullptr),
      fs_m(nullptr),
      ds_m(nullptr),
      phaseSpaceSinks_m(),
      isFollowupTrack_m(false),
      method_m(RunMethod::NONE){

    itsAttr[TRACKRUN::METHOD] = Attributes::makePredefinedString(
        "METHOD", "Name of tracking algorithm to use.", {"PARALLEL"});

    itsAttr[TRACKRUN::TURNS] = Attributes::makeReal(
        "TURNS",
        "Number of turns to be tracked; Number of neighboring bunches to be tracked in cyclotron.",
        1.0);

    itsAttr[TRACKRUN::FIELDSOLVER] =
        Attributes::makeString("FIELDSOLVER", "Field solver to be used.");

    itsAttr[TRACKRUN::BOUNDARYGEOMETRY] = Attributes::makeString(
        "BOUNDARYGEOMETRY", "Boundary geometry to be used NONE (default).", "NONE");

    itsAttr[TRACKRUN::TRACKBACK] =
        Attributes::makeBool("TRACKBACK", "Track in reverse direction, default: false.", false);

    registerOwnership(AttributeHandler::SUB_COMMAND);
    opal_m = OpalData::getInstance();
}

TrackRun::TrackRun(const std::string& name, TrackRun* parent)
    : Action(name, parent),
      itsTracker_m(nullptr),
      fs_m(nullptr),
      ds_m(nullptr),
      phaseSpaceSinks_m(),
      isFollowupTrack_m(false),
      method_m(RunMethod::NONE){
    /*
      the opal dictionary
    */

    opal_m = OpalData::getInstance();

    const Vector_t<int, 3> nr(8);

    ippl::NDIndex<3> domain;
    for (unsigned i = 0; i < 3; i++) {
        domain[i] = ippl::Index(nr[i]);
    }

    std::array<bool, 3> isParallel;

    for (unsigned d = 0; d < 3; ++d) {
        isParallel[d] = true;
    }
}

TrackRun::~TrackRun() {
}

TrackRun* TrackRun::clone(const std::string& name) {
    return new TrackRun(name, this);
}

void TrackRun::execute() {
   
    const int currentVersion = ((buildinfo::version_major * 100) + buildinfo::version_minor) * 100;

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
   
    // Follow-up behavior is still based on whether a bunch was allocated already.
    // Emission sources are resolved from the selected BEAM later.
    isFollowupTrack_m = opal_m->hasBunchAllocated();
    if (!itsAttr[TRACKRUN::FIELDSOLVER]) {
        throw OpalException("TrackRun::execute", "\"FIELDSOLVER\" must be set in \"RUN\" command.");
    }

    // Field solver commands are registry-owned by OpalData; TrackRun only borrows it.
    fs_m = FieldSolverCmd::find(Attributes::getString(itsAttr[TRACKRUN::FIELDSOLVER]));
    *gmsg << level1 << *fs_m << endl;
    if (fs_m->hasBinningCmd()) {
        *gmsg << level1 << *fs_m->getBinningCmd() << endl;
    }

    // Process BEAM object names
    std::vector<std::string> beamNames = Track::block->beamNames_m;
    if (beamNames.empty()) {
        throw OpalException("TrackRun::execute", 
            "No beam specified: set TRACK::BEAM or TRACK::BEAMS.");
    }

    // Create vector of BEAMs
    std::vector<Beam*> beams;
    beams.reserve(beamNames.size());
    for (const auto& name : beamNames) {
        if (name.empty()) {
            throw OpalException("TrackRun::execute", "Empty beam name in resolved beam list.");
        }
        beams.push_back(Beam::find(name));  // fail fast
    }
    for (const auto* b : beams) {
        if (b->isPhoton()) {
            throw OpalException(
                "TrackRun::execute",
                "TRACK does not support BEAM, PARTICLE=PHOTON yet. "
                "Photon beams may be defined for future OPALX features, but they are currently "
                "rejected during tracking.");
        }
    }
    *gmsg << level1 << "* RUN resolved beams: ";
    for (size_t i = 0; i < beamNames.size(); ++i) {
        *gmsg << beamNames[i] << (i + 1 < beamNames.size() ? ", " : "");
    }
    *gmsg << endl;
    // Print the BEAM banner for each resolved beam.
    for (Beam* b : beams) {
        *gmsg << level1 << *b << endl;
    }

    // Vectors for each species
    std::vector<double> macrocharges;
    std::vector<double> macromasses;
    std::vector<std::vector<EmissionSource*>> emissionSourcesLists;
    macrocharges.reserve(beams.size());
    macromasses.reserve(beams.size());
    emissionSourcesLists.reserve(beams.size());

    // Fill macro quantities and emissionSourceList per container (beam)
    for (size_t i = 0; i < beams.size(); ++i) {
        Beam* b = beams[i];

        const double macrocharge = b->getChargePerParticle();
        const double macromass   = b->getMassPerParticle();
        macrocharges.push_back(macrocharge);
        macromasses.push_back(macromass);

        const double part_per_macro_ratio = macrocharge / (b->getCharge() * Physics::q_e);
        *gmsg << level2 << "* Beam[" << i << "] " << beamNames[i]
              << " macro charge per particle [C]: " << macrocharge << endl;
        *gmsg << level2 << "* Beam[" << i << "] " << beamNames[i]
              << " macro mass per particle [GeV/c^2]: " << macromass << endl;
        *gmsg << level2 << "* Beam[" << i << "] " << beamNames[i]
              << " particles per macro particle: " << part_per_macro_ratio 
              << endl << endl;

        EmissionSourceList* esl = EmissionSourceList::find(b->getEmissionSourceListName());
        const auto& sources = esl->fetchSources();
        if (sources.empty()) {
            throw OpalException("TrackRun::execute",
                                "Emission sources list for beam '" + beamNames[i] +
                                "' must contain at least one EMISSIONSOURCE.");
        }
        emissionSourcesLists.emplace_back(sources.begin(), sources.end());
    }

    /*
    Need the following units for mass and charge:
    - Charge per macro particle in [C], this should be macrocharge_m or q_m in the bunch. 
      This will be used for the field calculations.
    - The pusher needs consistent units: eV for mass and elementary charges for charge. 
      This will (hopefully) be handled inside the pusher routines!
    */

    // ? Need to see how this interacts with multiple containers
    initDataSink(beams.size());

    // Set total particles per container (beam)
    std::vector<size_t> totalParticlesPerBeam(beams.size());
    for (size_t i = 0; i < beams.size(); ++i) {
        Beam* b = beams[i];
        totalParticlesPerBeam[i] = computeTotalAllocationForBunch(b, emissionSourcesLists[i]);
    }

    // Create PartBunch (PIC Manager) with multiple particle containers
    bunch_m = std::make_unique<bunch_type>(
        macrocharges,           // Macro charge [C]
        macromasses,            // Macro Mass [GeV]
        beams,                  // Beam objects per container
        totalParticlesPerBeam,  // Per-beam particle counts for allocation
        1.0,                    // lbt
        "LF2",                  // Integrator
        fs_m,                   // Fieldsolver
        ds_m);                  // Data sink

    // Validate container setup produced by constructor
    const auto& particleContainers = bunch_m->getParticleContainers();
    if (particleContainers.size() != beams.size()) {
        throw OpalException("TrackRun::execute",
                            "Mismatch between number of beams and particle containers.");
    }

    // BC handler
    *gmsg << level2 << *(bunch_m->getBCHandler()) << endl;

    setupBoundaryGeometry();

    // Get algorithm to use.
    setRunMethod();

    switch (method_m) {
        case RunMethod::PARALLEL: {
            break;
        }
        default: {
            throw OpalException("TrackRun::execute", "Unknown \"METHOD\" for the \"RUN\" command");
        }
    }


    //double deltaP = Attributes::getReal(itsAttr[Distribution::OFFSETP]);
    //if (inputMoUnits_m == InputMomentumUnits::EVOVERC) {
    //    deltaP = Util::convertMomentumEVoverCToBetaGamma(deltaP, beam->getM());
    //}

    if (ippl::Comm->rank() == 0) {
        long number_of_processors = sysconf(_SC_NPROCESSORS_ONLN);
        *gmsg << level5 << "sysconf(_SC_NPROCESSORS_ONLN)= " << number_of_processors << endl;

        // *gmsg << "omp_get_max_threads() " << omp_get_max_threads() << endl;

        int world_size;
        MPI_Comm_size( MPI_COMM_WORLD, &world_size );
        *gmsg << level5 << "MPI_Comm_size= " << world_size << endl;
    }

    // Setup all distributions and samplers, perform initial sampling (t0 == 0),
    // and prepare per-container emitting sampler lists for ParallelTracker.
    // Do this for each particle container
    std::vector<emittingSamplers_t> emittingSamplersList(particleContainers.size());
    for(size_t i=0; i<particleContainers.size(); ++i){
        setupDistributionsAndSamplers(
            emissionSourcesLists[i], 
            beams[i],
            emittingSamplersList[i],
            i);
    }

    // Calculate extents and update moments for each container
    bunch_m->bunchUpdate();
    bunch_m->print(*gmsg);

    // Set ZStart, ZStop, and dT
    if (bunch_m->getParticleContainer()->getTotalNum() > 0) {
        double spos = Track::block->zstart;
        auto& zstop = Track::block->zstop;
        auto it     = Track::block->dT.begin();

        unsigned int i = 0;
        while (i + 1 < zstop.size() && zstop[i + 1] < spos) {
            ++i;
            ++it;
        }

        bunch_m->setdT(*it);
    } else {
        Track::block->zstart = 0.0;
    }

    /* \todo this is also not unsed in the master.
       This needs to come back as soon as we have RF

       findPhasesForMaxEnergy();

    */
    itsTracker_m = std::make_unique<ParallelTracker>(
        *Track::block->use->fetchLine(), *bunch_m, ds_m, false,
        Track::block->localTimeSteps,
        Track::block->zstart, Track::block->zstop, Track::block->dT, emittingSamplersList);
    itsTracker_m->execute();

    /*
    opal_m->setRestartRun(false);

    opal_m->bunchIsAllocated();
    */

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


void TrackRun::initDataSink(size_t numParticleContainers) {
    phaseSpaceSinks_m.clear();
    phaseSpaceSinks_m.reserve(numParticleContainers);

    const std::string base = opal_m->getInputBasename();

    for (size_t i = 0; i < numParticleContainers; ++i) {
        const std::string stem = DataSink::diagnosticStemForContainer(base, numParticleContainers, i);
        const std::string dest   = stem + std::string(".h5");

        if (opal_m->inRestartRun()) {
            const std::string src = h5RestartSourceForContainer(
                OpalData::getInstance()->getRestartFileName(), dest, numParticleContainers);
            phaseSpaceSinks_m.push_back(
                std::make_unique<H5PartWrapperForPT>(dest, opal_m->getRestartStep(), src, H5_O_WRONLY));
        } else if (isFollowupTrack_m) {
            phaseSpaceSinks_m.push_back(
                std::make_unique<H5PartWrapperForPT>(dest, -1, stem + std::string(".h5"), H5_O_WRONLY));
        } else {
            phaseSpaceSinks_m.push_back(std::make_unique<H5PartWrapperForPT>(dest, H5_O_WRONLY));
        }
    }

    const std::vector<H5PartWrapper*> sinks = borrowedPhaseSpaceSinks();
    if (!opal_m->inRestartRun()) {
        if (!opal_m->hasDataSinkAllocated()) {
            opal_m->setDataSink(new DataSink(sinks, false, numParticleContainers));
        } else {
            DataSink* raw = opal_m->getDataSink();
            raw->changeH5Wrappers(sinks);
        }
    } else {
        opal_m->setDataSink(new DataSink(sinks, true, numParticleContainers));
    }

    // DataSink lifetime is managed by OpalData; TrackRun only borrows it.
    ds_m = opal_m->getDataSink();
}

std::vector<H5PartWrapper*> TrackRun::borrowedPhaseSpaceSinks() const {
    std::vector<H5PartWrapper*> sinks;
    sinks.reserve(phaseSpaceSinks_m.size());
    for (const auto& sink : phaseSpaceSinks_m) {
        sinks.push_back(sink.get());
    }
    return sinks;
}

void TrackRun::setupBoundaryGeometry() {
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

size_t TrackRun::computeTotalAllocationForBunch(
    Beam* beam,
    const std::vector<EmissionSource*>& sources) const {
    size_t beamAllocSize = beam->getNumAlloc();

    size_t totalFromDists = 0;
    for (EmissionSource* src : sources) {
        auto* dist = Distribution::find(src->getDistributionName());
        totalFromDists += dist->getNumParticles();
    }

    if (totalFromDists > 0) {
        *gmsg << level3
              << "* Sum of per-distribution NPARTDIST over all emission sources = "
              << totalFromDists << ", BEAM::NALLOC = " << beamAllocSize << endl;
        if (totalFromDists > beamAllocSize) {
            *gmsg << level1
                  << "* WARNING: Sum of NPARTDIST over all distributions ("
                  << totalFromDists
                  << ") exceeds BEAM::NALLOC (" << beamAllocSize
                  << "). Allocation baseline may be insufficient; "
                  << "macro-charge per particle is still derived from BEAM::NALLOC." << endl;
        }
        return totalFromDists;
    }

    return beamAllocSize;
}

void TrackRun::setupDistributionsAndSamplers(
    const std::vector<EmissionSource*>& sources,
    Beam* beam,
    emittingSamplers_t& emittingSamplers,
    size_t index) {
    static IpplTimings::TimerRef samplingTime = IpplTimings::getTimer("samplingTime");

    IpplTimings::startTimer(samplingTime);

    // Common containers / parameters used by all samplers.
    auto pc               = bunch_m->getParticleContainer(index);
    auto fc               = bunch_m->getFieldContainer();
    Vector_t<int, Dim> nr = bunch_m->nr_m;
    const double avrgpz   = beam->getMomentum() / beam->getMass();

    emittingSamplers.clear();
    distrs_m.clear();

    for (EmissionSource* src : sources) {
        // Distribution objects are registry-owned; samplers only borrow them.
        Distribution* opalDist = Distribution::find(src->getDistributionName());

        // Ensure distribution parameters and reference momentum are up to date.
        opalDist->setDistType();
        opalDist->setDist();
        opalDist->setAvrgPz(avrgpz);

        // FROMFILE distributions carry absolute momenta — the BEAM's
        // PC/ENERGY/GAMMA would be silently ignored, so forbid the combination.
        if (opalDist->getType() == DistributionType::FROMFILE) {
            if (beam->hasExplicitEnergy()) {
                throw OpalException(
                    "TrackRun::setupDistributionsAndSamplers()",
                    "FROMFILE distribution \"" + src->getDistributionName()
                    + "\" cannot be combined with PC/ENERGY/GAMMA on the BEAM. "
                      "Remove the energy attribute from the BEAM command — "
                      "particle momenta are read from the file.");
            }
        } else {
            if (!beam->hasExplicitEnergy()) {
                throw OpalException(
                    "TrackRun::setupDistributionsAndSamplers()",
                    "The energy hasn't been set. "
                    "Set either \"GAMMA\", \"ENERGY\" or \"PC\" on the BEAM command.");
            }
        }

        distrs_m.push_back(opalDist);

        // Build a sampler instance for this emission source.
        std::shared_ptr<SamplingBase> sampler;
        switch (opalDist->getType()) {
            case DistributionType::GAUSS:
                sampler = std::make_shared<Gaussian>(pc, fc, opalDist);
                break;
            case DistributionType::MULTIVARIATEGAUSS:
                sampler = std::make_shared<MultiVariateGaussian>(pc, fc, opalDist);
                break;
            case DistributionType::FLATTOP:
                sampler = std::make_shared<FlatTop>(pc, fc, opalDist);
                break;
            case DistributionType::FROMFILE:
                sampler = std::make_shared<FromFile>(pc, fc, opalDist);
                break;
            default:
                throw OpalException("Distribution::create",
                                    "Unknown \"TYPE\" of \"DISTRIBUTION\"");
        }

        sampler->setBunchStateHandler(bunch_m->getBunchStateHandler());

        // Per-source emission offsets, start time, and emission model.
        const auto  R0  = src->getR0();
        const auto  P0  = src->getP0();
        const double t0 = src->getT0();
        sampler->setEmissionOffsets(R0, P0, t0, src->getEmissionModel());

        const size_t Ndist = opalDist->getNumParticles();
        size_t       Nmutable = Ndist;
    
        // Always call generateParticles once per source; time-independent samplers
        // will internally early-return when t0 > 0, while FlatTop will set up its
        // emission structures irrespective of t0.
        sampler->generateParticles(Nmutable, nr);

        // Time-dependent (emitted) distributions (e.g. FlatTop) and delayed
        // one-shot injectors (t0 > 0) participate in emitParticles(t, dt)
        // during tracking.
        if (opalDist->emitting_m || src->getT0() > 0.0) {
            emittingSamplers.push_back(sampler);
            *gmsg << level2 << "* Configured emitting source of type "
                  << opalDist->getTypeofDistribution() << " with NPARTDIST = "
                  << Ndist << ", t0 = " << t0 << endl;
        }
    }

    *gmsg << level2 << "* Particle sampling / sampler setup for all emission sources done."
          << endl;
    IpplTimings::stopTimer(samplingTime);
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
       << "* MAXSTEPS                      = " << Track::block->localTimeSteps.front() << '\n';

    std::string primaryBeamName;
    if (Track::block && !Track::block->beamNames_m.empty()) {
        primaryBeamName = Track::block->beamNames_m.front();
    }

    if (!primaryBeamName.empty()) {
        Beam* beam = Beam::find(primaryBeamName);
        os << "* Mass of simulation particle   = " << beam->getMassPerParticle() << " [GeV/c^2]" << '\n'
           << "* Charge of simulation particle = " << beam->getChargePerParticle() << " [C]" << '\n';
    } else {
        os << "* Mass of simulation particle   = <unresolved>" << '\n'
           << "* Charge of simulation particle = <unresolved>" << '\n';
    }
    os << "* ********************************************************************************** ";
    return os;
}
