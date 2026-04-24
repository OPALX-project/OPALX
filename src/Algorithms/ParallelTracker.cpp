/**
 * @file ParallelTracker.cpp
 * @brief Implementation of ParallelTracker (OPAL-T time-based parallel tracker).
 *
 * @copyright Copyright (c) 200x - 2014, Christof Kraus, Paul Scherrer Institut, Villigen PSI,
 * Switzerland
 * @copyright 2015 - 2016, Christof Metzger-Kraus, Helmholtz-Zentrum Berlin, Germany
 * @copyright 2017 - 2020, Christof Metzger-Kraus
 *
 * All rights reserved. This file is part of OPAL.
 *
 * OPAL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * You should have received a copy of the GNU General Public License
 * along with OPAL. If not, see https://www.gnu.org/licenses/.
 */
#include "Algorithms/ParallelTracker.h"

#include <array>
#include <cfloat>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "Algorithms/Matrix.h"
#include "BasicActions/DumpEMFields.h"

#include "AbstractObjects/OpalData.h"
#include "Algorithms/CavityAutophaser.h"
#include "Algorithms/OrbitThreader.h"
#include "BasicActions/Option.h"

#include "Beamlines/Beamline.h"
#include "Beamlines/FlaggedBeamline.h"
#include "Distribution/Distribution.h"
#include "PartBunch/BinnedFieldSolver.h"
#include "Physics/Units.h"

#include "Structure/BoundaryGeometry.h"
#include "Structure/BoundingBox.h"
#include "Utilities/LogicalError.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utilities/Timer.h"
#include "Utilities/Util.h"
#include "ValueDefinitions/RealVariable.h"

#include "AbsBeamline/PluginElement.h"
#include "AbsBeamline/VerticalFFAMagnet.h"

extern Inform* gmsg;

// --- Constructors ---

/**
 * @brief Construct tracker from a beamline only (see class constructor overload).
 */
ParallelTracker::ParallelTracker(const Beamline& beamline, bool revBeam)
    : Tracker(beamline, revBeam, false),
      itsDataSink_m(),
      itsOpalBeamline_m(beamline.getOrigin3D(), beamline.getInitialDirection()),
      globalEOL_m(false),
      zstart_m(0.0),
      dtCurrentTrack_m(0.0),
      repartFreq_m(0),
      timeIntegrationTimer1_m(IpplTimings::getTimer("TIntegration1")),
      timeIntegrationTimer2_m(IpplTimings::getTimer("TIntegration2")),
      fieldEvaluationTimer_m(IpplTimings::getTimer("External field eval")),
      PluginElemTimer_m(IpplTimings::getTimer("PluginElements")),
      BinRepartTimer_m(IpplTimings::getTimer("Binaryrepart")),
      OrbThreader_m(IpplTimings::getTimer("OrbThreader")) {}

/**
 * @brief Construct tracker with bunch, data sink, z-segments, and optional emitters.
 */
ParallelTracker::ParallelTracker(
        const Beamline& beamline, PartBunch_t& bunch, DataSink* ds, bool revBeam,
        const std::vector<unsigned long long>& maxSteps, double zstart,
        const std::vector<double>& zstop, const std::vector<double>& dt,
        const std::vector<std::vector<std::shared_ptr<SamplingBase>>>& emittingSamplers)
    : Tracker(beamline, bunch, revBeam, false),
      itsDataSink_m(ds),
      itsOpalBeamline_m(beamline.getOrigin3D(), beamline.getInitialDirection()),
      globalEOL_m(false),
      zstart_m(zstart),
      dtCurrentTrack_m(0.0),
      repartFreq_m(0),
      emittingSamplers_m(emittingSamplers),
      timeIntegrationTimer1_m(IpplTimings::getTimer("TIntegration1")),
      timeIntegrationTimer2_m(IpplTimings::getTimer("TIntegration2")),
      fieldEvaluationTimer_m(IpplTimings::getTimer("External field eval")),
      BinRepartTimer_m(IpplTimings::getTimer("Binaryrepart")),
      OrbThreader_m(IpplTimings::getTimer("OrbThreader")) {
    for (unsigned int i = 0; i < zstop.size(); ++i) {
        stepSizes_m.push_back(dt[i], zstop[i], maxSteps[i]);
    }

    stepSizes_m.sortAscendingZStop();
    stepSizes_m.resetIterator();
}

/**
 * @copybrief ParallelTracker::~ParallelTracker
 */
ParallelTracker::~ParallelTracker() {}
// --- Visit functions ---

/**
 * @copybrief ParallelTracker::visitComponent
 */
void ParallelTracker::visitComponent(const Component& comp) {
    if (comp.getType() == ElementType::LASER) {
        throw LogicalError(
                "ParallelTracker::visitComponent()",
                "Tracking of the \"LASER\" element is not implemented yet.");
    }

    Tracker::visitComponent(comp);
}

/**
 * @brief Iterate beamline elements and dispatch into OpalBeamline.
 * @param bl Flagged beamline whose elements are visited.
 */
void ParallelTracker::visitBeamline(const Beamline& bl) {
    const FlaggedBeamline* fbl = static_cast<const FlaggedBeamline*>(&bl);
    /*
    if (fbl->getRelativeFlag()) {
        OpalBeamline stash(fbl->getOrigin3D(), fbl->getInitialDirection());
        stash.swap(itsOpalBeamline_m);
        fbl->iterate(*this, false);
        itsOpalBeamline_m.prepareSections();
        itsOpalBeamline_m.compute3DLattice();
        stash.merge(itsOpalBeamline_m);
        stash.swap(itsOpalBeamline_m);
    } else {
        fbl->iterate(*this, false);
    }
    */
    fbl->iterate(*this, false);
}

// --- execute() ---

/**
 * @copybrief ParallelTracker::execute
 */
void ParallelTracker::execute() {
    Inform m("ParallelTracker::execute");

    // PartBunch::resetPcActive() ran in the constructor while containers were still empty
    // (allocate-then-destroy for capacity). Initial particles are loaded later in TrackRun
    // without refreshing these flags, so reference updates must not skip all containers.
    itsBunch_m->resetPcActive();
    activateEmittingContainers(itsBunch_m->getT());

    // Initialize the Boris particle pusher
    BorisPusher pusher;
    m << level3 << "Initialized Boris pusher." << endl;

    // Ensure the time step is positive for the setup phase
    itsBunch_m->setdT(std::abs(itsBunch_m->getdT()));

    // Populate the OpalBeamline and calculate coordinate transformations
    prepareSections();

    // Select the minimal time step from the configuration
    double minTimeStep = stepSizes_m.getMinTimeStep();
    m << level3 << "Selected minimum time step from configuration: " << minTimeStep << endl;

    // Activate all beamline elements (sets Component::online_m = true)
    itsOpalBeamline_m.activateElements();
    m << level3 << "Activated all beamline elements." << endl;

    // Calculate the coordinate transformation from the beamline origin to the lab frame
    CoordinateSystemTrafo beamlineToLab = itsOpalBeamline_m.getCSTrafoLab2Local().inverted();

    // Per-container: lab transform and reference orbit state (each beam's PartData for P0
    // fallback).
    const auto& particleContainers = itsBunch_m->getParticleContainers();
    for (const auto& pc : particleContainers) {
        if (!pc) {
            continue;
        }
        if (!pc->getReference()) {
            throw OpalException(
                    "ParallelTracker::execute",
                    "Particle container has null PartData reference during lab-frame init.");
        }
        pc->setToLabTrafo(beamlineToLab);
        pc->getRefPartR() = beamlineToLab.transformTo(Vector_t<double, 3>(0, 0, 0));
        if (pc->getTotalNum() > 0) {
            pc->getRefPartP() = beamlineToLab.rotateTo(pc->getMeanP());
        } else {
            const PartData& pref = *pc->getReference();
            const double P0      = pref.getP() / pref.getM();  // beta*gamma from BEAM pc
            pc->getRefPartP()    = beamlineToLab.rotateTo(Vector_t<double, 3>(0.0, 0.0, P0));
        }
    }

    m << level4
      << "Transformed reference particle position and momentum to lab frame (all containers)."
      << endl;

    // Integrate reference orbits forward until all container path lengths reach zstart_m (when
    // needed).
    findStartPositions(pusher);

    // Advance through the time step configurations, skipping any that end
    // before the start position (zstart_m)
    stepSizes_m.advanceToPos(zstart_m);

    // Global spatial bounds: union over all containers
    Vector_t<double, 3> rmin(0.0), rmax(0.0);
    if (itsBunch_m->getTotalNumAllContainers() > 0) {
        computeInitialBounds(rmin, rmax);
        m << level4 << "Initial bunch bounds:\n"
          << "  rmin = " << rmin << "\n"
          << "  rmax = " << rmax << endl;
    }

    // Print reference particle information for all containers
    printInitialContainerRefs(m);

    // Start timing for the OrbitThreader section
    IpplTimings::startTimer(OrbThreader_m);

    // Flags to control phase space and statistics dumping for the initial step
    bool const psDump0   = 0;
    bool const statDump0 = 0;

    // Write initial phase space and statistics
    writePhaseSpace(0, psDump0, statDump0);
    m << level2 << "Dump initial phase space done." << endl;

    // Create an OrbitThreader object to handle orbit threading and element queries
    OrbitThreader oth(
            *itsBunch_m->getParticleContainer(0)->getReference(),
            itsBunch_m->getParticleContainer(0)->getRefPartR(),
            itsBunch_m->getParticleContainer(0)->getRefPartP(),
            itsBunch_m->getParticleContainer(0)->get_sPos(),
            -rmin(2),            // Negative minimum z bound
            itsBunch_m->getT(),  // Current bunch time
            minTimeStep,
            stepSizes_m,         // Step size configuration
            itsOpalBeamline_m);  // OpalBeamline object

    // Execute the orbit threading (queries elements and updates state)
    oth.execute();
    m << level4 << "Orbit threader execution done." << endl;

    // Stop timing for the OrbitThreader section
    IpplTimings::stopTimer(OrbThreader_m);

    // Get bounding box
    BoundingBox globalBoundingBox = oth.getBoundingBox();

    // Set the time view of the particle bunch
    setTime();
    m << level4 << "Set time view of particle bunch." << endl;

    // Reset the bunch time?
    double time = itsBunch_m->getT();
    itsBunch_m->setT(time);
    m << level4 << "Reset bunch time to " << time << "." << endl;

    // Get the current global tracking step
    unsigned long long step = itsBunch_m->getGlobalTrackStep();
    OPALTimer::Timer myt1;
    *gmsg << level1 << "* Track start at: " << myt1.time() << ", t= " << Util::getTimeString(time)
          << "; "
          << "zstart at: " << Util::getLengthString(itsBunch_m->getParticleContainer(0)->get_sPos())
          << endl
          << "* Initial dt = " << Util::getTimeString(itsBunch_m->getdT()) << endl
          << "* Max integration steps = " << stepSizes_m.getMaxSteps() << ", next step = " << step
          << endl;

    setOptionalVariables();

    globalEOL_m = false;
    // wakeStatus_m       = false;

    stepSizes_m.printDirect(*gmsg);

    // Before the tracker loop: bunch sanity checks.
    this->itsBunch_m->performBunchSanityChecks();

    // Handle any dump field requests
    DumpEMFields::writeFields(itsOpalBeamline_m.getElements());

    // Main tracking loop over step size configurations
    m << level5 << ">>>>>>>>>>>>>>>>>> Starting Tracking Loop >>>>>>>>>>>>>>>>>>" << endl;
    while (!stepSizes_m.reachedEnd()) {
        // Set the number of steps for the current track
        unsigned long long trackSteps = stepSizes_m.getNumSteps() + step;
        dtCurrentTrack_m              = stepSizes_m.getdT();

        // Select global dt from dtCurrentTrack_m and copy to all container dt views.
        changeDT();
        itsBunch_m->resetPcActive();
        activateEmittingContainers(itsBunch_m->getT());

        // Inner loop over the number of steps for the current configuration
        m << level2 << "Starting track with dt = " << Util::getTimeString(dtCurrentTrack_m)
          << ", track steps = " << step << " to " << trackSteps << "." << endl;
        for (; step < trackSteps; ++step) {
            if (!itsBunch_m->anyPcActive()) {
                m << level4 << "No active particle containers; ending inner track segment." << endl;
                break;
            }

            // Particle R and mesh are in REFERENCE frame for the whole step except inside
            // computeSpaceChargeFields (beam frame only during computeSelfFields).

            // Reset EOL flag each step: transient OutOfBounds (e.g. from invalid mesh
            // bounds immediately after first emission) must not persist across steps.
            globalEOL_m = false;

            // Get the bunch spatial bounds across all containers.
            Vector_t<double, 3> rmin(0.0), rmax(0.0);
            if (itsBunch_m->getTotalNumAllContainers() > 0) {
                computeInitialBounds(rmin, rmax);
            }

            // First half of the time integration
            timeIntegration1(pusher);
            m << level4 << "timeIntegration1 done at step " << step << "." << endl;
            itsBunch_m->bunchUpdate();
            m << level5 << "Bunch updated after timeIntegration1." << endl;

            // Reset E and B fields
            resetFields();
            m << level4 << "E and B fields reset at step " << step << "." << endl;

            // std::cout << "local num: " << itsBunch_m->getLocalNum() << std::endl;

            // Space charge field computation
            // if (itsBunch_m->getLocalNum() > 1) {
            // Otherwise no interaction, can skip (and for some reason seg-fault...)
            computeSpaceChargeFields(step);
            m << level4 << "Space charge field computation done at step " << step << "." << endl;
            //}

            // Emission is placed BETWEEN space-charge and external-field evaluation
            // to match the legacy OPAL ordering (ParallelTTracker): newly emitted
            // particles experience external fields in their creation step.
            selectDT();

            // Reset per-particle dt for all existing particles BEFORE emission, so that
            // newly emitted particles retain their fractional dt (sampled in generateUniformDisk).
            setTime();
            m << level5 << "Set time view of particle bunch to dt = "
              << Util::getTimeString(itsBunch_m->getdT()) << "." << endl;

            // Emit particles from time-dependent (emitting) sources (R set in REFERENCE frame).
            // New particles receive a fractional per-particle dt ∈ (0, dt) from the sampler.
            emitFromEmissionSources(itsBunch_m->getT(), itsBunch_m->getdT());
            m << level4 << "Emit particles from emission sources done at step " << step << "."
              << endl;
            itsBunch_m->bunchUpdate();
            m << level5 << "Bunch updated after emission." << endl;

            // selectDT(back_track);
            // m << level5 << "Selected new time step for next iteration, back_track = " <<
            // back_track << "." << endl;

            // External field computation
            computeExternalFields(oth);
            m << level4 << "External field computation done at step " << step << "." << endl;

            // Second half of the time integration
            timeIntegration2(pusher);
            m << level4 << "timeIntegration2 done at step " << step << "." << endl;
            itsBunch_m->bunchUpdate();
            m << level5 << "Bunch updated after timeIntegration2." << endl;

            if (itsBunch_m->getTotalNumAllContainers() == 0) {
                m << level5 << "WARNING: No particles in the bunch at step " << step << " on rank "
                  << ippl::Comm->rank() << ". This has no effect on the simulation." << endl;
            }

            // Update the bunch time
            itsBunch_m->incrementT();
            m << level5 << "Incremented bunch time to " << Util::getTimeString(itsBunch_m->getT())
              << "." << endl;

            // Update reference particle
            if (itsBunch_m->getT() > 0.0) {
                updateReference(pusher);
                m << level4 << "Updated reference particle at step " << step << "." << endl;
            }

            double sigmas = static_cast<double>(Options::boundpDestroy);
            // if (sigmas > 0.0) {
            size_t nDeleted = itsBunch_m->getParticleContainer()->deleteParticlesOutside(sigmas);
            if (nDeleted > 0) {
                m << level2 << "Deleted " << nDeleted << " particles outside " << sigmas
                  << "-sigma boundary, " << itsBunch_m->getTotalNumAllContainers() << " remaining."
                  << endl;
            }
            // }
            // deletedParticles_m = false;
            //}

            // Update the path length
            // itsBunch_m->set_sPos(pathLength_m);
            // m << level4 << "Updated path length to " << pathLength_m << "." << endl;

            // hasEndOfLineReached uses isOutside(RefPartR_m) with a lab-frame
            // bounding box, but RefPartR_m lives in the reference frame in OPALX
            // (near origin), so the check always fails.  Disabled until the
            // coordinate mismatch is resolved.
            // if (hasEndOfLineReached(globalBoundingBox)) break;

            // Dump phase space and statistics at configured intervals
            bool const psDump =
                    ((itsBunch_m->getGlobalTrackStep() % Options::psDumpFreq) + 1
                     == Options::psDumpFreq);
            bool const statDump =
                    ((itsBunch_m->getGlobalTrackStep() % Options::statDumpFreq) + 1
                     == Options::statDumpFreq);
            dumpStats(step, psDump, statDump);

            itsBunch_m->incTrackSteps();
            m << level5 << "Track steps incremented." << endl;

            for (size_t i = 0; i < particleContainers.size(); ++i) {
                const auto& pc = particleContainers[i];
                if (!pc || !itsBunch_m->isPcActive(i)) {
                    continue;
                }
                ippl::Vector<double, 3> pdivg =
                        pc->getRefPartP() / Util::getGamma(pc->getRefPartP());
                double beta             = euclidean_norm(pdivg);
                double driftPerTimeStep = std::abs(itsBunch_m->getdT()) * Physics::c * beta;
                m << level4 << "Calculated drift per time step (container " << i
                  << "): " << Util::getLengthString(driftPerTimeStep) << "." << endl;

                if (std::abs(stepSizes_m.getZStop() - pc->get_sPos()) < 0.5 * driftPerTimeStep) {
                    m << level2
                      << "Approaching end of current step size configuration for container " << i
                      << " (zstop = " << Util::getLengthString(stepSizes_m.getZStop())
                      << ", path length = " << Util::getLengthString(pc->get_sPos()) << ")."
                      << endl;
                    itsBunch_m->setPcAtZStop(i);
                }
            }
            if (!itsBunch_m->anyPcActive()) {
                m << level2
                  << "All active containers reached current zstop. Preparing to switch to "
                  << "next configuration." << endl;
                break;
            }
        }

        // globalEOL_m is reset at the start of each step, so if it is still true
        // here the last step genuinely ended with an out-of-bounds or reference
        // particle hitting an element. Synchronize across ranks before deciding.
        ippl::Comm->allreduce(globalEOL_m, 1, std::logical_and<bool>());

        if (globalEOL_m) break;
        ++stepSizes_m;
    }
    bool const psDump =
            (((itsBunch_m->getGlobalTrackStep() - 1) % Options::psDumpFreq) + 1
             != Options::psDumpFreq);
    bool const statDump =
            (((itsBunch_m->getGlobalTrackStep() - 1) % Options::statDumpFreq) + 1
             != Options::statDumpFreq);

    writePhaseSpace((step + 1), psDump, statDump);

    *gmsg << level2 << "* Dump phase space of last step" << endl;

    itsOpalBeamline_m.switchElementsOff();

    // Ensure all Kokkos operations are complete
    Kokkos::fence();

    if (itsBunch_m->hasFieldSolver()) {
        m << level2 << "Total FieldSolver calls: " << itsBunch_m->getFieldSolver()->getCallCounter()
          << endl;
    }

    OPALTimer::Timer myt3;
    *gmsg << level1 << endl
          << "* Done executing ParallelTracker at " << myt3.time() << endl
          << endl;
}

// --- PIC integration and fields ---

/**
 * @copybrief ParallelTracker::timeIntegration1
 */
void ParallelTracker::timeIntegration1(BorisPusher& pusher) {
    Inform m("ParallelTracker::timeIntegration1");
    IpplTimings::startTimer(timeIntegrationTimer1_m);
    const size_t n = itsBunch_m->getNumParticleContainers();
    for (size_t i = 0; i < n; ++i) {
        if (!itsBunch_m->isPcActive(i)) {
            continue;
        }
        auto pc = itsBunch_m->getParticleContainer(i);
        if (!pc) {
            continue;
        }
        pushParticles(pusher, *pc);
    }
    IpplTimings::stopTimer(timeIntegrationTimer1_m);
    m << level4 << "Push particles done for all containers." << endl;
}

/**
 * @copybrief ParallelTracker::timeIntegration2
 */
void ParallelTracker::timeIntegration2(BorisPusher& pusher) {
    // Legacy note: cathode transport/emission was sequenced after space charge so that
    // the first step of newborn particles omits self-fields; multi-container emission
    // is handled separately in execute().
    Inform m("ParallelTracker::timeIntegration2");

    IpplTimings::startTimer(timeIntegrationTimer2_m);
    const size_t n = itsBunch_m->getNumParticleContainers();
    for (size_t i = 0; i < n; ++i) {
        if (!itsBunch_m->isPcActive(i)) {
            continue;
        }
        auto pc = itsBunch_m->getParticleContainer(i);
        if (!pc) {
            continue;
        }
        kickParticles(pusher, *pc);
        pushParticles(pusher, *pc);
    }
    m << level4 << "Kick/push particles done for all containers." << endl;

    // Put the update directly into pushParticles instead!...
    // double newdT = itsBunch_m->getdT();
    // itsBunch_m->getParticleContainer()->dt = newdT;
    // m << "Update particle time step done." << endl;

    IpplTimings::stopTimer(timeIntegrationTimer2_m);
}

/**
 * @copybrief ParallelTracker::computeSpaceChargeFields
 *
 * @par Frame of reference
 * - Entry: @f$R@f$, @f$E@f$, @f$B@f$ in the reference (lab) frame.
 * - After transform to beam: @f$R@f$ in the beam frame (origin at reference, z along momentum).
 * - Inside computeSelfFields / bunchUpdate: mesh follows @f$R@f$, so mesh is in the beam frame.
 * - After transform back: @f$R@f$, @f$E@f$, @f$B@f$ in the reference frame again.
 * - After final bunchUpdate(): mesh matches reference-frame @f$R@f$.
 */
void ParallelTracker::computeSpaceChargeFields(unsigned long long step) {
    Inform m("ParallelTracker::computeSpaceChargeFields");
    // Current limitation: space-charge transform/scatter/gather is applied via the primary
    // container path only. Keep this behavior until the dedicated multi-container SC refactor.
    if (!itsBunch_m->hasFieldSolver()) {
        /*
        This should not happen, so when we do not have a field solve, we can
        throw an exception. If we have "no solver" and want to run it, we would
        choose the null solver.
        */
        *gmsg << level1 << "no solver available!" << endl;
        throw OpalException(
                "ParallelTracker::computeSpaceChargeFields",
                "Bunch has no field solver assigned! If you want to run without "
                "space charge effects, please use TYPE=NONE for the field solver.");
    }

    itsBunch_m->calcBeamParameters();
    m << level4 << "Calculate beam parameters done." << endl;

    // Use mean momentum for beam-frame alignment; with 0 or 1 particle get_pmean() can be
    // zero or negligible (e.g. rank with no particles), which would make getQuaternion throw.
    const double pmean_tol    = 1e-12;
    Vector_t<double, 3> pmean = itsBunch_m->getParticleContainer()->getMeanP();
    double pmean_len2         = dot(pmean, pmean);
    if (pmean_len2 < pmean_tol * pmean_tol) {
        pmean      = itsBunch_m->getParticleContainer()->getRefPartP();
        pmean_len2 = dot(pmean, pmean);
    }
    if (pmean_len2 < pmean_tol * pmean_tol) {
        pmean = Vector_t<double, 3>(0, 0, 1);
    }
    Quaternion alignment = getQuaternion(pmean, Vector_t<double, 3>(0, 0, 1));

    CoordinateSystemTrafo beamToReferenceCSTrafo(
            Vector_t<double, 3>(0, 0, itsBunch_m->getParticleContainer()->get_sPos()),
            alignment.conjugate());

    CoordinateSystemTrafo referenceToBeamCSTrafo = beamToReferenceCSTrafo.inverted();

    // Transform particle positions to the beam frame.
    referenceToBeamCSTrafo.transformBunchTo(
            itsBunch_m->getParticleContainer()->R.getView(),
            itsBunch_m->getParticleContainer()->getLocalNum());
    m << level4 << "Transform particle positions to beam coordinate system done." << endl;
    itsBunch_m->bunchUpdate();
    m << level5 << "Bunch updated for positions in beam coordinate system." << endl;

    // TODO: itsBunch_m->boundp() not implemented yet.
    // itsBunch_m->boundp();

    if (step % repartFreq_m + 1 == repartFreq_m) {
        doBinaryRepartition();
        m << level4 << "Binary repartition done." << endl;
    }

    // itsBunch_m->setGlobalMeanR(itsBunch_m->get_centroid());

    itsBunch_m->computeSelfFields();
    m << level3 << "Compute self fields done." << endl;

    // Transform positions back to the reference frame.
    const size_t nLocRef = itsBunch_m->getParticleContainer()->getLocalNum();
    beamToReferenceCSTrafo.transformBunchTo(
            itsBunch_m->getParticleContainer()->R.getView(), nLocRef);
    m << level5 << "Transform particle positions back to reference coordinate system done." << endl;
    // Rotate E and B back to the reference frame.
    beamToReferenceCSTrafo.rotateBunchTo(itsBunch_m->getParticleContainer()->E.getView(), nLocRef);
    m << level5 << "Rotate E fields back to reference coordinate system done." << endl;
    beamToReferenceCSTrafo.rotateBunchTo(itsBunch_m->getParticleContainer()->B.getView(), nLocRef);
    m << level5
      << "Rotate B fields back to reference coordinate system done. ComputeSelfFields done."
      << endl;
    // Rebuild mesh from reference-frame R (computeSelfFields used beam-frame bunchUpdate).
    itsBunch_m->bunchUpdate();
    m << level5 << "Bunch updated for positions in reference coordinate system." << endl;
}

/**
 * @copybrief ParallelTracker::computeExternalFields
 */
void ParallelTracker::computeExternalFields(OrbitThreader& oth) {
    IpplTimings::startTimer(fieldEvaluationTimer_m);
    Inform msg("ParallelTracker ", *gmsg);

    const size_t nContainers = itsBunch_m->getNumParticleContainers();
    for (size_t ci = 0; ci < nContainers; ++ci) {
        if (!itsBunch_m->isPcActive(ci)) {
            continue;
        }
        auto pc = itsBunch_m->getParticleContainer(ci);
        if (!pc) {
            continue;
        }

        // Bunch bounds for this container.
        Vector_t<double, 3> rmin(0.0), rmax(0.0);
        if (pc->getTotalNum() > 0) {
            pc->computeMinMaxR();
            rmin = pc->getMinR();
            rmax = pc->getMaxR();

            // get_bounds returns cached mesh extents (rmin_m/rmax_m). These are
            // sentinel values (DBL_MAX / DBL_MIN) when the mesh hasn't been
            // computed yet — e.g. immediately after the first particles are emitted
            // before calcBeamParameters() has run.  Fall back to zero so the query
            // uses pathLength_m as the centre with zero half-width.
            if (!std::isfinite(rmin(2)) || !std::isfinite(rmax(2))) {
                rmin = rmax = 0.0;
            }
        }

        // Get elements at bunch position.
        IndexMap::value_t elements;
        try {
            elements = oth.query(pc->get_sPos() + 0.5 * (rmax(2) + rmin(2)), rmax(2) - rmin(2));
        } catch (IndexMap::OutOfBounds& e) {
            globalEOL_m = true;
            continue;
        }

        // Iterate over all elements.
        IndexMap::value_t::const_iterator it        = elements.begin();
        const IndexMap::value_t::const_iterator end = elements.end();
        for (; it != end; ++it) {
            CoordinateSystemTrafo refToLocalCSTrafo =
                    (itsOpalBeamline_m.getMisalignment((*it))
                     * (itsOpalBeamline_m.getCSTrafoLab2Local((*it)) * pc->getToLabTrafo()));

            CoordinateSystemTrafo localToRefCSTrafo = refToLocalCSTrafo.inverted();

            (*it)->setCurrentSCoordinate(pc->get_sPos() + rmin(2));

            pc->transformBunch(refToLocalCSTrafo);

            // Apply element to this iteration's particle container.
            (*it)->apply(pc);

            pc->transformBunch(localToRefCSTrafo);
        }
    }

    IpplTimings::stopTimer(fieldEvaluationTimer_m);
}

/**
 * @copybrief ParallelTracker::emitFromEmissionSources
 */
void ParallelTracker::emitFromEmissionSources(double t, double dt) {
    const auto& containers = itsBunch_m->getParticleContainers();
    for (size_t ci = 0; ci < containers.size(); ++ci) {
        const auto& pc = containers[ci];
        if (!pc) {
            continue;
        }
        if (itsBunch_m->pcAtZStop(ci)) {
            continue;
        }

        // Record the extent of the position array and the current local particle count
        // before emission. If an internal resize (Kokkos::realloc) happens during
        // emission, the extent of R will change and we can flag this as an error.
        const size_t extentBeforeEmission = pc->R.getView().extent(0);

        CoordinateSystemTrafo refToGun =
                itsOpalBeamline_m.getCSTrafoLab2Local() * pc->getToLabTrafo();
        pc->transformBunch(refToGun);
        if (ci < emittingSamplers_m.size()) {
            for (const auto& sampler : emittingSamplers_m[ci]) {
                if (sampler) {
                    sampler->emitParticles(t, dt);
                }
            }
        }
        pc->transformBunch(refToGun.inverted());

        pc->setM(pc->getMassPerParticle());
        pc->setQ(pc->getChargePerParticle());
        // itsBunch_m->updateNumTotal(); // handled internally by ippl
        // itsBunch_m->bunchUpdate();

        // Sanity guard: the total number of macroparticles in the bunch must
        // never exceed the globally configured BEAM::NALLOC value. Overshooting
        // this limit would trigger internal reallocations in the particle
        // container and silently drop already-tracked particles/delete their data in the particle
        // attributes. This is only a check for the number of local particles.
        const size_t extentAfterEmission = pc->R.getView().extent(0);
        if (extentAfterEmission != extentBeforeEmission) {
            throw OpalException(
                    "ParallelTracker::emitFromEmissionSources",
                    "Local particle storage was resized during emission (likely due to "
                    "over-emission causing a Kokkos::realloc). This leads to loss of "
                    "previously tracked particles. Please increase the total number of "
                    "macroparticles or adjust NPARTDIST / the emission profile. If you are "
                    "using emission sources, please check the emission profile and adjust "
                    "the number of particles emitted.");
        }
    }
    itsBunch_m->refreshPcActiveAfterEmit();
}

/**
 * @copybrief ParallelTracker::resetFields
 */
void ParallelTracker::resetFields() {
    const size_t n = itsBunch_m->getNumParticleContainers();
    for (size_t i = 0; i < n; ++i) {
        if (!itsBunch_m->isPcActive(i)) {
            continue;
        }
        auto pc = itsBunch_m->getParticleContainer(i);
        if (!pc) {
            continue;
        }
        pc->E = 0;
        pc->B = 0;
    }
}

/**
 * @brief Boris position push in unitless coordinates (per-particle dt via pusher internals).
 * @param pusher Boris pusher.
 * @param pc     Non-null target particle container.
 */
void ParallelTracker::pushParticles(
        const BorisPusher& pusher, PartBunch_t::ParticleContainer_t& pc) {
    // Per-particle dt is used so that newly emitted particles with fractional dt
    // (sampled during emission) are pushed proportionally to their sub-timestep fraction.
    pc.switchToUnitlessPositions();

    auto Rview = pc.R.getView();
    auto Pview = pc.P.getView();
    // auto dtview = pc.dt.getView();

    Kokkos::parallel_for(
            "pushParticles", pc.getLocalNum(), KOKKOS_LAMBDA(const size_t i) {
                // Half drift: x_{n+1/2} = x_n + (dt/2) * v; unitless form via pusher.push(..., 0).
                // TODO: verify sign convention for half-step push.
                Vector_t<double, 3> x = Rview(i);
                pusher.push(
                        x, Pview(i),
                        0);  // this 0 is "dt" that is not used with unitless positions!
                Rview(i) = x;
            });

    pc.switchOffUnitlessPositions();
    // TODO: pc->update() changes results on single rank; keep disabled until investigated.
    // itsBunch_m->getParticleContainer()->update();
    Kokkos::fence();
    ippl::Comm->barrier();
    pc.markMomentsDirty();
    // itsBunch_m->bunchUpdate();
}

/**
 * @brief Boris velocity kick from E and B using per-particle dt.
 * @param pusher Boris pusher.
 * @param pc     Non-null target particle container.
 */
void ParallelTracker::kickParticles(
        const BorisPusher& pusher, PartBunch_t::ParticleContainer_t& pc) {
    Inform m("ParallelTracker::kickParticles");

    // auto Rview  = pc.R.getView();
    auto Pview  = pc.P.getView();
    auto dtview = pc.dt.getView();
    auto Efview = pc.E.getView();
    auto Bfview = pc.B.getView();
    m << level5 << "Got particle views for kick operation." << endl;
    // Mass (eV) and charge (proton charges) from this container's reference particle,
    // passed explicitly into BorisPusher::kick for GPU-safe kernels.

    const PartData& ref = *pc.getReference();
    const double mass   = ref.getM();
    const double charge = ref.getQ();
    Kokkos::parallel_for(
            "kickParticles", pc.getLocalNum(), KOKKOS_LAMBDA(const size_t i) {
                // Boris kick: Birdsall & Langdon (1985), ch. 4-4 (non-relativistic) and ch. 15-4
                // (relativistic); scaled P = (v/c)*gamma, mass in rest energy units.
                Vector_t<double, 3> p = Pview(i);
                // TODO: consider dropping unused R/dt arguments from kick when API allows.
                pusher.kick(0, p, Efview(i), Bfview(i), dtview(i), mass, charge);
                Pview(i) = p;
            });
    /*
    Wait until everyone completed the kick operation before proceeding. For now,
    this is just a precaution and could be removed for a small performance gain
    (technically, these shouldn't be necessary!).
    */
    Kokkos::fence();
    ippl::Comm->barrier();
    pc.markMomentsDirty();

    m << level5 << "Completed parallel kick operation." << endl;
}

// --- Helpers (beamline, dt, bounds, I/O) ---

/**
 * @copybrief ParallelTracker::prepareSections
 */
void ParallelTracker::prepareSections() {
    // Calls ParallelTracker::visitBeamline() -> TBeamline::iterate() ->
    // For each element:
    // FlaggedElemPtr::accept() -> DefaultVisitor::visitFlaggedElmPtr() ->
    // ElementBase::accept() -> ParallelTracker::visit[ElemName]() ->
    // OpalBeamline::visit([ElemName], bunch):
    // This initialises the FieldList elements_m object of OpalBeamline
    // with clones of all the elements
    itsBeamline_m.accept(*this);

    // Sorts the elements in OpalBeamline::elements_m by starting position
    itsOpalBeamline_m.prepareSections();

    // Computes the coordinate transformations for non-straight sections
    itsOpalBeamline_m.compute3DLattice();

    // Write 3D Lattice
    itsOpalBeamline_m.save3DLattice();
    itsOpalBeamline_m.save3DInput();
}

/**
 * @copybrief ParallelTracker::selectDT
 */
void ParallelTracker::selectDT() { itsBunch_m->setdT(dtCurrentTrack_m); }

/**
 * @copybrief ParallelTracker::changeDT
 */
void ParallelTracker::changeDT() {
    Inform m("ParallelTracker::changeDT");
    selectDT();
    double newdT = itsBunch_m->getdT();
    for (const auto& pc : itsBunch_m->getParticleContainers()) {
        if (pc) {
            pc->dt = newdT;
        }
    }
    m << level5 << "Changed particle container time step to " << newdT << "." << endl;
}

/**
 * @copybrief ParallelTracker::activateEmittingContainers
 */
void ParallelTracker::activateEmittingContainers(double t) {
    for (size_t i = 0; i < emittingSamplers_m.size(); ++i) {
        for (const auto& sampler : emittingSamplers_m[i]) {
            if (sampler && !sampler->isEmissionDone(t)) {
                itsBunch_m->setPcActive(i);
                break;
            }
        }
    }
}

/**
 * @copybrief ParallelTracker::doBinaryRepartition
 */
void ParallelTracker::doBinaryRepartition() {
    Inform m("ParallelTracker::doBinaryRepartition");
    if (itsBunch_m->hasFieldSolver()) {
        m << level4 << "*****************************************************************" << endl;
        m << level4 << "do repartition because of repartFreq_m" << endl;
        m << level4 << "*****************************************************************" << endl;
        itsBunch_m->do_binaryRepart();
        ippl::Comm->barrier();
        IpplTimings::stopTimer(BinRepartTimer_m);
        m << level4 << "*****************************************************************" << endl;
        m << level3 << "do repartition done" << endl;
        m << level4 << "*****************************************************************" << endl;
    }
}

/**
 * @copybrief ParallelTracker::computeInitialBounds
 */
void ParallelTracker::computeInitialBounds(Vector_t<double, 3>& rmin, Vector_t<double, 3>& rmax) {
    const auto& particleContainers = itsBunch_m->getParticleContainers();
    bool hasNonEmpty               = false;
    ippl::Vector<double, 3> rminLoc(0.0), rmaxLoc(0.0);

    for (const auto& pc : particleContainers) {
        if (!pc || pc->getTotalNum() == 0) {
            continue;
        }
        pc->computeMinMaxR();
        const ippl::Vector<double, 3> mn = pc->getMinR();
        const ippl::Vector<double, 3> mx = pc->getMaxR();
        if (!hasNonEmpty) {
            rminLoc     = mn;
            rmaxLoc     = mx;
            hasNonEmpty = true;
        } else {
            for (int i = 0; i < 3; ++i) {
                rminLoc[i] = std::min(rminLoc[i], mn[i]);
                rmaxLoc[i] = std::max(rmaxLoc[i], mx[i]);
            }
        }
    }

    if (!hasNonEmpty) {
        if (particleContainers.empty() || !particleContainers[0]) {
            throw OpalException(
                    "ParallelTracker::computeInitialBounds",
                    "No valid particle container for initial bounds.");
        }
        particleContainers[0]->computeMinMaxR();
        rminLoc = particleContainers[0]->getMinR();
        rmaxLoc = particleContainers[0]->getMaxR();
    }

    rmax = rmaxLoc;
    rmin = rminLoc;
    ippl::Comm->allreduce(rmax, 1, std::greater<ippl::Vector<double, 3>>());
    ippl::Comm->allreduce(rmin, 1, std::less<ippl::Vector<double, 3>>());
    ippl::Comm->barrier();
}

/**
 * @copybrief ParallelTracker::printInitialContainerRefs
 */
void ParallelTracker::printInitialContainerRefs(Inform& m) const {
    const auto& particleContainers = itsBunch_m->getParticleContainers();
    for (size_t i = 0; i < particleContainers.size(); ++i) {
        const auto& pc = particleContainers[i];
        if (!pc) {
            m << level3 << "ParallelTrack: container " << i << " is null." << endl;
            continue;
        }
        if (pc->getTotalNum() > 0) {
            m << level3 << "ParallelTrack (container " << i
              << "): momentum z = " << pc->getMeanP()(2) << "\n"
              << "RefPartR (container " << i << ") = " << pc->getRefPartR() << "\n"
              << "RefPartP (container " << i << ") = " << pc->getRefPartP() << endl;
        } else {
            m << level3 << "ParallelTrack: container " << i
              << " empty; total particles = " << itsBunch_m->getTotalNumAllContainers() << "\n"
              << "RefPartR (container " << i << ") = " << pc->getRefPartR() << "\n"
              << "RefPartP (container " << i << ") = " << pc->getRefPartP() << endl;
        }
    }
}

/**
 * @copybrief ParallelTracker::updateReference
 */
void ParallelTracker::updateReference(const BorisPusher& pusher) {
    Inform m("ParallelTracker::updateReference");
    updateReferenceParticles(pusher);
    updateRefToLabCSTrafo();
    m << level5 << "Updated reference particles." << endl;
}

/**
 * @copybrief ParallelTracker::updateReferenceParticles
 */
void ParallelTracker::updateReferenceParticles(const BorisPusher& pusher) {
    const double dt          = std::min(itsBunch_m->getT(), itsBunch_m->getdT());
    const double scaleFactor = Physics::c * dt;

    const size_t n = itsBunch_m->getNumParticleContainers();
    for (size_t i = 0; i < n; ++i) {
        if (!itsBunch_m->isPcActive(i)) {
            continue;
        }
        auto pcPtr = itsBunch_m->getParticleContainer(i);
        if (!pcPtr) {
            continue;
        }
        auto& pc                = *pcPtr;
        const PartData& refKick = *pc.getReference();
        Vector_t<double, 3> Ef(0.0), Bf(0.0);

        pc.getRefPartR() /= scaleFactor;
        pusher.push(pc.getRefPartR(), pc.getRefPartP(), dt);
        pc.getRefPartR() *= scaleFactor;

        IndexMap::value_t elements           = itsOpalBeamline_m.getElements(pc.getRefPartR());
        IndexMap::value_t::const_iterator it = elements.begin();
        const IndexMap::value_t::const_iterator end = elements.end();

        for (; it != end; ++it) {
            const CoordinateSystemTrafo& refToLocalCSTrafo =
                    itsOpalBeamline_m.getCSTrafoLab2Local((*it));

            Vector_t<double, 3> localR = refToLocalCSTrafo.transformTo(pc.getRefPartR());
            Vector_t<double, 3> localP = refToLocalCSTrafo.rotateTo(pc.getRefPartP());
            Vector_t<double, 3> localE(0.0), localB(0.0);

            if ((*it)->applyToReferenceParticle(
                        localR, localP, itsBunch_m->getT() - 0.5 * dt, localE, localB)) {
                *gmsg << level1 << "The reference particle hit an element" << endl;
                globalEOL_m = true;
            }

            Ef += refToLocalCSTrafo.rotateFrom(localE);
            Bf += refToLocalCSTrafo.rotateFrom(localB);
        }

        pusher.kick(pc.getRefPartR(), pc.getRefPartP(), Ef, Bf, dt, refKick.getM(), refKick.getQ());

        pc.getRefPartR() /= scaleFactor;
        pusher.push(pc.getRefPartR(), pc.getRefPartP(), dt);
        pc.getRefPartR() *= scaleFactor;
    }
}

/**
 * @copybrief ParallelTracker::updateRefToLabCSTrafo
 */
void ParallelTracker::updateRefToLabCSTrafo() {
    // Transform reference position to lab, but only rotate the momentum vector.
    // Momentum is a direction/axis and must not be translated.
    const double bunchDT = itsBunch_m->getdT();
    const size_t n       = itsBunch_m->getNumParticleContainers();
    for (size_t i = 0; i < n; ++i) {
        if (!itsBunch_m->isPcActive(i)) {
            continue;
        }
        auto pc = itsBunch_m->getParticleContainer(i);
        if (pc) {
            pc->updateRefToLabCSTrafo(bunchDT);
        }
    }
}

/**
 * @copybrief ParallelTracker::findStartPositions
 */
void ParallelTracker::findStartPositions(const BorisPusher& pusher) {
    // Primary container (index 0) drives segment advances.
    auto primary = itsBunch_m->getParticleContainer(0);

    if (zstart_m <= primary->get_sPos()) {
        return;
    }

    StepSizeConfig stepSizesCopy(stepSizes_m);

    double t = 0.0;
    itsBunch_m->setT(t);

    dtCurrentTrack_m = stepSizesCopy.getdT();
    selectDT();

    const auto& containers = itsBunch_m->getParticleContainers();

    while (true) {
        autophaseCavities(pusher);

        t += itsBunch_m->getdT();
        itsBunch_m->setT(t);

        std::vector<Vector_t<double, 3>> oldRs(containers.size());
        for (size_t i = 0; i < containers.size(); ++i) {
            if (containers[i]) {
                oldRs[i] = containers[i]->getRefPartR();
            }
        }

        updateReferenceParticles(pusher);

        for (size_t i = 0; i < containers.size(); ++i) {
            if (!containers[i]) {
                continue;
            }
            Vector_t<double, 3> dR = containers[i]->getRefPartR() - oldRs[i];
            containers[i]->set_sPos(containers[i]->get_sPos() + euclidean_norm(dR));
        }

        Vector_t<double, 3> tmp =
                primary->getRefPartP() * Physics::c / Util::getGamma(primary->getRefPartP());
        double speed = euclidean_norm(tmp);

        if (primary->get_sPos() > stepSizesCopy.getZStop()) {
            ++stepSizesCopy;

            if (stepSizesCopy.reachedEnd()) {
                --stepSizesCopy;
                const double zTarget    = stepSizesCopy.getZStop();
                const double tauPrimary = (zTarget - primary->get_sPos()) / speed;
                itsBunch_m->setT(itsBunch_m->getT() + tauPrimary);
                for (const auto& pc : containers) {
                    if (!pc) {
                        continue;
                    }
                    Vector_t<double, 3> pv =
                            pc->getRefPartP() * Physics::c / Util::getGamma(pc->getRefPartP());
                    double speed_i = euclidean_norm(pv);
                    double tau_i   = (zTarget - pc->get_sPos()) / speed_i;
                    pc->applyFractionalStep(pusher, tau_i, zstart_m);
                }

                break;
            }

            dtCurrentTrack_m = stepSizesCopy.getdT();
            selectDT();
        }

        if (std::abs(primary->get_sPos() - zstart_m) <= 0.5 * itsBunch_m->getdT() * speed) {
            const double zTarget    = zstart_m;
            const double tauPrimary = (zTarget - primary->get_sPos()) / speed;
            itsBunch_m->setT(itsBunch_m->getT() + tauPrimary);
            for (const auto& pc : containers) {
                if (!pc) {
                    continue;
                }
                Vector_t<double, 3> pv =
                        pc->getRefPartP() * Physics::c / Util::getGamma(pc->getRefPartP());
                double speed_i = euclidean_norm(pv);
                double tau_i   = (zTarget - pc->get_sPos()) / speed_i;
                pc->applyFractionalStep(pusher, tau_i, zstart_m);
            }

            break;
        }
    }

    changeDT();
}

/**
 * @copybrief ParallelTracker::dumpStats
 */
void ParallelTracker::dumpStats(long long step, bool psDump, bool statDump) {
    OPALTimer::Timer myt2;
    const size_t totalAll = itsBunch_m->getTotalNumAllContainers();

    if (totalAll == 0) {
        *gmsg << level1 << "* " << myt2.time() << " "
              << "Step " << std::setw(6) << itsBunch_m->getGlobalTrackStep() << "; "
              << "   -- no emission yet --     "
              << "t= " << Util::getTimeString(itsBunch_m->getT()) << endl;
        return;
    }

    bool anyLogged         = false;
    const auto& containers = itsBunch_m->getParticleContainers();
    for (size_t ci = 0; ci < containers.size(); ++ci) {
        const auto& pc = containers[ci];
        if (!pc || pc->getTotalNum() == 0) {
            continue;
        }
        pc->updateMoments();
        const double sPos = pc->get_sPos();
        if (std::isnan(sPos) || std::isinf(sPos)) {
            throw OpalException(
                    "ParallelTracker::dumpStats()",
                    "invalid path length s for particle container " + std::to_string(ci));
        }
        *gmsg << level1 << "* " << myt2.time() << " "
              << "Step " << std::setw(6) << itsBunch_m->getGlobalTrackStep() << " "
              << "container[" << ci << "] "
              << "at " << Util::getLengthString(sPos) << ", "
              << "t= " << Util::getTimeString(itsBunch_m->getT()) << ", "
              << "E=" << Util::getEnergyString(pc->getMeanKineticEnergy()) << endl;
        anyLogged = true;
    }

    if (anyLogged) {
        writePhaseSpace(step, psDump, statDump);
    }
}

/**
 * @copybrief ParallelTracker::setOptionalVariables
 */
void ParallelTracker::setOptionalVariables() {
    /*
    minStepforReBin_m = Options::minStepForRebin;
    RealVariable* br =
        dynamic_cast<RealVariable*>(OpalData::getInstance()->find("MINSTEPFORREBIN"));
    if (br)
        minStepforReBin_m = static_cast<int>(br->getReal());
    msg << level2 << "MINSTEPFORREBIN " << minStepforReBin_m << endl;
    */
    Inform m("ParallelTracker::setOptionalVariables");

    // there is no point to do repartitioning with one node
    if (ippl::Comm->size() == 1) {
        repartFreq_m = std::numeric_limits<unsigned long long>::max();
    } else {
        repartFreq_m = Options::repartFreq * 100;
        RealVariable* rep =
                dynamic_cast<RealVariable*>(OpalData::getInstance()->find("REPARTFREQ"));
        if (rep) repartFreq_m = static_cast<int>(rep->getReal());
        m << level3 << "REPARTFREQ set to " << repartFreq_m << "." << endl;
    }
}

/**
 * @copybrief ParallelTracker::hasEndOfLineReached
 */
bool ParallelTracker::hasEndOfLineReached(const BoundingBox& globalBoundingBox) {
    // Old IPPL used OpBitwiseAndAssign() via reduce(); new IPPL needs allreduce
    // so that all ranks receive the result (reduce sends only to root).
    // In-place allreduce avoids the aliased-buffer error of MPI_Reduce.
    ippl::Comm->allreduce(globalEOL_m, 1, std::logical_and<bool>());
    globalEOL_m = globalEOL_m
                  || globalBoundingBox.isOutside(itsBunch_m->getParticleContainer()->getRefPartR());
    return globalEOL_m;
}

/**
 * @copybrief ParallelTracker::setTime
 */
void ParallelTracker::setTime() {
    double newdT = itsBunch_m->getdT();
    for (const auto& pc : itsBunch_m->getParticleContainers()) {
        if (pc) {
            pc->dt = newdT;
        }
    }
}

/**
 * @copybrief ParallelTracker::writePhaseSpace
 */
void ParallelTracker::writePhaseSpace(const long long /*step*/, bool psDump, bool statDump) {
    Inform m("ParallelTracker::writePhaseSpace");
    Vector_t<double, 3> externalE, externalB;

    Vector_t<double, 3> rmin, rmax;
    itsBunch_m->get_bounds(rmin, rmax);
    m << level5 << "Bunch bounds in REFERENCE frame: rmin = " << rmin << ", rmax = " << rmax
      << endl;

    const size_t nContainers = itsBunch_m->getNumParticleContainers();
    std::vector<std::array<Vector_t<double, 3>, 2>> fdByContainer(nContainers);

    if (psDump || statDump) {
        for (size_t i = 0; i < nContainers; ++i) {
            auto pc = itsBunch_m->getParticleContainer(i);
            if (!pc || pc->getTotalNum() == 0) {
                fdByContainer[i][0] = Vector_t<double, 3>(0.0);
                fdByContainer[i][1] = Vector_t<double, 3>(0.0);
                continue;
            }
            externalB = Vector_t<double, 3>(0.0);
            externalE = Vector_t<double, 3>(0.0);
            itsOpalBeamline_m.getFieldAt(
                    pc->getRefPartR(), pc->getRefPartP(),
                    itsBunch_m->getT() - 0.5 * itsBunch_m->getdT(), externalE, externalB);
            fdByContainer[i][0] = externalB;
            fdByContainer[i][1] = externalE * Units::Vpm2MVpm;
            m << level5 << "External fields (container " << i << "): externalE = " << externalE
              << ", externalB = " << externalB << endl;
        }
    }

    if (statDump) {
        itsDataSink_m->dumpSDDS(*itsBunch_m, fdByContainer, -1.0);
        *gmsg << level2 << "* Wrote beam statistics." << endl;
    }

    if (psDump && itsBunch_m->getTotalNumAllContainers() > 0) {
        itsDataSink_m->dumpH5(*itsBunch_m, fdByContainer);

        /*
        // Write fields to .h5 file.
        const size_t localNum    = itsBunch_m->getLocalNum();
        double distToLastStop    = stepSizes_m.getFinalZStop() -
        itsBunch_m->getParticleContainer()->get_sPos(); Vector_t<double, 3> beta =
        itsBunch_m->RefPartP_m / Util::getGamma(itsBunch_m->RefPartP_m); Vector_t<double, 3>
        driftPerTimeStep = itsBunch_m->getdT()
            * Physics::c;  // \todo  * itsBunch_m->toLabTrafo_m.rotateFrom(beta);
        bool driftToCorrectPosition =
            std::abs(distToLastStop) < 0.5 * euclidean_norm(driftPerTimeStep);
        // \todo Ppos_t stashedR;
        Vector_t<double, 3> stashedR;
        Vector_t<double, 3> stashedRefPartR;

        if (driftToCorrectPosition) {
            const double tau =
                distToLastStop / euclidean_norm(driftPerTimeStep) * itsBunch_m->getdT();

            if (localNum > 0) {

                stashedR.create(localNum);
                stashedR        = itsBunch_m->R;
                stashedRefPartR = itsBunch_m->RefPartR_m;

                for (size_t i = 0; i < localNum; ++i) {
                    itsBunch_m->R[i] +=
                        tau
                        * (Physics::c * itsBunch_m->P[i] / Util::getGamma(itsBunch_m->P[i])
                           - driftPerTimeStep / itsBunch_m->getdT());
                }

            }

            driftPerTimeStep = itsBunch_m->toLabTrafo_m.rotateTo(driftPerTimeStep);
            itsBunch_m->RefPartR_m =
                itsBunch_m->RefPartR_m + tau * driftPerTimeStep / itsBunch_m->getdT();
            CoordinateSystemTrafo update(
                tau * driftPerTimeStep / itsBunch_m->getdT(), Quaternion(1.0, 0.0, 0.0, 0.0));
            itsBunch_m->toLabTrafo_m = itsBunch_m->toLabTrafo_m * update.inverted();

            itsBunch_m->set_sPos(stepSizes_m.getFinalZStop());

            itsBunch_m->calcBeamParameters();
        }
        if (!statDump && !driftToCorrectPosition)
            itsBunch_m->calcBeamParameters();

            msg << *itsBunch_m << endl;

        itsDataSink_m->dumpH5(*itsBunch_m, FDext);
        */

        /*

        if (driftToCorrectPosition) {
            if (localNum > 0) {
                itsBunch_m->R = stashedR;
                stashedR.destroy(localNum, 0);
            }

            itsBunch_m->RefPartR_m = stashedRefPartR;
            itsBunch_m->getParticleContainer()->set_sPos(
                itsBunch_m->getParticleContainer()->get_sPos());

            itsBunch_m->calcBeamParameters();
        }
        */

        *gmsg << level3 << "* Wrote beam phase space." << endl;
    }
}

// --- Autophasing (RF and traveling-wave cavities) ---

/**
 * @copybrief ParallelTracker::updateRFElement
 */
void ParallelTracker::updateRFElement(std::string elName, double maxPhase) {
    Inform m("ParallelTracker::updateRFElement");
    FieldList cavities       = itsOpalBeamline_m.getElementByType(ElementType::RFCAVITY);
    FieldList travelingwaves = itsOpalBeamline_m.getElementByType(ElementType::TRAVELINGWAVE);
    cavities.insert(cavities.end(), travelingwaves.begin(), travelingwaves.end());
    m << level5 << "Got cavities and traveling waves." << endl;

    for (FieldList::iterator fit = cavities.begin(); fit != cavities.end(); ++fit) {
        if ((*fit).getElement()->getName() == elName) {
            RFCavity* element = static_cast<RFCavity*>((*fit).getElement().get());

            element->setPhasem(maxPhase);
            element->setAutophaseVeto();

            m << level3 << "Restored cavity phase from the h5 file. Name: " << element->getName()
              << ", phase: " << maxPhase << " rad" << endl;
            return;
        }
    }
}

/**
 * @copybrief ParallelTracker::saveCavityPhases
 */
void ParallelTracker::saveCavityPhases() { itsDataSink_m->storeCavityInformation(); }

/**
 * @copybrief ParallelTracker::restoreCavityPhases
 */
void ParallelTracker::restoreCavityPhases() {
    typedef std::vector<MaxPhasesT>::iterator iterator_t;

    if (OpalData::getInstance()->hasPriorTrack() || OpalData::getInstance()->inRestartRun()) {
        iterator_t it  = OpalData::getInstance()->getFirstMaxPhases();
        iterator_t end = OpalData::getInstance()->getLastMaxPhases();
        for (; it < end; ++it) {
            updateRFElement((*it).first, (*it).second);
        }
    }
}

/**
 * @copybrief ParallelTracker::autophaseCavities
 */
void ParallelTracker::autophaseCavities(const BorisPusher& pusher) {
    const PartData& ref = *itsBunch_m->getParticleContainer()->getReference();
    double t            = itsBunch_m->getT();
    Vector_t<double, 3> nextR =
            itsBunch_m->getParticleContainer()->getRefPartR() / (Physics::c * itsBunch_m->getdT());
    pusher.push(nextR, itsBunch_m->getParticleContainer()->getRefPartP(), itsBunch_m->getdT());
    nextR *= Physics::c * itsBunch_m->getdT();

    auto elementSet = itsOpalBeamline_m.getElements(nextR);
    for (auto element : elementSet) {
        if (element->getType() == ElementType::TRAVELINGWAVE) {
            const TravelingWave* TWelement = static_cast<const TravelingWave*>(element.get());
            if (!TWelement->getAutophaseVeto()) {
                CavityAutophaser ap(ref, element);
                ap.getPhaseAtMaxEnergy(
                        itsOpalBeamline_m.transformToLocalCS(
                                element, itsBunch_m->getParticleContainer()->getRefPartR()),
                        itsOpalBeamline_m.rotateToLocalCS(
                                element, itsBunch_m->getParticleContainer()->getRefPartP()),
                        t, itsBunch_m->getdT());
            }

        } else if (element->getType() == ElementType::RFCAVITY) {
            const RFCavity* RFelement = static_cast<const RFCavity*>(element.get());
            if (!RFelement->getAutophaseVeto()) {
                CavityAutophaser ap(ref, element);
                ap.getPhaseAtMaxEnergy(
                        itsOpalBeamline_m.transformToLocalCS(
                                element, itsBunch_m->getParticleContainer()->getRefPartR()),
                        itsOpalBeamline_m.rotateToLocalCS(
                                element, itsBunch_m->getParticleContainer()->getRefPartP()),
                        t, itsBunch_m->getdT());
            }
        }
    }
}

// --- Commented-out legacy / not-in-use implementations (reference only) ---
/*
void ParallelTracker::visitScalingFFAMagnet(const ScalingFFAMagnet& //bend) {
    *gmsg << level4 << "Adding ScalingFFAMagnet" << endl;
    *gmsg << level4 << "passed ScalingFFAMagnet argument not used in
ParallelTracker::visitScalingFFAMagnet" << endl;
}

void ParallelTracker::visitRing(const Ring& ring) {
    *gmsg << level4 << "* ----------------------------- Ring -------------------------------------
*" << endl;

    delete opalRing_m;
    opalRing_m = dynamic_cast<Ring*>(ring.clone());
    myElements.push_back(opalRing_m);
    opalRing_m->initialise(itsBunch_m);

    referenceR     = opalRing_m->getBeamRInit();
    referencePr    = opalRing_m->getBeamPRInit();
    referenceTheta = opalRing_m->getBeamPhiInit();

    if (referenceTheta <= -180.0 || referenceTheta > 180.0) {
        throw OpalException(
            "Error in ParallelTracker::visitRing", "PHIINIT is out of [-180, 180)!");
    }

    referenceZ  = 0.0;
    referencePz = 0.0;

    referencePtot = itsReference.getGamma() * itsReference.getBeta();
    referencePt   = std::sqrt(referencePtot * referencePtot - referencePr * referencePr);

    if (referencePtot < 0.0)
        referencePt *= -1.0;

    sinRefTheta_m = std::sin(referenceTheta * Units::deg2rad);
    cosRefTheta_m = std::cos(referenceTheta * Units::deg2rad);

    double BcParameter[8] = {};
    buildupFieldList(BcParameter, ElementType::RING, opalRing_m);
}

void ParallelTracker::visitVerticalFFAMagnet(const VerticalFFAMagnet& mag) {
    *gmsg << level4 << "Adding Vertical FFA Magnet" << endl;
    if (opalRing_m != nullptr)
        opalRing_m->appendElement(mag);
    else
        throw OpalException(
            "ParallelCyclotronTracker::visitVerticalFFAMagnet",
            "Need to define a RINGDEFINITION to use VerticalFFAMagnet element");
}

void ParallelTracker::buildupFieldList(
    double BcParameter[], ElementType elementType, Component* elptr) {
    beamline_list::iterator sindex;
    type_pair* localpair = new type_pair();
    localpair->first     = elementType;

    for (int i = 0; i < 8; i++)
        *(((localpair->second).first) + i) = *(BcParameter + i);

    (localpair->second).second = elptr;
    if (elementType == ElementType::RING) {
        sindex = FieldDimensions.begin();
    } else {
        sindex = FieldDimensions.end();
    }
    FieldDimensions.insert(sindex, localpair);
}

bool ParallelTracker::applyPluginElements(const double dt) {
    IpplTimings::startTimer(PluginElemTimer_m);

    bool flag = false;
    for (PluginElement* element : pluginElements_m) {
        bool tmp = element->check(itsBunch_m, turnnumber_m, itsBunch_m->getT(), dt);
        flag |= tmp;
    }

    IpplTimings::stopTimer(PluginElemTimer_m);
    return flag;
}
*/
